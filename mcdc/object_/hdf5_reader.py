"""Read mcdc-data representations into transport objects without unit conversion.

The HDF5 format uses eV and zero-based inclusive interpolation boundaries.
DataTable stores exclusive boundaries, so this module performs that index
conversion at the file boundary. Source-library metadata stays in the files.
"""

import numpy as np

from mcdc.constant import INTERPOLATION_LINEAR, REFERENCE_FRAME_COM, REFERENCE_FRAME_LAB

from mcdc.object_.data import DataPolynomial, DataTable, encode_interpolation
from mcdc.object_.distribution import (
    DistributionNone,
    DistributionTabulated,
    DistributionMultiTable,
    DistributionLevelScattering,
    DistributionEvaporation,
    DistributionMaxwellian,
    DistributionKalbachMann,
    DistributionTabulatedEnergyAngle,
    DistributionNBody,
)

# ======================================================================================
# File identity and interpolation
# ======================================================================================


def text(value):
    """Decode an HDF5 string regardless of the h5py string representation."""
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def check_format(file):
    """Reject files whose format identity or schema version is unsupported."""
    name = text(file.attrs.get("format", ""))
    version = text(file.attrs.get("format-version", ""))
    if name != "mcdc-data" or version != "0.1.0":
        raise ValueError(
            f"{file.filename}: expected mcdc-data 0.1.0, found {name!r} {version!r}"
        )


def read_interpolation(group):
    """Return packed laws and exclusive boundaries from an interpolation group."""
    interpolation = group["interpolation"]
    laws = np.array(
        [encode_interpolation(text(law)) for law in interpolation["laws"][()]],
        dtype=np.int64,
    )
    boundaries = np.asarray(interpolation["boundaries"][()], dtype=np.int64) + 1
    coordinate = next(name for name in ("value", "grid", "energy") if name in group)
    size = len(group[coordinate])
    if (
        laws.ndim != 1
        or boundaries.ndim != 1
        or len(laws) != len(boundaries)
        or not len(laws)
    ):
        raise ValueError(
            f"{group.name}: interpolation laws and boundaries must be nonempty vectors of equal length"
        )
    if boundaries[0] < 2 or boundaries[-1] != size or np.any(np.diff(boundaries) <= 0):
        raise ValueError(f"{group.name}: invalid interpolation region boundaries")
    return laws, boundaries


def numbered_groups(group, prefix=""):
    """Read contiguous zero-based child indices in numeric order."""
    names = [name for name in group if name.startswith(prefix)]
    names.sort(key=lambda name: int(name[len(prefix) :]))
    if names != [f"{prefix}{i}" for i in range(len(names))]:
        raise ValueError(
            f"{group.name}: expected contiguous zero-based {prefix!r} indices"
        )
    return [group[name] for name in names]


# ======================================================================================
# Shared-grid cross sections and reaction fields
# ======================================================================================


def read_xs_energy_grid(group):
    """Read a projectile's eV grid, retaining repeated coordinates at jumps."""
    grid = np.asarray(group["xs-energy-grid"][()], dtype=np.float64)
    if (
        grid.ndim != 1
        or len(grid) < 2
        or not np.all(np.isfinite(grid))
        or np.any(np.diff(grid) < 0)
        or grid[-1] <= grid[0]
    ):
        raise ValueError(
            f"{group.name}: cross-section energies must be finite and nondecreasing with at least two distinct values"
        )
    return grid


def read_reaction_properties(group):
    """Read source MT, cross-section segment, offset, and reference frame.

    Reaction category determines transport dispatch; ENDF-MT is optional source
    metadata. Cross-product data must not be silently discarded by a reader.
    """
    if len(group.get("cross-products", ())):
        raise ValueError(f"{group.name}: cross-product transport is unsupported")
    MT = int(group.attrs.get("ENDF-MT", 0))
    xs = np.asarray(group["xs/value"][()], dtype=np.float64)
    offset = int(group["xs/offset"][()])
    frame = text(group["reference-frame"][()])
    frames = {"LAB": REFERENCE_FRAME_LAB, "COM": REFERENCE_FRAME_COM}
    if frame not in frames:
        raise ValueError(f"{group.name}: unknown reference frame {frame!r}")
    return MT, xs, offset, frames[frame]


def set_xs_interpolation(reaction, group):
    """Validate shared-grid cross sections and attach their interpolation metadata.

    Neutron and photon reactions retain their own interpolation regions while
    sharing the projectile energy grid. File boundaries are converted to the
    exclusive point indices used by the transport objects.
    """
    laws, boundaries = read_interpolation(group["xs"])
    if (
        reaction.xs.ndim != 1
        or reaction.xs_offset_ < 0
        or np.any(reaction.xs < 0)
        or not np.all(np.isfinite(reaction.xs))
    ):
        raise ValueError(f"{group.name}: cross sections must be finite and nonnegative")
    reaction.xs_interpolations = laws
    reaction.xs_interpolation_boundaries = boundaries
    return reaction


def xs_is_linear(reaction):
    """Whether summing grid values before interpolation preserves this reaction.

    A nonzero first value at a positive offset opens abruptly at its threshold.
    Such a segment must be evaluated separately to avoid interpolating the
    summed array into the interval below that threshold.
    """
    return bool(
        np.all(reaction.xs_interpolations == INTERPOLATION_LINEAR)
        and (reaction.xs_offset_ == 0 or reaction.xs[0] == 0.0)
    )


# ======================================================================================
# Scalar data
# ======================================================================================


def read_data(group):
    """Load constant, polynomial, or tabulated scalar data in canonical units."""
    kind = text(group.attrs["type"])
    if kind == "constant":
        return DataPolynomial(np.array([float(group["value"][()])]))
    if kind == "polynomial":
        return DataPolynomial(group["coefficient"][()])
    if kind == "tabulated":
        return DataTable(
            group["grid"][()], group["value"][()], *read_interpolation(group)
        )
    raise ValueError(f"{group.name}: unsupported Data type {kind!r}")


# ======================================================================================
# Probability distributions
# ======================================================================================


def read_distribution(group):
    """Load an angle, energy, or correlated energy-angle distribution.

    Nested tables retain their own density interpolation. Incident-energy
    interpolation is carried separately, so angular sampling can interpolate
    probabilities while outgoing energies additionally use unit-base scaling.
    """
    kind = text(group.attrs["type"])
    if kind == "isotropic":
        return DistributionNone()
    if kind == "tabulated":
        laws, boundaries = read_interpolation(group)
        return DistributionTabulated(
            group["value"][()],
            pdf=group["pdf"][()],
            interpolations=laws,
            interpolation_boundaries=boundaries,
        )
    if kind == "incident-energy-dependent":
        laws, boundaries = read_interpolation(group)
        return DistributionMultiTable(
            group["energy"][()],
            tables=[
                read_distribution(child)
                for child in numbered_groups(group["distributions"])
            ],
            interpolations=laws,
            interpolation_boundaries=boundaries,
        )
    if kind == "level-scattering":
        return DistributionLevelScattering(
            float(group["C1"][()]), float(group["C2"][()])
        )
    if kind in ("evaporation", "maxwellian"):
        cls = (
            DistributionEvaporation if kind == "evaporation" else DistributionMaxwellian
        )
        return cls(
            read_data(group["temperature"]),
            restriction_energy=float(group["restriction-energy"][()]),
        )
    if kind == "n-body":
        density = read_distribution_density(group)
        result = DistributionNBody(density.x, density.y)
        result.pdf = density
        return result
    if kind in ("kalbach-mann", "energy-angle-tabulated"):
        return read_correlated_distribution(group, kind)
    raise ValueError(f"{group.name}: unsupported Distribution type {kind!r}")


def read_distribution_density(group):
    """Return a normalized density table including its cumulative probabilities."""
    laws, boundaries = read_interpolation(group)
    return DistributionTabulated(
        group["value"][()],
        pdf=group["pdf"][()],
        interpolations=laws,
        interpolation_boundaries=boundaries,
    ).pdf


def segment_laws(table):
    """Expand interpolation regions to one code per point for flattened tables."""
    result = np.empty(table.N, dtype=np.int64)
    start = 0
    for law, boundary in zip(table.interpolations, table.interpolation_boundaries):
        result[start:boundary] = law
        start = boundary - 1
    return result


def read_correlated_distribution(group, kind):
    """Flatten conditional tables while preserving local indices and interpolation."""
    children = numbered_groups(group["tables"])
    tables = [read_distribution_density(child["outgoing-energy"]) for child in children]
    offset = np.cumsum([0] + [table.N for table in tables[:-1]], dtype=np.int64)
    values = np.concatenate([table.x for table in tables])
    pdf = np.concatenate([table.y for table in tables])
    laws = np.concatenate([segment_laws(table) for table in tables])
    if len(children) != len(group["energy"]):
        raise ValueError(
            f"{group.name}: incident energies and tables must have equal lengths"
        )
    if kind == "kalbach-mann":
        factors = [read_data(child["precompound-factor"]) for child in children]
        slopes = [read_data(child["angular-slope"]) for child in children]
        # Parameters retain their evaluable Data representations at each incident energy.
        result = DistributionKalbachMann(
            group["energy"][()],
            offset,
            values,
            pdf,
            np.zeros(len(values)),
            np.zeros(len(values)),
        )
        result.precompound_tables = factors
        result.angular_slope_tables = slopes
    else:
        angles = []
        for child, table in zip(children, tables):
            row = numbered_groups(child["cosine-distributions"])
            if len(row) != table.N:
                raise ValueError(
                    f"{child.name}: each outgoing energy requires a cosine distribution"
                )
            angles.extend(read_distribution_density(item) for item in row)
        cosine_offset = np.cumsum(
            [0] + [table.N for table in angles[:-1]], dtype=np.int64
        )
        result = DistributionTabulatedEnergyAngle(
            group["energy"][()],
            offset,
            values,
            pdf,
            cosine_offset,
            np.concatenate([table.x for table in angles]),
            np.concatenate([table.y for table in angles]),
        )
        result.cosine_pdf = np.concatenate([table.y for table in angles])
        result.cosine_cdf = np.concatenate([table.aux[0] for table in angles])
        result.cosine_interpolations = np.concatenate(
            [segment_laws(table) for table in angles]
        )
    result.pdf = pdf
    result.cdf = np.concatenate([table.aux[0] for table in tables])
    result.interpolations = laws
    return result
