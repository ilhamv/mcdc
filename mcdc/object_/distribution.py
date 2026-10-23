import numpy as np

from numbers import Integral
from numpy import float64, int64
from numpy.typing import ArrayLike, NDArray

from mcdc.constant import (
    DISTRIBUTION_NONE,
    DISTRIBUTION_TABLE,
    DISTRIBUTION_MULTITABLE,
    DISTRIBUTION_LEVEL_SCATTERING,
    DISTRIBUTION_EVAPORATION,
    DISTRIBUTION_MAXWELLIAN,
    DISTRIBUTION_KALBACH_MANN,
    DISTRIBUTION_TABLE_ENERGY_ANGLE,
    DISTRIBUTION_N_BODY,
    INTERPOLATION_HISTOGRAM,
    INTERPOLATION_LINEAR,
)
from mcdc.object_.base import MCDCPolymorphic
from mcdc.object_.data import DataBase, DataNone, DataTable, decode_interpolation
from mcdc.object_.util import (
    cdf_from_pdf,
    pdf_from_cdf,
    multi_cdf_from_pdf,
)
from mcdc.print_ import print_1d_array, print_error

# ======================================================================================
# Distribution base class
# ======================================================================================


class DistributionBase(MCDCPolymorphic):
    """Base class for probability distributions sampled during transport."""

    # MC/DC framework metadata
    label = "distribution"
    sub_type = -1  # Polymorphic base


# ======================================================================================
# None
# ======================================================================================
# Placeholder for a distribution that does not need to store data:
#   - Isotropic distributions
#   - Energy-correlated angles stored in the energy distribution


class DistributionNone(DistributionBase):
    """Placeholder for an implicit or externally stored distribution."""

    # MC/DC framework metadata
    label = "none_distribution"
    sub_type = DISTRIBUTION_NONE


# ======================================================================================
# Table
# ======================================================================================


class DistributionTable(DistributionBase):
    """One-dimensional continuous, discrete, or mixed distribution.

    ``value`` with ``pdf`` (or a continuous ``cdf``) describes a density;
    ``discrete_value`` with ``pmf`` describes point masses. Either component
    may be omitted. The PDF integral and PMF sum are normalized together,
    preserving their relative probabilities. Interpolation applies only to
    the continuous density. CDF input describes a piecewise-linear CDF.
    """

    label = "table_distribution"
    sub_type = DISTRIBUTION_TABLE

    pdf: DataBase
    discrete_value: NDArray[float64]
    pmf: NDArray[float64]
    cmf: NDArray[float64]
    continuous_probability: float

    def __init__(
        self,
        value=None,
        pdf=None,
        cdf=None,
        *,
        discrete_value=None,
        pmf=None,
        interpolations=INTERPOLATION_LINEAR,
        interpolation_boundaries=None,
    ):
        super().__init__()
        continuous = value is not None
        discrete = discrete_value is not None
        if continuous != (pdf is not None or cdf is not None):
            raise ValueError("Continuous values require exactly one PDF or CDF")
        if pdf is not None and cdf is not None:
            raise ValueError("Supply either PDF or CDF, not both")
        if discrete != (pmf is not None) or not (continuous or discrete):
            raise ValueError(
                "Supply a continuous density, discrete values with PMF, or both"
            )
        self.discrete_value = np.asarray(
            [] if discrete_value is None else discrete_value, dtype=float64
        )
        masses = np.asarray([] if pmf is None else pmf, dtype=float64)
        if (
            self.discrete_value.ndim != 1
            or masses.ndim != 1
            or self.discrete_value.shape != masses.shape
            or (discrete and len(masses) == 0)
            or not np.all(np.isfinite(self.discrete_value))
            or not np.all(np.isfinite(masses))
            or np.any(masses < 0)
        ):
            raise ValueError(
                "Discrete values and nonnegative finite masses must be matching vectors"
            )
        self.pdf = DataNone()
        continuous_mass = 0.0
        if continuous:
            x = np.asarray(value, dtype=float64)
            if cdf is not None:
                raw_cdf = np.asarray(cdf, dtype=float64)
                density, cumulative = pdf_from_cdf(x, raw_cdf)
                continuous_mass = float(raw_cdf[-1] - raw_cdf[0])
                density, cumulative = (
                    density * continuous_mass,
                    cumulative * continuous_mass,
                )
                interpolations, interpolation_boundaries = INTERPOLATION_HISTOGRAM, None
            else:
                density, cumulative = cdf_from_pdf(
                    x, pdf, interpolations, interpolation_boundaries, normalize=False
                )
                continuous_mass = float(cumulative[-1])
        total = float(np.sum(masses)) + continuous_mass
        if not np.isfinite(total) or total <= 0:
            raise ValueError("Combined probability must have positive finite mass")
        self.pmf = masses / total
        self.cmf = np.r_[0.0, np.cumsum(self.pmf)]
        self.continuous_probability = continuous_mass / total
        if continuous_mass > 0:
            self.pdf = DataTable(
                x,
                density / total,
                interpolations,
                interpolation_boundaries,
                aux=cumulative / total,
            )

    def __repr__(self):
        text = super().__repr__()
        text += f"  - continuous probability: {self.continuous_probability}\n"
        text += f"  - discrete values: {print_1d_array(self.discrete_value)}\n"
        text += f"  - probability masses: {print_1d_array(self.pmf)}\n"
        if self.continuous_probability > 0:
            text += f"  - continuous values: {print_1d_array(self.pdf.x)}\n"
            text += f"  - probability density: {print_1d_array(self.pdf.y)}\n"
        return text


# ======================================================================================
# Multi-table
# ======================================================================================


class DistributionMultiTable(DistributionBase):
    """Conditional continuous PDFs stored in independent packed arrays.

    ``value`` and ``pdf`` (or ``cdf``) may be sequences of numerical vectors,
    one per indexing-grid point. Alternatively supply concatenated vectors
    with N + 1 ``offsets``. No DistributionTable objects are constructed.
    Outer interpolation is histogram or linear. Packed PDF interpolation uses
    local exclusive boundaries and N + 1 ``pdf_interpolation_offsets``;
    a scalar PDF law applies to every table. Every PDF is normalized separately.
    """

    label = "multi_table_distribution"
    sub_type = DISTRIBUTION_MULTITABLE

    grid: NDArray[float64]
    offsets: NDArray[int64]
    value: NDArray[float64]
    pdf: NDArray[float64]
    cdf: NDArray[float64]
    interpolations: NDArray[int64]
    interpolation_boundaries: NDArray[int64]
    pdf_interpolations: NDArray[int64]
    pdf_interpolation_boundaries: NDArray[int64]
    pdf_interpolation_offsets: NDArray[int64]

    def __init__(
        self,
        grid,
        value,
        pdf=None,
        *,
        cdf=None,
        offsets=None,
        interpolations=INTERPOLATION_LINEAR,
        interpolation_boundaries=None,
        pdf_interpolations=INTERPOLATION_LINEAR,
        pdf_interpolation_boundaries=None,
        pdf_interpolation_offsets=None,
    ):
        super().__init__()
        self.grid = np.asarray(grid, dtype=float64)
        if (
            self.grid.ndim != 1
            or not len(self.grid)
            or not np.all(np.isfinite(self.grid))
            or np.any(np.diff(self.grid) < 0)
        ):
            raise ValueError(
                "Multi-table grid must be finite, nonempty, and nondecreasing"
            )
        count = len(self.grid)
        self.interpolations, self.interpolation_boundaries = _distribution_regions(
            count, interpolations, interpolation_boundaries, outer=True
        )
        if offsets is not None:
            raw_offsets = np.asarray(offsets)
            if (
                raw_offsets.ndim != 1
                or not np.issubdtype(raw_offsets.dtype, np.integer)
                or len(raw_offsets) != count + 1
                or raw_offsets[0] != 0
                or np.any(np.diff(raw_offsets) < 2)
                or raw_offsets[-1] != len(value)
            ):
                raise ValueError(
                    "Packed PDF offsets must delimit one table per grid point"
                )
            slices = [slice(a, b) for a, b in zip(raw_offsets[:-1], raw_offsets[1:])]
            values = [np.asarray(value)[s] for s in slices]
            pdfs = [None if pdf is None else np.asarray(pdf)[s] for s in slices]
            cdfs = [None if cdf is None else np.asarray(cdf)[s] for s in slices]
            if (pdf is not None and len(pdf) != len(value)) or (
                cdf is not None and len(cdf) != len(value)
            ):
                raise ValueError("Packed density arrays must have equal lengths")
        else:
            values = list(value)
            pdfs = [None] * count if pdf is None else list(pdf)
            cdfs = [None] * count if cdf is None else list(cdf)
        if len(values) != count or len(pdfs) != count or len(cdfs) != count:
            raise ValueError("Grid and PDF table counts must match")
        if not isinstance(pdf_interpolations, Integral):
            raw_rules = np.asarray(pdf_interpolation_offsets)
            if (
                raw_rules.ndim != 1
                or not np.issubdtype(raw_rules.dtype, np.integer)
                or len(raw_rules) != count + 1
                or raw_rules[0] != 0
                or np.any(np.diff(raw_rules) < 1)
                or raw_rules[-1] != len(pdf_interpolations)
                or pdf_interpolation_boundaries is None
                or len(pdf_interpolation_boundaries) != len(pdf_interpolations)
            ):
                raise ValueError("Invalid packed PDF interpolation offsets")
        xs, ys, cs, rules, ends = [], [], [], [], []
        self.offsets = np.r_[0, np.cumsum([len(x) for x in values])].astype(int64)
        roffsets = [0]
        for i, (x, y, c) in enumerate(zip(values, pdfs, cdfs)):
            if (y is None) == (c is None):
                raise ValueError(
                    "Each continuous table requires exactly one PDF or CDF"
                )
            x = np.asarray(x, dtype=float64)
            if isinstance(pdf_interpolations, Integral):
                law, boundary = pdf_interpolations, None
            else:
                a, b = raw_rules[i : i + 2]
                law = np.asarray(pdf_interpolations)[a:b]
                boundary = np.asarray(pdf_interpolation_boundaries)[a:b]
            if c is not None:
                y, c = pdf_from_cdf(x, c)
                law, boundary = INTERPOLATION_HISTOGRAM, None
            else:
                y, c = cdf_from_pdf(x, y, law, boundary)
            law, boundary = _distribution_regions(len(x), law, boundary)
            xs.append(x)
            ys.append(y)
            cs.append(c)
            rules.append(law)
            ends.append(boundary)
            roffsets.append(roffsets[-1] + len(law))
        self.value, self.pdf, self.cdf = map(np.concatenate, (xs, ys, cs))
        self.pdf_interpolations = np.concatenate(rules)
        self.pdf_interpolation_boundaries = np.concatenate(ends)
        self.pdf_interpolation_offsets = np.asarray(roffsets, dtype=int64)

    def __repr__(self):
        return (
            super().__repr__()
            + f"  - grid: {print_1d_array(self.grid)}\n  - PDF tables: {len(self.grid)}\n"
        )


def _distribution_regions(count, laws, boundaries, outer=False):
    """Return validated exclusive interpolation regions for one local grid."""
    if isinstance(laws, Integral):
        laws, boundaries = np.array([laws], dtype=int64), np.array([count], dtype=int64)
    else:
        laws = np.asarray(laws, dtype=int64)
        if boundaries is None:
            raise ValueError("Interpolation boundaries are required for multiple laws")
        boundaries = np.asarray(boundaries, dtype=int64)
    if (
        laws.ndim != 1
        or boundaries.ndim != 1
        or not len(laws)
        or len(laws) != len(boundaries)
        or boundaries[0] < min(2, count)
        or boundaries[-1] != count
        or np.any(np.diff(boundaries) <= 0)
    ):
        raise ValueError("Invalid interpolation regions")
    for law in laws:
        decode_interpolation(law)
    if outer and np.any(
        ~np.isin(laws, [INTERPOLATION_HISTOGRAM, INTERPOLATION_LINEAR])
    ):
        raise ValueError("Multi-table grid interpolation must be histogram or linear")
    return laws, boundaries


# ======================================================================================
# Level scattering
# ======================================================================================


class DistributionLevelScattering(DistributionBase):
    """Discrete-level inelastic-scattering energy distribution.

    Parameters
    ----------
    C1, C2 : float
        Level-scattering law coefficients.
    """

    # MC/DC framework metadata
    label = "level_scattering_distribution"
    sub_type = DISTRIBUTION_LEVEL_SCATTERING

    C1: float
    C2: float

    def __init__(self, C1: float, C2: float) -> None:
        super().__init__()

        self.C1 = C1
        self.C2 = C2

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - C1: {self.C1} [/eV^l]\n"
        text += f"  - C2: {self.C2}\n"
        return text


# ======================================================================================
# Evaporation
# ======================================================================================


class DistributionEvaporation(DistributionBase):
    """Evaporation spectrum with temperature supplied as an evaluable Data object."""

    # MC/DC framework metadata
    label = "evaporation_distribution"
    sub_type = DISTRIBUTION_EVAPORATION

    nuclear_temperature: DataBase
    restriction_energy: float

    def __init__(
        self,
        nuclear_temperature: DataBase,
        restriction_energy: float = 0.0,
    ) -> None:
        super().__init__()

        self.restriction_energy = restriction_energy
        self.nuclear_temperature = nuclear_temperature

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - Restriction energy: {self.restriction_energy} [eV]\n"
        text += f"  - Nuclear temperature data: {self.nuclear_temperature}\n"
        return text


# ======================================================================================
# Maxwellian
# ======================================================================================


class DistributionMaxwellian(DistributionBase):
    """Maxwellian spectrum with temperature supplied as an evaluable Data object."""

    # MC/DC framework metadata
    label = "maxwellian_distribution"
    sub_type = DISTRIBUTION_MAXWELLIAN

    nuclear_temperature: DataBase
    restriction_energy: float

    def __init__(
        self,
        nuclear_temperature: DataBase,
        restriction_energy: float = 0.0,
    ) -> None:
        super().__init__()

        self.restriction_energy = restriction_energy
        self.nuclear_temperature = nuclear_temperature

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - Restriction energy: {self.restriction_energy} [eV]\n"
        text += f"  - Nuclear temperature data: {self.nuclear_temperature}\n"
        return text


# ======================================================================================
# Kalbach-Mann
# ======================================================================================


class DistributionKalbachMann(DistributionBase):
    """Correlated Kalbach-Mann outgoing energy-angle distribution.

    Incident-energy tables are stored in flattened arrays delimited by
    ``offset``. Probability densities are normalized per table.
    """

    # MC/DC framework metadata
    label = "kalbach_mann_distribution"
    sub_type = DISTRIBUTION_KALBACH_MANN

    energy: NDArray[float64]
    offset: NDArray[int64]
    energy_out: NDArray[float64]
    pdf: NDArray[float64]
    cdf: NDArray[float64]
    precompound_factor: NDArray[float64]
    angular_slope: NDArray[float64]
    interpolations: NDArray[int64]
    precompound_tables: list[DataBase]
    angular_slope_tables: list[DataBase]

    def __init__(
        self,
        energy: ArrayLike,
        offset: ArrayLike,
        energy_out: ArrayLike,
        pdf: ArrayLike,
        precompound_factor: ArrayLike,
        angular_slope: ArrayLike,
    ) -> None:
        super().__init__()

        self.energy = np.asarray(energy, dtype=float64)
        self.offset = np.asarray(offset, dtype=int64)
        self.energy_out = np.asarray(energy_out, dtype=float64)
        pdf_array = np.array(pdf, dtype=float64, copy=True)
        self.precompound_factor = np.asarray(
            precompound_factor,
            dtype=float64,
        )
        self.angular_slope = np.asarray(angular_slope, dtype=float64)
        self.interpolations = np.full(
            len(self.energy_out), INTERPOLATION_LINEAR, dtype=int64
        )
        self.precompound_tables = []
        self.angular_slope_tables = []

        self.pdf, self.cdf = multi_cdf_from_pdf(
            self.offset,
            self.energy_out,
            pdf_array,
        )

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - grid {print_1d_array(self.energy)} [eV]\n"
        text += f"  - offset {print_1d_array(self.offset)}\n"
        text += f"  - energy {print_1d_array(self.energy_out)} [eV]\n"
        text += f"  - energy-pdf {print_1d_array(self.pdf)} [/eV]\n"
        text += "  - precompound factor " f"{print_1d_array(self.precompound_factor)}\n"
        text += f"  - angular slope {print_1d_array(self.angular_slope)}\n"
        return text


# ======================================================================================
# Table energy-angle
# ======================================================================================


class DistributionTableEnergyAngle(DistributionBase):
    """Correlated tabulated outgoing energy and scattering-angle distribution.

    ``offset`` delimits outgoing-energy tables and ``cosine_offset`` delimits
    conditional cosine tables in the flattened arrays.
    """

    # MC/DC framework metadata
    label = "table_energy_angle_distribution"
    sub_type = DISTRIBUTION_TABLE_ENERGY_ANGLE

    energy: NDArray[float64]
    offset: NDArray[int64]
    energy_out: NDArray[float64]
    pdf: NDArray[float64]
    cdf: NDArray[float64]
    cosine_offset_: NDArray[int64]
    cosine: NDArray[float64]
    cosine_pdf: NDArray[float64]
    cosine_cdf: NDArray[float64]
    interpolations: NDArray[int64]
    cosine_interpolations: NDArray[int64]

    def __init__(
        self,
        energy: ArrayLike,
        offset: ArrayLike,
        energy_out: ArrayLike,
        pdf: ArrayLike,
        cosine_offset: ArrayLike,
        cosine: ArrayLike,
        cosine_pdf: ArrayLike,
    ) -> None:
        super().__init__()

        self.energy = np.asarray(energy, dtype=float64)
        self.offset = np.asarray(offset, dtype=int64)
        self.energy_out = np.asarray(energy_out, dtype=float64)
        pdf_array = np.array(pdf, dtype=float64, copy=True)
        self.cosine_offset_ = np.asarray(cosine_offset, dtype=int64)
        self.cosine = np.asarray(cosine, dtype=float64)
        cosine_pdf_array = np.asarray(cosine_pdf, dtype=float64)
        self.interpolations = np.full(
            len(self.energy_out), INTERPOLATION_LINEAR, dtype=int64
        )
        self.cosine_interpolations = np.full(
            len(self.cosine), INTERPOLATION_LINEAR, dtype=int64
        )

        self.pdf, self.cdf = multi_cdf_from_pdf(
            self.offset,
            self.energy_out,
            pdf_array,
        )

        self.cosine_pdf = cosine_pdf_array.copy()
        self.cosine_cdf = np.zeros_like(self.cosine_pdf)

        for i in range(len(self.offset)):
            energy_start = self.offset[i]

            if i + 1 < len(self.offset):
                energy_stop = self.offset[i + 1]
            else:
                energy_stop = len(self.energy_out)

            inner_offset = self.cosine_offset_[energy_start:energy_stop]

            if len(inner_offset) == 0:
                print_error(
                    "Each incident-energy table must reference at least one "
                    "cosine distribution."
                )

            cosine_start = inner_offset[0]

            if i + 1 < len(self.offset):
                cosine_stop = self.cosine_offset_[energy_stop]
            else:
                cosine_stop = len(self.cosine)

            inner_offset_local = inner_offset - cosine_start

            (
                self.cosine_pdf[cosine_start:cosine_stop],
                self.cosine_cdf[cosine_start:cosine_stop],
            ) = multi_cdf_from_pdf(
                inner_offset_local,
                self.cosine[cosine_start:cosine_stop],
                cosine_pdf_array[cosine_start:cosine_stop],
            )

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - grid {print_1d_array(self.energy)} [eV]\n"
        text += f"  - offset {print_1d_array(self.offset)}\n"
        text += f"  - energy {print_1d_array(self.energy_out)} [eV]\n"
        text += f"  - energy-pdf {print_1d_array(self.pdf)} [/eV]\n"
        text += "  - cosine-offset " f"{print_1d_array(self.cosine_offset_)}\n"
        text += f"  - cosine {print_1d_array(self.cosine)}\n"
        text += f"  - cosine-pdf {print_1d_array(self.cosine_pdf)}\n"
        return text


# ======================================================================================
# N-body
# ======================================================================================


class DistributionNBody(DistributionBase):
    """N-body phase-space outgoing-energy distribution.

    Parameters
    ----------
    values, probabilities : array_like
        Outgoing values and their piecewise-linear relative density.
    """

    # MC/DC framework metadata
    label = "nbody_distribution"
    sub_type = DISTRIBUTION_N_BODY

    pdf: DataTable

    def __init__(
        self,
        values: ArrayLike,
        probabilities: ArrayLike,
    ) -> None:
        super().__init__()

        value_array = np.asarray(values, dtype=float64)
        probability_array = np.asarray(probabilities, dtype=float64)

        if value_array.ndim != 1:
            print_error("values must be one-dimensional.")
        if probability_array.ndim != 1:
            print_error("probabilities must be one-dimensional.")
        if len(value_array) != len(probability_array):
            print_error("values and probabilities must have the same length.")

        pdf_normalized, cdf_normalized = cdf_from_pdf(
            value_array,
            probability_array,
        )

        self.pdf = DataTable(
            value_array,
            pdf_normalized,
            INTERPOLATION_LINEAR,
            aux=cdf_normalized,
        )

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - value {print_1d_array(self.pdf.x)}\n"
        text += f"  - pdf {print_1d_array(self.pdf.y)}\n"
        return text
