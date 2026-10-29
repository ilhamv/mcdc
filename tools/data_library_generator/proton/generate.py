"""Prepare MC/DC proton data from ACE and NIST PSTAR inputs.

Writes temperature-independent ``proton-{nuclide}.h5`` files in the mcdc-data
format. Source parsing and distribution conversion are shared with the neutron
and photon generators; this module handles proton reaction selection and
secondary-product records.

Environment variables
    MCDC_ACELIB_PROTON  directory with proton ACE files (*.ace)
    MCDC_PSTAR_LIB      directory with NIST PSTAR stopping-power tables (<Symbol>.txt)
    MCDC_LIB_PROTON     output directory for the MC/DC HDF5 files

Usage
    python generate.py                 # convert only nuclides with no output yet
    python generate.py --rewrite       # regenerate everything
    python generate.py --verbose       # detailed per-reaction output
    python generate.py --only p-Al027.ace [more.ace ...]
"""

import argparse
import os
from pathlib import Path
import sys
from dataclasses import dataclass

import ACEtk
import h5py
import numpy as np
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from data_library_generator.common.ace_converter import (
    EV_PER_MEV,
    read_source_libraries,
    write_angular_distribution as write_proton_angular_distribution,
    write_energy_distribution as write_proton_energy_distribution,
)
from data_library_generator.common.batch import ConversionErrors, add_keep_going_option
from data_library_generator.common.hdf5_writer import (
    atomic_file,
    interpolation,
    prepared_file_exists,
    quantity,
)

# ==============================================================================
# Constants
# ==============================================================================

ELEMENT_SYMBOLS = (  # index 0 -> Z=1
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)

# Radiation length [g/cm2] for Z = 1..118 (same order as ELEMENT_SYMBOLS)
RADIATION_LENGTH = (
    63.04,
    94.32,
    82.77,
    65.19,
    52.68,
    42.70,
    37.99,
    34.24,
    32.93,
    28.93,
    27.74,
    25.03,
    24.01,
    21.82,
    21.21,
    19.50,
    19.28,
    19.55,
    17.32,
    16.14,
    16.55,
    16.16,
    15.84,
    14.94,
    14.64,
    13.84,
    13.62,
    12.68,
    12.86,
    12.43,
    12.47,
    12.25,
    11.94,
    11.91,
    11.42,
    11.37,
    11.03,
    10.76,
    10.41,
    10.20,
    9.92,
    9.80,
    9.58,
    9.48,
    9.27,
    9.20,
    8.97,
    9.00,
    8.85,
    8.82,
    8.73,
    8.83,
    8.48,
    8.48,
    8.31,
    8.31,
    8.14,
    7.96,
    7.76,
    7.71,
    7.51,
    7.57,
    7.44,
    7.48,
    7.36,
    7.32,
    7.23,
    7.14,
    7.03,
    7.02,
    6.92,
    6.89,
    6.82,
    6.76,
    6.69,
    6.68,
    6.59,
    6.54,
    6.46,
    6.44,
    6.42,
    6.37,
    6.29,
    6.16,
    6.07,
    6.28,
    6.19,
    6.15,
    6.06,
    6.07,
    5.93,
    6.00,
    5.87,
    5.93,
    5.80,
    5.79,
    5.69,
    5.68,
    5.61,
    5.62,
    5.55,
    5.48,
    5.45,
    5.47,
    5.40,
    5.34,
    5.27,
    5.17,
    5.26,
    5.24,
    5.18,
    5.16,
    5.10,
    5.08,
    5.01,
    5.00,
    4.95,
    4.88,
)

SYMBOL_TO_Z = {symbol: z for z, symbol in enumerate(ELEMENT_SYMBOLS, start=1)}

DEUTERIUM_RADIATION_LENGTH = 125.98  # g/cm2, for (Z, A) == (1, 2)

# Isotopes to generate stopping-power-only HDF5 files for when no ACE file exists
# (TENDL excludes H and He because TALYS does not apply).
# Format: (symbol, A, atomic_weight_ratio); AWR = atomic mass / neutron mass.
NEUTRON_MASS_U = 1.008664916
PSTAR_ONLY_ISOTOPES = (
    ("H", 1, 1.00794 / NEUTRON_MASS_U),  # natural H ~ H-1
    ("H", 2, 2.01410 / NEUTRON_MASS_U),  # deuterium
    ("He", 3, 3.01603 / NEUTRON_MASS_U),
    ("He", 4, 4.00260 / NEUTRON_MASS_U),
)

# --- Reaction bookkeeping -------------------------------------------------------
# MTs that are never written (total, non-elastic, "(z,anything)", disappearance).
ALWAYS_REDUNDANT_MTS = frozenset({1, 3, 10, 101})

# A "sum" MT is redundant ONLY if (some of) its partial reactions are present in the
# same file; otherwise the sum is the only record of that channel and must be kept.
# (TENDL-2021 proton files contain no partials at all, so every sum is kept there.)
SUM_MT_PARTIALS = {
    4: range(50, 92),  # (z,n) <- levels
    16: range(875, 892),  # (z,2n) <- levels
    103: range(600, 650),  # (z,p)
    104: range(650, 700),  # (z,d)
    105: range(700, 750),  # (z,t)
    106: range(750, 800),  # (z,He3)
    107: range(800, 850),  # (z,alpha)
}
MAX_MT = 891  # reactions above this (e.g. 9xx partial capture) are not written

# --- Secondary particles --------------------------------------------------------
# ACE particle-type id (PTYPE block) -> (product name, ZAP = 1000*Z + A)
ACE_PARTICLE_TO_PRODUCT = {
    1: ("neutron", 1),
    9: ("proton", 1001),
    31: ("deuteron", 1002),
    32: ("triton", 1003),
    33: ("helion", 2003),
    34: ("alpha", 2004),
}

# Photons are stored in the separate MTRP/SIGP/DLWP/ANDP blocks.
PHOTON_NAME = "photon"

# ACE photon-production reactions are numbered MTRP = 1000*MT + photon_index.
PHOTON_MT_FACTOR = 1000


class ConversionError(Exception):
    """The ACE data are not what this converter supports."""


def fail(message):
    raise ConversionError(message)


def note(message):
    print(f"  [note] {message}")


# ==============================================================================
# Small helpers
# ==============================================================================


def decode_ace_zaid(zaid):
    """Return (Z, A, S) from an ACE ZAID string such as '13027.00h'."""
    za = int(zaid.strip().split(".")[0])
    S = 0
    if za >= 600000:  # metastable convention: A + 300 + 100*S
        S = (za % 1000) // 400
        za = za - S * 400
    return za // 1000, za % 1000, S


def radiation_length(Z, A):
    if (Z, A) == (1, 2):
        return DEUTERIUM_RADIATION_LENGTH
    return RADIATION_LENGTH[Z - 1] if 1 <= Z <= len(RADIATION_LENGTH) else None


def frame_name(frame):
    if frame == ACEtk.ReferenceFrame.Laboratory:
        return "LAB"
    if frame == ACEtk.ReferenceFrame.CentreOfMass:
        return "COM"
    fail(f"Unknown reference frame {frame!r}")


def load_pstar_file(filepath):
    """Two-column PSTAR file (Energy MeV, Stopping power MeV cm2/g)."""
    energies, stopping = [], []
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            try:
                energies.append(float(parts[0]))
                stopping.append(float(parts[1]))
            except ValueError:
                continue
    return np.array(energies), np.array(stopping)


def write_stopping_power(h5, pstar_dir, symbol, verbose=False):
    """Write PSTAR mass stopping power using the mcdc-data representation."""
    path = Path(pstar_dir) / f"{symbol}.txt"
    if not path.is_file():
        if verbose:
            print(f"  [warn] No PSTAR file for {symbol}")
        return False
    energy, stopping = load_pstar_file(path)
    if (
        len(energy) < 2
        or energy.shape != stopping.shape
        or not np.all(np.isfinite(energy))
        or not np.all(np.isfinite(stopping))
        or np.any(energy <= 0)
        or np.any(np.diff(energy) <= 0)
        or np.any(stopping < 0)
    ):
        raise ValueError(f"{path}: invalid PSTAR stopping-power table")
    group = h5.create_group("stopping-power")
    group.attrs["type"] = "table"
    quantity(group, "grid", energy * EV_PER_MEV, "eV")
    quantity(group, "value", stopping * EV_PER_MEV, "eV cm2/g")
    interpolation(group, len(energy))
    return True


# ==============================================================================
# Tabulated yields (multiplicity as a function of incident energy)
# ==============================================================================


def read_tabulated_yield(data, what):
    """
    Read an ACEtk TabulatedSecondaryParticleMultiplicity (SIGH/SIGP entries) or
    TabulatedMultiplicity (primary DLW) into (energy [MeV], value) arrays.
    Only linear-linear interpolation is supported (all TENDL-2021 files use it).
    """
    offset = (
        2
        if isinstance(data, ACEtk.continuous.TabulatedSecondaryParticleMultiplicity)
        else 0
    )
    if not data.interpolation_data.is_linear_linear:
        regions = int(data.NB)
        boundaries = [int(data.xss(offset + index)) for index in range(2, regions + 2)]
        interpolants = [
            int(data.xss(offset + index))
            for index in range(regions + 2, 2 * regions + 2)
        ]
        fail(
            f"{what}: yield table is not linear-linear "
            f"(NBT={boundaries}, INT={interpolants})"
        )
    regions = int(data.NB)
    count = int(data.NE)
    energy_start = offset + 2 * regions + 3
    value_start = energy_start + count
    energy = np.array(
        [data.xss(index) for index in range(energy_start, value_start)],
        dtype=np.float64,
    )
    value = np.array(
        [data.xss(index) for index in range(value_start, value_start + count)],
        dtype=np.float64,
    )
    if len(energy) != len(value) or len(energy) < 2:
        fail(
            f"{what}: malformed yield table ({len(energy)} energies, {len(value)} values)"
        )
    if np.any(np.diff(energy) < 0.0):
        fail(f"{what}: yield energies are not non-decreasing")
    if np.any(value < 0.0):
        fail(f"{what}: negative yield")
    return energy, value


def write_yield_table(group, name, energy, value):
    """Write an incident-energy-dependent yield as a standard Data table."""
    energy = np.asarray(energy, dtype=np.float64)
    value = np.asarray(value, dtype=np.float64)
    if (
        energy.ndim != 1
        or energy.shape != value.shape
        or len(energy) < 2
        or not np.all(np.isfinite(energy))
        or not np.all(np.isfinite(value))
        or np.any(np.diff(energy) < 0)
        or np.any(value < 0)
    ):
        fail(f"{name}: invalid yield table")
    same_as_previous = np.zeros(len(energy), dtype=bool)
    same_as_previous[1:] = (energy[1:] == energy[:-1]) & (value[1:] == value[:-1])
    energy, value = energy[~same_as_previous], value[~same_as_previous]
    if len(energy) < 2 or not np.any(np.diff(energy) > 0):
        fail(f"{name}: yield table must contain at least two distinct energies")

    table = group.create_group(name)
    table.attrs["type"] = "table"
    quantity(table, "grid", energy * EV_PER_MEV, "eV")
    quantity(table, "value", value)
    interpolation(table, len(energy))
    return table


# ==============================================================================
# Angular and energy distributions
# ==============================================================================


def write_cosine_distribution(data, group):
    """Write an angular Distribution; return False when energy supplies the angle."""
    if isinstance(data, ACEtk.continuous.DistributionGivenElsewhere):
        return False
    write_proton_angular_distribution(group, data)
    return True


def write_energy_distribution(data, group):
    """Write an ACE energy or joint energy-angle Distribution."""
    write_proton_energy_distribution(group, data)
    return group.attrs["type"]


def is_joint_energy_distribution(data):
    """Return whether an ACE spectrum supplies both outgoing energy and angle."""
    C = ACEtk.continuous
    joint_types = (
        C.KalbachMannDistributionData,
        C.EnergyAngleDistributionData,
        C.NBodyPhaseSpaceDistribution,
    )
    if isinstance(data, C.MultiDistributionData):
        joint = [
            is_joint_energy_distribution(data.distribution(index))
            for index in range(1, data.number_distributions + 1)
        ]
        if any(joint) and not all(joint):
            fail(
                "A product mixture cannot combine joint and independent angle-energy laws"
            )
        return all(joint)
    return isinstance(data, joint_types)


# ==============================================================================
# Secondary particles (neutrons, d, t, He3, alpha and photons)
# ==============================================================================


def write_secondary_product(
    parent, name, *, frame, yield_energy, yield_value, energy, angle
):
    """Write one additional ACE product as a standard CrossProduct record."""
    C = ACEtk.continuous
    joint = is_joint_energy_distribution(energy)
    index = str(len(parent))
    product = parent.create_group(index)
    quantity(product, "particle", name)
    quantity(product, "reference-frame", frame)
    write_yield_table(product, "yield", yield_energy, yield_value)
    quantity(product, "decay-rate", np.inf, "1/s")

    distribution = product.create_group("angle-energy-distribution")
    if joint:
        write_energy_distribution(energy, distribution)
    else:
        if isinstance(angle, C.DistributionGivenElsewhere):
            fail(
                f"{parent.parent.name}/{name}: angle is not included in its energy distribution"
            )
        distribution.attrs["type"] = "independent-angle-energy"
        write_cosine_distribution(
            angle, distribution.create_group("angular-cosine-distribution")
        )
        write_energy_distribution(
            energy, distribution.create_group("energy-distribution")
        )
    return product


def _products_group(MT, reaction_groups, skipped_MTs, what, verbose):
    """Return a reaction's cross-products group or skip a redundant sum."""
    if MT in reaction_groups:
        return reaction_groups[MT]["cross-products"]
    if MT in skipped_MTs:
        if verbose:
            note(f"{what}: MT-{MT:03} is a redundant sum in this file; skipped")
        return None
    fail(f"{what}: MT-{MT:03} has no reaction group in the output")


def write_secondary_particles(table, reaction_groups, skipped_MTs, verbose=False):
    """
    Write all non-proton products of all reactions.
      reaction_groups: {MT: HDF5 group of that reaction}
      skipped_MTs:     MTs deliberately left out of the output (redundant sums)
    Returns the number of products written.
    """
    C = ACEtk.continuous
    written = 0

    # ---- Particles with their own blocks: neutrons, d, t, He3, alpha ----------
    for p in range(1, table.NTYPE + 1):
        ace_id = table.PTYPE.particle_type(p)
        if ace_id not in ACE_PARTICLE_TO_PRODUCT:
            fail(f"Unknown ACE secondary-particle id {ace_id}")
        if ace_id == 9:
            fail(
                "Proton listed as a secondary particle type; protons are the primary here"
            )
        name, zap = ACE_PARTICLE_TO_PRODUCT[ace_id]

        n_rx = table.NTRO.NP(p)
        mtrh, tyrh, sigh = table.MTRH(p), table.TYRH(p), table.SIGH(p)
        dlwh, andh = table.DLWH(p), table.ANDH(p)
        sizes = (
            mtrh.number_reactions, tyrh.number_reactions, len(sigh.data),
            dlwh.number_reactions, andh.number_reactions,
        )  # fmt: skip
        if any(s != n_rx for s in sizes):
            fail(f"{name}: inconsistent block sizes {sizes}, expected {n_rx} each")

        for j in range(1, n_rx + 1):
            MT = mtrh.MT(j)
            what = f"{name} from MT-{MT:03}"
            if tyrh.multiplicity(j) != 1:
                fail(f"{what}: TYRH multiplicity {tyrh.multiplicity(j)} (expected 1)")
            yield_data = sigh.data[j - 1]
            if yield_data.MT != MT:
                fail(f"{what}: SIGH entry is for MT-{yield_data.MT:03}, not MT-{MT:03}")
            if yield_data.MFTYPE != 12:
                fail(
                    f"{what}: SIGH MFTYPE {yield_data.MFTYPE} (only 12, a yield, is supported)"
                )
            group = _products_group(MT, reaction_groups, skipped_MTs, what, verbose)
            if group is None:
                continue

            y_energy, y_value = read_tabulated_yield(yield_data, what)
            write_secondary_product(
                group, name,
                frame=frame_name(tyrh.reference_frame(j)),
                yield_energy=y_energy, yield_value=y_value,
                energy=dlwh.energy_distribution_data(j),
                angle=andh.angular_distribution_data(j),
            )  # fmt: skip
            written += 1
            if verbose:
                print(
                    f"      MT-{MT:03} {name}: yield {y_value.min():g}..{y_value.max():g}"
                )

    # ---- Photons (separate blocks MTRP / SIGP / DLWP / ANDP) ------------------
    n_photon = table.NTRP
    if n_photon:
        mtrp, sigp, dlwp, andp = table.MTRP, table.SIGP, table.DLWP, table.ANDP
        sizes = (
            mtrp.number_reactions, len(sigp.data),
            dlwp.number_reactions, andp.number_reactions,
        )  # fmt: skip
        if any(s != n_photon for s in sizes):
            fail(f"photon: inconsistent block sizes {sizes}, expected {n_photon} each")

        for k in range(1, n_photon + 1):
            MT, index = divmod(mtrp.MT(k), PHOTON_MT_FACTOR)
            name = PHOTON_NAME if index == 1 else f"{PHOTON_NAME}_{index}"
            what = f"{name} from MT-{MT:03}"
            yield_data = sigp.data[k - 1]
            if not isinstance(yield_data, C.TabulatedSecondaryParticleMultiplicity):
                fail(
                    f"{what}: SIGP entry is a {type(yield_data).__name__}, not a yield table"
                )
            if yield_data.MT != MT or yield_data.MFTYPE not in (12, 16):
                fail(
                    f"{what}: unexpected SIGP entry (MT={yield_data.MT}, MFTYPE={yield_data.MFTYPE})"
                )
            group = _products_group(MT, reaction_groups, skipped_MTs, what, verbose)
            if group is None:
                continue

            # ACE photon production has no reference-frame flag (ACEtk gives None):
            # photons are emitted in the laboratory frame.
            frame = dlwp.reference_frame(k)
            y_energy, y_value = read_tabulated_yield(yield_data, what)
            write_secondary_product(
                group, name,
                frame="LAB" if frame is None else frame_name(frame),
                yield_energy=y_energy, yield_value=y_value,
                energy=dlwp.energy_distribution_data(k),
                angle=andp.angular_distribution_data(k),
            )  # fmt: skip
            written += 1
            if verbose:
                print(
                    f"      MT-{MT:03} {name}: yield {y_value.min():g}..{y_value.max():g}"
                )
    return written


# ==============================================================================
# Primary reactions
# ==============================================================================


@dataclass
class Reaction:
    MT: int
    index: int  # 1-based position in the MTR / TYR / LQR / SIG blocks
    group: str  # "capture" (no outgoing proton) or "inelastic_scattering"
    multiplicity: int = 0  # fixed number of outgoing protons
    proton_yield: tuple = (
        None  # (energy, value) when the proton yield depends on energy
    )


def classify_reactions(table):
    """
    Sort the file's reactions into capture / inelastic and find the redundant ones.
    Returns (reactions, skipped_MTs).

    A reaction with outgoing-proton multiplicity 0 is "capture" (the proton is
    absorbed; any secondaries come from the product yields).  A positive multiplicity
    is "inelastic".  |TYR| > 100 means the proton yield is energy dependent; that
    table lives in the primary DLW block (MT-5 in TENDL files).
    """
    mtr = table.reaction_number_block
    tyr = table.frame_and_multiplicity_block
    dlw = table.energy_distribution_block
    n = tyr.number_reactions
    present = {mtr.MT(i) for i in range(1, n + 1)}

    skipped = set(ALWAYS_REDUNDANT_MTS)
    for sum_MT, partials in SUM_MT_PARTIALS.items():
        if sum_MT in present and present.intersection(partials):
            skipped.add(sum_MT)
    skipped |= {MT for MT in present if MT > MAX_MT}

    if (
        present.intersection((18, 19, 20, 21, 38))
        or table.fission_multiplicity_block is not None
    ):
        fail("Fission data in a proton ACE file are not supported by this converter")

    reactions = []
    for i in range(1, n + 1):
        MT = mtr.MT(i)
        if MT == 2 or MT in skipped:  # elastic comes from the principal xs block
            continue
        raw = tyr.multiplicity(i)
        if raw == 0:
            reactions.append(Reaction(MT, i, "capture"))
        elif raw > 100:
            if i > dlw.number_reactions:
                fail(f"MT-{MT:03}: energy-dependent multiplicity but no DLW entry")
            data = dlw.multiplicity_data(i)
            if not isinstance(data, ACEtk.continuous.TabulatedMultiplicity):
                fail(
                    f"MT-{MT:03}: TYR={raw} but multiplicity data is {type(data).__name__}"
                )
            reactions.append(
                Reaction(MT, i, "inelastic_scattering",
                         proton_yield=read_tabulated_yield(data, f"proton yield of MT-{MT:03}"))
            )  # fmt: skip
        else:
            reactions.append(Reaction(MT, i, "inelastic_scattering", multiplicity=raw))
    return reactions, skipped


def write_reaction(group, table, rx, verbose):
    """Write one proton reaction in the standard mcdc-data layout."""
    tyr = table.frame_and_multiplicity_block
    xs_block = table.cross_section_block
    group.attrs["ENDF-MT"] = rx.MT
    group.attrs["category"] = rx.group.replace("_", "-")

    xs = group.create_group("xs")
    quantity(xs, "value", xs_block.cross_sections(rx.index), "barn")
    quantity(xs, "offset", xs_block.energy_index(rx.index) - 1)
    interpolation(xs, len(xs["value"]))
    quantity(
        group,
        "Q-value",
        table.reaction_qvalue_block.q_value(rx.index) * EV_PER_MEV,
        "eV",
    )
    quantity(group, "reference-frame", frame_name(tyr.reference_frame(rx.index)))
    group.create_group("cross-products")

    if rx.group == "capture":
        return

    if rx.proton_yield is None:
        quantity(group, "multiplicity", rx.multiplicity)
    else:
        write_yield_table(group, "multiplicity", *rx.proton_yield)

    angle = table.angular_distribution_block.angular_distribution_data(rx.index)
    if not isinstance(angle, ACEtk.continuous.DistributionGivenElsewhere):
        write_cosine_distribution(
            angle, group.create_group("angular-cosine-distribution")
        )

    spectra = table.energy_distribution_block.energy_distribution_data(rx.index)
    write_energy_distribution(spectra, group.create_group("energy-distribution"))


# ==============================================================================
# Per-file processing
# ==============================================================================


def source_libraries(table, include_pstar, *, source_library=None):
    """Build ACE provenance through the common reader and add PSTAR when used."""
    sources = read_source_libraries(table, source_library=source_library)
    if include_pstar:
        if "PSTAR" in sources:
            raise ValueError("ACE source-library label conflicts with PSTAR provenance")
        sources["PSTAR"] = {
            "title": "NIST PSTAR stopping-power tables",
            "version": "",
            "release-date": "",
            "reference": "NIST PSTAR",
            "comments": "",
            "format": "PSTAR",
            "format-version": "",
            "processing-date": "",
        }
    return sources


def _write_target_properties(file, table, Z, A, state):
    quantity(file, "atomic-number", Z)
    quantity(file, "mass-number", A)
    quantity(file, "atomic-weight-ratio", table.atomic_weight_ratio)
    quantity(file, "excitation-level", state)
    quantity(file, "fissionable", False)
    quantity(file, "radiation-length", radiation_length(Z, A), "g/cm2")


def _write_reactions(file, table, verbose):
    reactions, skipped = classify_reactions(table)
    energy = np.asarray(table.principal_cross_section_block.energies) * EV_PER_MEV
    if (
        energy.ndim != 1
        or len(energy) < 2
        or not np.all(np.isfinite(energy))
        or np.any(np.diff(energy) < 0)
    ):
        fail("Invalid principal proton energy grid")
    quantity(file, "xs-energy-grid", energy, "eV")
    root = file.create_group("reactions", track_order=True)

    elastic = root.create_group("MT-002", track_order=True)
    elastic.attrs["ENDF-MT"] = 2
    elastic.attrs["category"] = "elastic-scattering"
    xs = elastic.create_group("xs")
    quantity(xs, "value", table.principal_cross_section_block.elastic, "barn")
    quantity(xs, "offset", 0)
    interpolation(xs, len(energy))
    quantity(elastic, "reference-frame", "COM")
    quantity(elastic, "Q-value", 0.0, "eV")
    elastic.create_group("cross-products")
    angle = table.angular_distribution_block.angular_distribution_data(0)
    if isinstance(angle, ACEtk.continuous.DistributionGivenElsewhere):
        fail("MT-002 elastic scattering has no explicit angular distribution")
    write_cosine_distribution(
        angle, elastic.create_group("angular-cosine-distribution")
    )
    reaction_groups = {2: elastic}

    for reaction in reactions:
        group = root.create_group(f"MT-{reaction.MT:03}", track_order=True)
        write_reaction(group, table, reaction, verbose)
        reaction_groups[reaction.MT] = group
    count = write_secondary_particles(table, reaction_groups, skipped, verbose)
    if verbose:
        print(f"  Reactions written: {len(reactions) + 1}")
        print(f"  Redundant or unsupported MTs skipped: {sorted(skipped)}")
        print(f"  Cross-products written: {count}")


def write_nuclide(file, table, Z, A, state, pstar_path, verbose):
    """Write target properties, optional stopping power, and proton reactions."""
    _write_target_properties(file, table, Z, A, state)
    if pstar_path is not None:
        write_stopping_power(file, pstar_path.parent, ELEMENT_SYMBOLS[Z - 1], verbose)
    _write_reactions(file, table, verbose)


def process_ace_file(
    ace_path,
    output_dir,
    pstar_dir=None,
    verbose=False,
    rewrite=False,
    *,
    source_library=None,
):
    """Convert one ACE evaluation into ``proton-{nuclide}.h5``."""
    ace_path, output_dir = Path(ace_path), Path(output_dir)
    table = ACEtk.ContinuousEnergyTable.from_file(str(ace_path))
    header = table.header
    Z, A, S = decode_ace_zaid(header.zaid)
    if not 1 <= Z <= len(ELEMENT_SYMBOLS):
        fail(f"Unsupported Z={Z} in {ace_path}")
    symbol = ELEMENT_SYMBOLS[Z - 1]
    nuclide = f"{symbol}{A}" if S == 0 else f"{symbol}{A}m{S}"
    out_name = f"proton-{nuclide}.h5"
    out_path = output_dir / out_name
    pstar_path = Path(pstar_dir) / f"{symbol}.txt" if pstar_dir is not None else None
    if pstar_path is not None and not pstar_path.is_file():
        pstar_path = None
        if verbose:
            print(f"  [warn] No PSTAR file for {symbol}")
    sources = source_libraries(
        table, pstar_path is not None, source_library=source_library
    )
    fields = (
        "atomic-number",
        "mass-number",
        "atomic-weight-ratio",
        "excitation-level",
        "fissionable",
        "radiation-length",
        "xs-energy-grid",
        "reactions",
    )
    selected, _ = classify_reactions(table)
    required = list(fields)
    for reaction in selected:
        base = f"reactions/MT-{reaction.MT:03}"
        required.extend(
            f"{base}/{field}"
            for field in (
                "xs/value",
                "xs/offset",
                "reference-frame",
                "Q-value",
                "cross-products",
            )
        )
        if reaction.group == "inelastic":
            required.extend((f"{base}/multiplicity", f"{base}/energy-distribution"))
    required.extend(
        (
            "reactions/MT-002/xs/value",
            "reactions/MT-002/xs/offset",
            "reactions/MT-002/reference-frame",
            "reactions/MT-002/Q-value",
            "reactions/MT-002/cross-products",
            "reactions/MT-002/angular-cosine-distribution",
        )
    )
    if pstar_path is not None:
        required.extend(("stopping-power/grid", "stopping-power/value"))
    if not rewrite and prepared_file_exists(out_path, required, sources):
        if verbose:
            print(f"Skip {out_name}: output file already exists")
        return out_name
    if verbose:
        print(f"Create {out_name} from {ace_path.name}")
    with atomic_file(out_path, sources) as file:
        write_nuclide(file, table, Z, A, S, pstar_path, verbose)
    return out_name


def process_pstar_only_file(
    symbol, A, awr, output_dir, pstar_dir, verbose=False, rewrite=False
):
    """Write a stopping-power-only proton file when ACE data are unavailable."""
    pstar_path = Path(pstar_dir) / f"{symbol}.txt"
    if not pstar_path.is_file():
        return None
    Z = SYMBOL_TO_Z[symbol]
    out_name = f"proton-{symbol}{A}.h5"
    out_path = Path(output_dir) / out_name
    if out_path.is_file():
        with h5py.File(out_path, "r") as file:
            if "reactions" in file:
                return out_name
    sources = {
        "PSTAR": {
            "title": "NIST PSTAR stopping-power tables",
            "version": "",
            "release-date": "",
            "reference": "NIST PSTAR",
            "comments": "",
            "format": "PSTAR",
            "format-version": "",
            "processing-date": "",
        }
    }
    required = (
        "atomic-number",
        "mass-number",
        "atomic-weight-ratio",
        "excitation-level",
        "fissionable",
        "radiation-length",
        "stopping-power/grid",
        "stopping-power/value",
    )
    if not rewrite and prepared_file_exists(out_path, required, sources):
        return out_name
    if verbose:
        print(f"Create {out_name} from PSTAR")
    with atomic_file(out_path, sources) as file:
        quantity(file, "atomic-number", Z)
        quantity(file, "mass-number", A)
        quantity(file, "atomic-weight-ratio", awr)
        quantity(file, "excitation-level", 0)
        quantity(file, "fissionable", False)
        quantity(file, "radiation-length", radiation_length(Z, A), "g/cm2")
        write_stopping_power(file, pstar_path.parent, symbol, verbose)
    return out_name


# ==============================================================================
# Main
# ==============================================================================


def main():
    """Prepare selected ACE and optional PSTAR inputs with common error handling.

    Stop on the first conversion failure unless --keep-going is supplied.
    In continuation mode both ACE and PSTAR-only failures are collected and
    reported together with a nonzero exit status after all inputs are tried.
    """
    parser = argparse.ArgumentParser(
        description="Prepare MC/DC proton data from ACE and NIST PSTAR"
    )
    parser.add_argument(
        "--rewrite", action="store_true", help="regenerate existing files"
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Report files and reaction details"
    )
    parser.add_argument(
        "--only",
        nargs="+",
        metavar="ACE_FILE",
        help="convert only these ACE file names",
    )
    parser.add_argument(
        "--source-library",
        metavar="NAME",
        help="Label ACE source provenance; default: recognize the header or use unknown",
    )
    add_keep_going_option(parser)
    args = parser.parse_args()

    for variable in ("MCDC_ACELIB_PROTON", "MCDC_PSTAR_LIB", "MCDC_LIB_PROTON"):
        if not os.getenv(variable):
            parser.error(f"Environment variable {variable} is not set")
    ace_dir = Path(os.environ["MCDC_ACELIB_PROTON"])
    pstar_dir = Path(os.environ["MCDC_PSTAR_LIB"])
    output_dir = Path(os.environ["MCDC_LIB_PROTON"])
    if not ace_dir.is_dir():
        parser.error(f"ACE directory does not exist: {ace_dir}")
    if not pstar_dir.is_dir():
        parser.error(f"PSTAR directory does not exist: {pstar_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    ace_files = sorted(path for path in ace_dir.iterdir() if path.suffix == ".ace")
    if args.only:
        available = {path.name for path in ace_files}
        missing = [name for name in args.only if name not in available]
        if missing:
            parser.error(f"ACE filenames not found in {ace_dir}: {', '.join(missing)}")
        selected = set(args.only)
        ace_files = [path for path in ace_files if path.name in selected]

    print(
        f"ACE directory: {ace_dir}\nPSTAR directory: {pstar_dir}\n"
        f"Output directory: {output_dir}\n"
    )
    failures = ConversionErrors(parser, keep_going=args.keep_going)
    progress = tqdm(ace_files, disable=args.verbose, unit="file")
    for path in progress:
        progress.set_postfix_str(path.name)
        try:
            process_ace_file(
                path,
                output_dir,
                pstar_dir,
                args.verbose,
                args.rewrite,
                source_library=args.source_library,
            )
        except (ConversionError, ValueError, OSError, AttributeError) as error:
            failures.handle(path.name, error, progress=progress)

    for symbol, A, awr in PSTAR_ONLY_ISOTOPES:
        filename = f"proton-{symbol}{A}.h5"
        try:
            out = process_pstar_only_file(
                symbol, A, awr, output_dir, pstar_dir, args.verbose, args.rewrite
            )
            if out is None and args.verbose:
                print(f"Skip {filename}: no PSTAR data")
        except (ValueError, OSError) as error:
            failures.handle(filename, error)

    failures.finish()
    print(f"\nDone. {len(ace_files)} ACE files processed.")


if __name__ == "__main__":
    main()
