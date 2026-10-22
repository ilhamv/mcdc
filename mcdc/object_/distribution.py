import numpy as np

from collections.abc import Sequence
from numbers import Integral
from numpy import float64, int64
from numpy.typing import ArrayLike, NDArray

from mcdc.constant import (
    DISTRIBUTION_NONE,
    DISTRIBUTION_PMF,
    DISTRIBUTION_TABULATED,
    DISTRIBUTION_MULTITABLE,
    DISTRIBUTION_LEVEL_SCATTERING,
    DISTRIBUTION_EVAPORATION,
    DISTRIBUTION_MAXWELLIAN,
    DISTRIBUTION_KALBACH_MANN,
    DISTRIBUTION_TABULATED_ENERGY_ANGLE,
    DISTRIBUTION_N_BODY,
    INTERPOLATION_HISTOGRAM,
    INTERPOLATION_LINEAR,
)
from mcdc.object_.base import MCDCPolymorphic
from mcdc.object_.data import DataBase, DataTable
from mcdc.object_.util import (
    cdf_from_pdf,
    pdf_from_cdf,
    multi_cdf_from_pdf,
    cmf_from_pmf,
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
# Probability mass function
# ======================================================================================


class DistributionPMF(DistributionBase):
    """Discrete probability mass function.

    Parameters
    ----------
    value : array_like
        Values that may be sampled.
    pmf : array_like
        Nonnegative relative masses, normalized internally.
    """

    # MC/DC framework metadata
    label = "pmf_distribution"
    sub_type = DISTRIBUTION_PMF

    value: NDArray[float64]
    pmf: NDArray[float64]
    cmf: NDArray[float64]

    def __init__(self, value: ArrayLike, pmf: ArrayLike) -> None:
        super().__init__()

        self.value = np.asarray(value, dtype=float64)
        pmf_array = np.asarray(pmf, dtype=float64)

        if self.value.ndim != 1:
            print_error("value must be one-dimensional.")
        if pmf_array.ndim != 1:
            print_error("pmf must be one-dimensional.")
        if len(self.value) == 0:
            print_error("value and pmf must contain at least one entry.")
        if len(self.value) != len(pmf_array):
            print_error("value and pmf must have the same length.")

        self.pmf, self.cmf = cmf_from_pmf(pmf_array)

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - value {print_1d_array(self.value)}\n"
        text += f"  - pmf {print_1d_array(self.pmf)}\n"
        return text


# ======================================================================================
# Tabulated
# ======================================================================================


class DistributionTabulated(DistributionBase):
    """One-dimensional continuous tabulated distribution.

    Exactly one of ``pdf`` or ``cdf`` must be supplied alongside ``value``.
    PDF input uses the supplied interpolation regions, defaulting to linear.
    CDF input produces a histogram density. The distribution is normalized
    internally by integrating each segment under its declared law.
    """

    # MC/DC framework metadata
    label = "tabulated_distribution"
    sub_type = DISTRIBUTION_TABULATED

    pdf: DataTable

    def __init__(
        self,
        value: ArrayLike,
        pdf: ArrayLike | None = None,
        cdf: ArrayLike | None = None,
        interpolations=INTERPOLATION_LINEAR,
        interpolation_boundaries=None,
    ) -> None:
        super().__init__()

        if (pdf is None) == (cdf is None):
            print_error("Exactly one of pdf or cdf must be provided.")

        value_array = np.asarray(value, dtype=float64)

        if value_array.ndim != 1:
            print_error("value must be one-dimensional.")
        if len(value_array) == 0:
            print_error("value must contain at least one entry.")

        if pdf is not None:
            pdf_array = np.asarray(pdf, dtype=float64)

            if pdf_array.ndim != 1:
                print_error("pdf must be one-dimensional.")
            if len(pdf_array) != len(value_array):
                print_error("value and pdf must have the same length.")

            pdf_normalized, cdf_normalized = cdf_from_pdf(
                value_array, pdf_array, interpolations, interpolation_boundaries
            )
            interpolation = interpolations
        else:
            cdf_array = np.asarray(cdf, dtype=float64)

            if cdf_array.ndim != 1:
                print_error("cdf must be one-dimensional.")
            if len(cdf_array) != len(value_array):
                print_error("value and cdf must have the same length.")

            interpolation = INTERPOLATION_HISTOGRAM
            pdf_normalized, cdf_normalized = pdf_from_cdf(
                value_array,
                cdf_array,
            )

        self.pdf = DataTable(
            value_array,
            pdf_normalized,
            interpolation,
            interpolation_boundaries,
            aux=cdf_normalized,
        )

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - value {print_1d_array(self.pdf.x)}\n"
        text += f"  - probability density {print_1d_array(self.pdf.y)}\n"
        text += "  - cumulative distribution " f"{print_1d_array(self.pdf.aux[0])}\n"
        return text


# ======================================================================================
# Multi-table
# ======================================================================================


class DistributionMultiTable(DistributionBase):
    """Family of tabulated distributions indexed by another grid.

    ``tables`` supplies one DistributionTabulated object per grid point.

    ``interpolations`` and ``interpolation_boundaries`` describe interpolation
    between neighboring tables on ``grid``. Boundaries are exclusive upper
    point indices, with the final boundary equal to ``len(grid)``. A scalar
    law applies to the full grid; supported laws are histogram and linear.
    Each child table separately owns its PDF interpolation rules.
    """

    # MC/DC framework metadata
    label = "multi_table_distribution"
    sub_type = DISTRIBUTION_MULTITABLE

    grid: NDArray[float64]
    interpolations: NDArray[int64]
    interpolation_boundaries: NDArray[int64]
    tables: list[DistributionTabulated]

    def __init__(
        self,
        grid: ArrayLike,
        tables: list[DistributionTabulated],
        *,
        interpolations: int | Sequence[int] = INTERPOLATION_LINEAR,
        interpolation_boundaries: Sequence[int] | None = None,
    ) -> None:
        super().__init__()

        self.grid = np.asarray(grid, dtype=float64)
        if self.grid.ndim != 1 or len(self.grid) == 0:
            raise ValueError(
                "The multi-table grid must be a nonempty one-dimensional array"
            )
        if not np.all(np.isfinite(self.grid)) or np.any(np.diff(self.grid) < 0):
            raise ValueError("Grid values must be finite and nondecreasing")

        # Interpolation between neighboring grid tables.
        if isinstance(interpolations, Integral):
            self.interpolations = np.array([interpolations], dtype=int64)
            self.interpolation_boundaries = np.array([len(self.grid)], dtype=int64)
        else:
            self.interpolations = np.asarray(interpolations, dtype=int64)
            if interpolation_boundaries is None:
                raise ValueError(
                    "Interpolation boundaries are required for a sequence of laws"
                )
            self.interpolation_boundaries = np.asarray(
                interpolation_boundaries, dtype=int64
            )
        if (
            self.interpolations.ndim != 1
            or self.interpolation_boundaries.ndim != 1
            or len(self.interpolations) == 0
            or len(self.interpolations) != len(self.interpolation_boundaries)
        ):
            raise ValueError(
                "Interpolation laws and boundaries must be nonempty vectors of equal length"
            )
        if (
            self.interpolation_boundaries[0] < min(2, len(self.grid))
            or self.interpolation_boundaries[-1] != len(self.grid)
            or np.any(np.diff(self.interpolation_boundaries) <= 0)
        ):
            raise ValueError("Invalid grid interpolation boundaries")
        if np.any(
            ~np.isin(
                self.interpolations, [INTERPOLATION_HISTOGRAM, INTERPOLATION_LINEAR]
            )
        ):
            raise ValueError(
                "Multi-table interpolation requires histogram or linear interpolation"
            )

        # PDF interpolation belongs to the individual child distributions.
        if len(tables) != len(self.grid):
            raise ValueError(
                "Grid points and distribution tables must have equal lengths"
            )
        self.tables = tables

    def __repr__(self) -> str:
        text = super().__repr__()

        text += f"  - grid: {print_1d_array(self.grid)}\n"
        text += f"  - tables: {len(self.tables)}\n"

        for i, table in enumerate(self.tables):
            text += f"    - table[{i}] " f"(grid={self.grid[i]:.6g}, N={table.pdf.N})\n"

        return text


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
# Tabulated energy-angle
# ======================================================================================


class DistributionTabulatedEnergyAngle(DistributionBase):
    """Correlated tabulated outgoing energy and scattering-angle distribution.

    ``offset`` delimits outgoing-energy tables and ``cosine_offset`` delimits
    conditional cosine tables in the flattened arrays.
    """

    # MC/DC framework metadata
    label = "tabulated_energy_angle_distribution"
    sub_type = DISTRIBUTION_TABULATED_ENERGY_ANGLE

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
