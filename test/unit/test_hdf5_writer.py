"""Packed PDF tables preserve local interpolation regions and physical units."""

from pathlib import Path
import sys

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from data_library_generator.common import hdf5_writer as writer


def test_packed_tables_round_trip_regions_and_discontinuities(tmp_path):
    tables = [
        ([0.0, 1.0, 1.0, 2.0], [0.2, 0.2, 0.8, 0.8], ["histogram", "linear"], [2, 4]),
        ([1.0, 3.0], [0.25, 0.75], ["log"], [2]),
        ([0.0, 0.5, 1.0], [0.0, 2.0, 0.0], None, None),
    ]
    path = tmp_path / "packed.h5"
    with h5py.File(path, "w") as file:
        writer.multi_table(
            file,
            [1.0, 1.0, 10.0],
            tables,
            laws=["histogram", "linear"],
            boundaries=[2, 3],
        )
    with h5py.File(path) as file:
        assert file.attrs["type"] == "multi-table"
        assert file["grid"].attrs["unit"] == "eV"
        np.testing.assert_array_equal(file["grid"], [1, 1, 10])
        np.testing.assert_array_equal(file["interpolation/boundaries"], [2, 3])
        packed = file["tables"]
        np.testing.assert_array_equal(packed["offsets"], [0, 4, 6, 9])
        np.testing.assert_array_equal(packed["interpolation/offsets"], [0, 2, 3, 4])
        assert packed["value"].attrs["unit"] == "eV"
        assert packed["pdf"].attrs["unit"] == "1/eV"
        assert "unit" not in packed["offsets"].attrs
        for i, (value, pdf, laws, boundaries) in enumerate(tables):
            start, end = packed["offsets"][i : i + 2]
            np.testing.assert_array_equal(packed["value"][start:end], value)
            np.testing.assert_array_equal(packed["pdf"][start:end], pdf)
            start, end = packed["interpolation/offsets"][i : i + 2]
            assert packed["interpolation/laws"].asstr()[start:end].tolist() == (
                laws or ["linear"]
            )
            np.testing.assert_array_equal(
                packed["interpolation/boundaries"][start:end],
                boundaries if boundaries is not None else [len(value)],
            )


def test_dimensionless_grid_and_cosines_have_no_unit(tmp_path):
    with h5py.File(tmp_path / "cosines.h5", "w") as file:
        writer.multi_table(
            file,
            [-1, 1],
            [([-1, 1], [0.5, 0.5], None, None)] * 2,
            grid_unit=None,
            unit=None,
        )
        for name in ("grid", "tables/value", "tables/pdf"):
            assert "unit" not in file[name].attrs
        writer.pdf_table(file.create_group("single"), [-1, 1], [0.5, 0.5], unit=None)
        assert file["single"].attrs["type"] == "table"


@pytest.mark.parametrize(
    "tables",
    [
        [],
        [([0, 1], [1], None, None)],
        [([0], [1], None, None)],
        [([0, 1, 2], [1, 1, 1], ["linear", "log"], [2, 1])],
    ],
)
def test_invalid_packed_tables_do_not_write_partial_datasets(tmp_path, tables):
    with h5py.File(tmp_path / "invalid.h5", "w") as file:
        with pytest.raises(ValueError):
            writer.pdf_tables(file, tables)
        assert not len(file)


@pytest.mark.parametrize("grid", [[1], [2, 1], [1, np.inf]])
def test_multi_table_rejects_invalid_indexing_grid(tmp_path, grid):
    with h5py.File(tmp_path / "invalid-grid.h5", "w") as file:
        with pytest.raises(ValueError):
            writer.multi_table(file, grid, [([0, 1], [1, 1], None, None)] * 2)


def test_exclusive_regions_cover_each_interval_once(tmp_path):
    with h5py.File(tmp_path / "regions.h5", "w") as file:
        writer.interpolation(file, 6, ["linear", "log"], [3, 6])
        laws = file["interpolation/laws"].asstr()[:]
        boundaries = file["interpolation/boundaries"][:]
        # An interval belongs to the first region whose exclusive end exceeds
        # its upper-point index; the regions share their boundary coordinate.
        regions = np.searchsorted(boundaries, np.arange(1, 6), side="right")
        assert laws[regions].tolist() == ["linear", "linear", "log", "log", "log"]
        single = file.create_group("single")
        writer.interpolation(single, 6)
        np.testing.assert_array_equal(single["interpolation/boundaries"], [6])


@pytest.mark.parametrize("boundaries", [[1, 6], [3, 5], [3, 3, 6], [3.0, 6.0], []])
def test_invalid_region_boundaries_are_rejected_before_writing(tmp_path, boundaries):
    with h5py.File(tmp_path / "invalid-regions.h5", "w") as file:
        with pytest.raises(ValueError):
            writer.interpolation(file, 6, ["linear"] * len(boundaries), boundaries)
        assert not len(file)


@pytest.mark.parametrize(
    "interpolants,boundaries,expected",
    [([2, 5], [3, 6], [3, 6]), ([], [], [6])],
)
def test_ace_regions_retain_source_point_counts(
    tmp_path, interpolants, boundaries, expected
):
    pytest.importorskip("ACEtk")
    from types import SimpleNamespace
    from data_library_generator.common.ace_converter import write_interpolation

    source = SimpleNamespace(interpolants=interpolants, boundaries=boundaries)
    with h5py.File(tmp_path / "ace-regions.h5", "w") as file:
        write_interpolation(file, source, 6)
        np.testing.assert_array_equal(file["interpolation/boundaries"], expected)
