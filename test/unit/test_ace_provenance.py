"""ACE library identity and temperature metadata survive HDF5 preparation."""

from pathlib import Path
import sys
from types import SimpleNamespace

import h5py
import pytest

pytest.importorskip("ACEtk")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from data_library_generator.common.ace_converter import (
    BOLTZMANN_K,
    read_source_libraries,
    read_target,
)
from data_library_generator.common.hdf5_writer import stamp


def source(library, *, label_field="title", suffix="10c", kT=2.53e-8):
    header = SimpleNamespace(
        zaid=f"11022.{suffix}",
        title="Na22 Ref. source-reference",
        comments="",
        temperature=kT,
        version="2.0.0",
        date="10/10/24",
    )
    setattr(header, label_field, f"{getattr(header, label_field)} {library}".strip())
    return SimpleNamespace(
        header=header, atom_number=11, mass_number=22, isomeric_state=0
    )


@pytest.mark.parametrize(
    "library", ["Lib80", "Lib80x", "Lib81", "Lib81x", "Lib100", "Lib100x"]
)
@pytest.mark.parametrize("label_field", ["title", "comments"])
@pytest.mark.parametrize("with_temperature", [False, True])
def test_library_identity_and_metadata_are_preserved(
    tmp_path, library, label_field, with_temperature
):
    ace = source(library, label_field=label_field)
    target = read_target(ace)
    sources = read_source_libraries(ace, target, with_temperature=with_temperature)
    with h5py.File(tmp_path / "provenance.h5", "w") as file:
        stamp(file, sources)
        assert list(file["source-libraries"]) == [library]
        record = file[f"source-libraries/{library}"]
        assert record["title"].asstr()[()] == ace.header.title
        assert record["comments"].asstr()[()] == ace.header.comments
        assert record["format"].asstr()[()] == "ACE"
        assert record["processing-date"].asstr()[()] == ace.header.date
        assert ("temperature" in record) == with_temperature
        assert ("kT" in record) == with_temperature
        if with_temperature:
            assert record["temperature"][()] == target.temperature
            assert record["temperature"].attrs["unit"] == "K"
            assert record["kT"][()] == pytest.approx(0.0253)
            assert record["kT"].attrs["unit"] == "eV"


@pytest.mark.parametrize("library", ["Lib80", "Lib80x", "Lib100", "Lib100x"])
def test_other_lanl_releases_use_header_temperature(library):
    ace = source(library)
    target = read_target(ace)
    assert list(read_source_libraries(ace)) == [library]
    assert target.temperature == ace.header.temperature * 1e6 / BOLTZMANN_K
    assert target.temperature != 293.6


@pytest.mark.parametrize("library", ["Lib81", "Lib81x"])
@pytest.mark.parametrize(
    "suffix,temperature,kT",
    [
        ("10c", 293.6, 2.53e-8),
        ("11c", 600.0, 5.1704e-8),
        ("12c", 900.0, 7.7556e-8),
        ("13c", 1200.0, 1.0341e-7),
        ("14c", 2500.0, 2.1543e-7),
        ("15c", 0.1, 8.6173e-12),
        ("16c", 233.15, 2.0091e-8),
        ("17c", 273.15, 2.3538e-8),
    ],
)
def test_library_suffix_selects_nominal_temperature(library, suffix, temperature, kT):
    ace = source(library, suffix=suffix, kT=kT)
    target = read_target(ace)
    assert target.temperature == temperature
    assert target.filename == f"neutron-Na22-{temperature}K.h5"
    metadata = read_source_libraries(ace, target, with_temperature=True)[library]
    assert metadata["temperature"] == temperature
    assert metadata["kT"] == kT * 1e6


@pytest.mark.parametrize(
    "library",
    ["", "OtherLibrary", "Lib81-extra", "Lib81x-extra", "Lib", "Lib81xx", "MyLib81"],
)
@pytest.mark.parametrize("suffix", ["10c", "99c"])
def test_unrecognized_library_uses_header_temperature_and_retains_metadata(
    library, suffix
):
    ace = source(library, suffix=suffix)
    target = read_target(ace)
    assert target.temperature == ace.header.temperature * 1e6 / BOLTZMANN_K
    assert target.temperature != 293.6
    metadata = read_source_libraries(ace, target, with_temperature=True)["unknown"]
    assert metadata["title"] == ace.header.title
    assert metadata["temperature"] == target.temperature
    assert metadata["kT"] == pytest.approx(0.0253)


@pytest.mark.parametrize("library", ["Lib81", "Lib81x"])
def test_known_library_with_unmapped_suffix_uses_header_temperature(library):
    ace = source(library, suffix="99c")
    assert read_target(ace).temperature == ace.header.temperature * 1e6 / BOLTZMANN_K


@pytest.mark.parametrize("library", ["", "OtherLibrary", "Lib81", "Lib81x"])
def test_source_library_override_preserves_header_and_temperature(library):
    ace = source(library)
    target = read_target(ace)
    sources = read_source_libraries(
        ace, target, with_temperature=True, source_library="CustomLibrary"
    )
    assert list(sources) == ["CustomLibrary"]
    assert sources["CustomLibrary"]["title"] == ace.header.title
    assert sources["CustomLibrary"]["temperature"] == target.temperature
    assert read_target(ace) == target


@pytest.mark.parametrize("library", ["TENDL-2021", "TENDL 2021"])
def test_tendl_provenance_uses_shared_metadata_reader(library):
    from data_library_generator.proton.generate import source_libraries

    ace = source(library)
    metadata = read_source_libraries(ace)["TENDL-2021"]
    assert metadata["version"] == "2021"
    assert metadata["title"] == ace.header.title
    proton_sources = source_libraries(ace, True)
    assert proton_sources["TENDL-2021"] == metadata
    assert proton_sources["PSTAR"]["format"] == "PSTAR"


@pytest.mark.parametrize("library", ["Lib80", "Lib80x", "Lib81", "Lib81x"])
def test_lanl_header_with_underlying_tendl_reference_keeps_lanl_identity(library):
    metadata = read_source_libraries(source(f"{library} TENDL-2021"))[library]
    assert metadata["version"] == ""


def test_proton_provenance_accepts_an_unrecognized_ace_library():
    from data_library_generator.proton.generate import source_libraries

    ace = source("OtherLibrary")
    assert list(source_libraries(ace, False)) == ["unknown"]
    assert list(source_libraries(ace, False, source_library="MyEvaluation")) == [
        "MyEvaluation"
    ]


def test_ace_label_cannot_overwrite_pstar_provenance():
    from data_library_generator.proton.generate import source_libraries

    with pytest.raises(ValueError, match="conflicts with PSTAR provenance"):
        source_libraries(source("OtherLibrary"), True, source_library="PSTAR")


@pytest.mark.parametrize("label", ["", " ", ".", "..", "not/a/group", 42])
def test_invalid_source_library_override_is_rejected(label):
    with pytest.raises(ValueError, match="Invalid source-library identifier"):
        read_source_libraries(source("OtherLibrary"), source_library=label)


@pytest.mark.parametrize("kT", [-1.0, float("nan"), float("inf")])
def test_invalid_header_thermal_energy_is_rejected(kT):
    with pytest.raises(
        ValueError, match="ACE header kT must be finite and nonnegative"
    ):
        read_target(source("OtherLibrary", kT=kT))
