"""Merge particle-specific contributions into the common native library."""

import os
from pathlib import Path
import shutil
import tempfile

import h5py
import numpy as np

SHARED_PROPERTIES = (
    "nuclide_name",
    "excitation_level",
    "temperature",
    "atomic_number",
    "mass_number",
    "atomic_weight_ratio",
    "radiation_length",
)


def has_data(path, group):
    """Check for a particle contribution rather than just its nuclide file."""
    if not Path(path).exists():
        return False
    with h5py.File(path, "r") as file:
        return group in file


def write_contribution(source, path, particle):
    """Validate and atomically merge one generated particle contribution.

    Existing shared properties are retained. Floating-point properties must
    agree within a relative tolerance of 1e-5; identity fields must match exactly.
    The legacy root fissionable flag describes neutron-induced fission.
    """
    if particle not in ("neutron", "proton"):
        raise ValueError(f"Unsupported library contribution: {particle}")
    path = Path(path)
    owned = [f"{particle}_reactions"]
    if particle == "proton":
        owned.append("stopping_power")

    # Update a sibling temporary file so a failed merge preserves the library.
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, suffix=".h5")
    os.close(descriptor)
    try:
        if path.exists():
            shutil.copyfile(path, temporary)
        with h5py.File(temporary, "a" if path.exists() else "w") as target:
            for name in SHARED_PROPERTIES:
                if name not in source:
                    continue
                if name not in target:
                    source.copy(name, target)
                    continue
                old, new = target[name][()], source[name][()]
                if name in ("temperature", "atomic_weight_ratio", "radiation_length"):
                    matches = np.isclose(old, new, rtol=1e-5, atol=1e-12)
                else:
                    matches = old == new
                if not matches:
                    raise ValueError(
                        f"{path.name}: conflicting {name}: existing {old!r}, "
                        f"{particle} contribution {new!r}"
                    )

            # Preserve legacy ACE provenance under its original particle owner.
            provenance = target.require_group("provenance")
            for species in ("neutron", "proton"):
                if f"{species}_reactions" in target and species not in provenance:
                    record = provenance.create_group(species)
                    for name, value in target.attrs.items():
                        if name.startswith("source_"):
                            record.attrs[name] = value

            # Replace only groups supplied by this conversion.
            for name in owned:
                if name in source:
                    if name in target:
                        del target[name]
                    source.copy(name, target)

            # Keep ACE and stopping-power-only provenance distinct.
            key = (
                particle
                if f"{particle}_reactions" in source
                else "proton_stopping_power"
            )
            if key in provenance:
                del provenance[key]
            record = provenance.create_group(key)
            for name, value in source.attrs.items():
                record.attrs[name] = value
            record.attrs["fissionable"] = bool(source["fissionable"][()])

            if particle == "neutron":
                if "fissionable" in target:
                    del target["fissionable"]
                source.copy("fissionable", target)
            elif "neutron_reactions" not in target:
                if "fissionable" in target:
                    del target["fissionable"]
                target.create_dataset("fissionable", data=False)

        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)
