"""Write MC/DC HDF5 representations and coordinate output-file handling.

Source-format converters supply already interpreted values in the output
units. This module handles datasets, source-library provenance, interpolation
metadata, format identification, companion-file merging, and staged publication
of temperature files. It does not read ACE records or convert physical units.

``FORMAT`` and ``VERSION`` identify the output specification independently of
the MC/DC software version and source-library release. Compatibility checks
require exact matches. ``library_exists`` coordinates the three
neutron output filenames; the lower-level writers operate on caller-provided
files and groups.
"""

from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
import os
import tempfile

import h5py
import numpy as np

# ======================================================================================
# Format identifiers and source-library metadata schema
# ======================================================================================


FORMAT = "mcdc-data"
VERSION = "0.1.0"
INTERPOLATION_LAWS = {"histogram", "linear", "semilog-x", "semilog-y", "log"}
SOURCE_FIELDS = (
    "title",
    "version",
    "release-date",
    "reference",
    "comments",
    "format",
    "format-version",
    "processing-date",
)
SOURCE_UNITS = {"temperature": "K", "kT": "eV"}


# ======================================================================================
# File identity and source-library provenance
# ======================================================================================


def stamp(file, sources):
    """Write format identifiers and merge the supplied source-library records.

    Parameters
    ----------
    file : h5py.File
        Writable destination. When merging an existing nonempty file, the
        caller checks its format first; this function writes the identifiers
        but does not check the previous version.
    sources : dict of str to dict
        Library subgroup names mapped to scalar metadata. ``SOURCE_FIELDS``
        identifies the string fields; missing strings become empty. Optional
        numeric ``temperature`` and ``kT`` values are already in K and eV.
        The caller includes those fields only for temperature-specific files.

    Raises
    ------
    ValueError
        No sources are supplied; a library name is empty, ``.`` or ``..``,
        or contains a slash; a string field has a non-string value; or
        incoming nonempty metadata conflicts with an existing source record.

    Notes
    -----
    Root attributes ``format`` and ``format-version`` use the module constants.
    Each supplied library is stored under ``source-libraries/{library}`` as
    scalar datasets. Existing nonempty strings survive missing incoming
    information. Conflicting nonempty strings raise an error; numeric source
    temperatures are compared with NumPy's default ``isclose`` tolerances.
    Other source libraries and unrelated file contents remain intact.

    A UTC ISO 8601 ``conversion-time`` is generated for this write and applied
    only to the supplied libraries. Calling this function updates that time
    even if all other values already match, so skipped conversions must avoid
    calling it. Metadata conflicts are checked before writing; HDF5 I/O
    failures during publication can still leave partial changes in an open
    file. Existing named metadata datasets are recreated when written.
    """
    # Validate and merge source records before modifying the HDF5 file
    if not sources:
        raise ValueError("At least one source library is required")
    updates = {}
    for library, source in sources.items():
        if not library or library in (".", "..") or "/" in library:
            raise ValueError(f"Invalid source-library identifier: {library}")
        group = file.get(f"source-libraries/{library}")
        values = {}
        for name in SOURCE_FIELDS:
            value = source.get(name, "")
            if not isinstance(value, str):
                raise ValueError(f"{library}/{name}: source metadata must be a string")
            previous = (
                group[name].asstr()[()] if group is not None and name in group else ""
            )
            if previous and value and previous != value:
                raise ValueError(f"{file.filename}: conflicting {name} for {library}")
            values[name] = previous or value
        # Compare temperature-specific quantities in their declared output units
        for name in SOURCE_UNITS:
            if name in source:
                value = source[name]
                if (
                    group is not None
                    and name in group
                    and not np.isclose(group[name][()], value)
                ):
                    raise ValueError(
                        f"{file.filename}: conflicting {name} for {library}"
                    )
                values[name] = value
        updates[library] = values
    # Publish format identifiers and stamp the supplied library contributions
    file.attrs["format"] = FORMAT
    file.attrs["format-version"] = VERSION
    conversion_time = (
        datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")
    )
    for library, values in updates.items():
        group = file.require_group(f"source-libraries/{library}")
        values["conversion-time"] = conversion_time
        for name, value in values.items():
            if name in group:
                del group[name]
            quantity(group, name, value, SOURCE_UNITS.get(name))


# ======================================================================================
# Datasets and reusable tabulated representations
# ======================================================================================


def quantity(group, name, data, unit=None):
    """Create a scalar or array dataset and attach its supplied unit metadata.

    Parameters
    ----------
    group : h5py.Group or h5py.File
        Writable parent object.
    name : str
        Dataset name or relative path accepted by ``create_dataset``.
    data : object
        Scalar or array-like data accepted by h5py. Shape and datatype are
        inferred from this object.
    unit : str or sequence of str, optional
        Unit attribute to attach. A sequence is used for polynomial
        coefficients with different units. ``None`` omits the attribute,
        indicating dimensionless data.

    Returns
    -------
    h5py.Dataset
        Newly created dataset, available for further metadata assignment.

    Notes
    -----
    Values and unit strings are stored as supplied; no unit conversion or
    dimensional validation is performed. An existing object at ``name`` is
    not overwritten, and h5py creation errors propagate to the caller.
    """
    dataset = group.create_dataset(name, data=data)
    if unit is not None:
        dataset.attrs["unit"] = unit
    return dataset


def interpolation(group, count, laws=None, boundaries=None):
    """Create interpolation metadata for a table of known length.

    Parameters
    ----------
    group : h5py.Group
        Parent representation in which to create an ``interpolation`` group.
    count : int
        Number of local tabulated points, at least two.
    laws : sequence of str, optional
        One name from ``INTERPOLATION_LAWS`` per region. The default is one
        ``linear`` region.
    boundaries : array_like of int, optional
        Zero-based upper-point index of each region. The default is a single
        endpoint, ``count - 1``. The first region starts at point zero;
        subsequent regions start at the previous region's upper point.

    Raises
    ------
    ValueError
        There are fewer than two points; boundaries are empty, noninteger,
        or not a one-dimensional strictly increasing sequence beginning
        above zero and ending at ``count - 1``; counts of laws and boundaries
        differ; or a law name is unsupported.

    Notes
    -----
    For six points, boundaries ``[2, 5]`` describe intervals 0--1 and 1--2
    in the first region, followed by 2--3, 3--4, and 4--5 in the second.
    Indices are local to the described table, including a cross-section
    tail stored with a nonzero shared-grid offset.

    The child datasets are ``laws`` (fixed-width ``S16`` strings) and
    ``boundaries`` (integers). The function validates this metadata before
    creating the group, but does not inspect coordinates or values, evaluate
    interpolation, or check positivity domains for logarithmic laws.
    """
    # Validate that the interpolation regions cover every local table interval
    if count < 2:
        raise ValueError(f"{group.name}: tabulated data requires at least two points")
    if laws is None:
        laws = ["linear"]
    boundaries = np.asarray([count - 1] if boundaries is None else boundaries)
    if (
        boundaries.ndim != 1
        or len(boundaries) == 0
        or not np.issubdtype(boundaries.dtype, np.integer)
        or len(laws) != len(boundaries)
        or boundaries[-1] != count - 1
        or np.any(np.diff(np.r_[0, boundaries]) <= 0)
        or not set(laws) <= INTERPOLATION_LAWS
    ):
        raise ValueError(f"{group.name}: invalid interpolation regions")
    # Store explicit interpolation laws and their local upper-point indices
    rule = group.create_group("interpolation")
    rule.create_dataset("laws", data=np.asarray(laws, dtype="S16"))
    rule.create_dataset("boundaries", data=boundaries)


def pdf_table(group, values, pdf, law="linear", unit="eV"):
    """Write a one-dimensional tabulated probability-density representation.

    Parameters
    ----------
    group : h5py.Group
        Empty destination Distribution group.
    values : array_like
        Sample-coordinate grid, already expressed in the output unit.
    pdf : array_like
        Probability-density values aligned with the coordinate grid, already
        expressed in the inverse output unit.
    law : str, optional
        Interpolation law for all intervals, defaulting to ``linear``.
    unit : str, optional
        Coordinate unit, defaulting to ``eV``; the density unit becomes
        ``1/{unit}``. Pass ``None`` for dimensionless cosine coordinates
        and densities, so both datasets omit the unit attribute.

    Raises
    ------
    ValueError
        Coordinate and density lengths differ, there are fewer than two
        points, or the interpolation law is unsupported.

    Notes
    -----
    Writes ``type = tabulated``, ``value``, ``pdf``, and a single-region
    ``interpolation`` group. The caller supplies correctly converted,
    normalized data: this function does not normalize the PDF, construct
    a CDF, or validate monotonicity, nonnegativity, or total probability.
    A failure after dataset creation can leave a partial destination group.
    """
    if len(values) != len(pdf):
        raise ValueError(
            f"{group.name}: PDF coordinates and densities differ in length"
        )
    group.attrs["type"] = "tabulated"
    quantity(group, "value", values, unit)
    quantity(group, "pdf", pdf, f"1/{unit}" if unit is not None else None)
    interpolation(group, len(values), [law])


# ======================================================================================
# Format compatibility and prepared-file reuse
# ======================================================================================


def check_format(file):
    """Require exact supported format identifiers on an open HDF5 file.

    Parameters
    ----------
    file : h5py.File
        Open file whose root attributes are checked without modification.

    Raises
    ------
    ValueError
        Either ``format`` or ``format-version`` is missing or differs from
        ``FORMAT`` or ``VERSION`` respectively.

    Notes
    -----
    This is an identification check, not a schema or physical-data validator.
    It does not accept a version merely because its major number matches,
    and does not migrate existing files.
    """
    if (
        file.attrs.get("format") != FORMAT
        or file.attrs.get("format-version") != VERSION
    ):
        raise ValueError(f"{file.filename}: incompatible format or format version")


def library_exists(directory, filename, nuclide, symbol, sources_by_file):
    """Check whether the three neutron output files can be reused as prepared.

    Parameters
    ----------
    directory : str or pathlib.Path
        Output library directory, read without modification.
    filename : str
        Nuclide-temperature basename, such as ``H2-293.6K.h5``.
    nuclide : str
        Nuclide identifier selecting the companion ``{nuclide}.h5``.
    symbol : str
        Element symbol selecting the companion ``{symbol}.h5``.
    sources_by_file : dict of str to dict
        Expected source-library records keyed by each of the three basenames.
        Each value has the ``sources`` structure accepted by ``stamp``.
        Per-file records allow an element title to be empty while nuclide
        titles and temperature-only provenance are retained where applicable.

    Returns
    -------
    bool
        True if every file has the required identity/temperature paths,
        expected source groups, scalar string metadata, and a parseable UTC
        conversion timestamp ending in ``Z``. Supplied source temperature
        fields must also have matching scalar values and unit attributes.
        Missing or incomplete output returns False so conversion can proceed.

    Raises
    ------
    ValueError
        An existing file has incompatible format identifiers, as checked by
        ``check_format``; this is an error rather than a reusable-file miss.
    OSError
        An existing HDF5 file cannot be opened or read.
    KeyError
        An expected filename is absent from ``sources_by_file``.

    Notes
    -----
    Nonempty expected strings must match; empty expectations permit known
    metadata already recorded by another contribution. Additional source
    libraries and projectile groups are allowed. Numeric source temperatures
    use NumPy's default ``isclose`` tolerances. The timestamp is checked for
    valid storage, not compared with the current time.

    This is an incremental-conversion completeness check. It does not
    traverse and validate every reaction, compare cross-section sums, or
    confirm every identity value against an ACE table. No timestamps or file
    contents are updated by the check.
    """
    # Define the minimum required contents of the three coordinated output files
    required = (
        (filename, ("temperature", "kT", "projectiles/neutron/reactions")),
        (f"{symbol}.h5", ("atomic-number",)),
        (
            f"{nuclide}.h5",
            ("mass-number", "atomic-weight-ratio", "excitation-level", "fissionable"),
        ),
    )
    complete = True
    for filename, fields in required:
        path = Path(directory) / filename
        if not path.is_file():
            complete = False
            continue
        with h5py.File(path, "r") as file:
            # Check file identity and required paths before examining provenance
            check_format(file)
            complete = complete and all(field in file for field in fields)
            for library, source in sources_by_file[filename].items():
                group = file.get(f"source-libraries/{library}")
                if group is None:
                    complete = False
                    continue
                # Check scalar source metadata and the stored UTC conversion time
                for name in (*SOURCE_FIELDS, "conversion-time"):
                    if (
                        name not in group
                        or group[name].shape != ()
                        or h5py.check_string_dtype(group[name].dtype) is None
                    ):
                        complete = False
                        continue
                    value = group[name].asstr()[()]
                    if name == "conversion-time":
                        try:
                            timestamp = datetime.fromisoformat(
                                value.replace("Z", "+00:00")
                            )
                            complete = (
                                complete
                                and value.endswith("Z")
                                and timestamp.tzinfo is not None
                            )
                        except ValueError:
                            complete = False
                    elif source.get(name) and value != source[name]:
                        complete = False
                # Match supplied source temperatures and their unit attributes
                for name, unit in SOURCE_UNITS.items():
                    if name in source:
                        complete = (
                            complete
                            and name in group
                            and group[name].shape == ()
                            and np.isclose(group[name][()], source[name])
                            and group[name].attrs.get("unit") == unit
                        )
    return complete


# ======================================================================================
# Companion-file merging and staged file publication
# ======================================================================================


def write_companion(path, fields, sources):
    """Create or merge an element/nuclide companion's identity and provenance.

    Parameters
    ----------
    path : str or pathlib.Path
        Companion filename opened in append mode; its parent directory must
        already exist.
    fields : dict of str to object
        Dimensionless numeric or boolean identity datasets to create if
        absent, such as atomic number or mass number and fissionable flag.
    sources : dict of str to dict
        Source-library records applicable to this file, passed to ``stamp``.

    Raises
    ------
    ValueError
        Existing nonempty content has incompatible format identifiers, an
        identity value conflicts, or source provenance conflicts.
    OSError
        Opening or modifying the file fails.

    Notes
    -----
    Existing identity values must agree under ``allclose`` with ``rtol=1e-10``
    and zero absolute tolerance; matching values are retained. Named identity
    datasets omit unit attributes, including when such an attribute was
    previously attached. Other projectile data and unrelated objects remain
    intact. Source records are merged and their conversion times updated.

    Identity comparisons occur before provenance or dataset writes. The
    subsequent append operation is not atomic, and an I/O failure can leave
    partial updates. A surrounding temperature-file ``atomic_file`` context
    does not roll back changes to companion files.
    """
    with h5py.File(path, "a") as file:
        # Validate existing file and target identities before applying updates
        if len(file) or len(file.attrs):
            check_format(file)
        for name, value in fields.items():
            if name in file and not np.allclose(
                file[name][()], value, rtol=1e-10, atol=0
            ):
                raise ValueError(f"{path}: inconsistent {name}")
        # Merge provenance and fill missing dimensionless identity datasets
        stamp(file, sources)
        for name, value in fields.items():
            if name not in file:
                quantity(file, name, value)
            if "unit" in file[name].attrs:
                del file[name].attrs["unit"]


@contextmanager
def atomic_file(path, sources):
    """Stage one HDF5 output and replace its destination on successful exit.

    Parameters
    ----------
    path : str or pathlib.Path
        Final output path. The parent directory must already exist.
    sources : dict of str to dict
        Source-library records written by ``stamp`` before control is yielded.

    Yields
    ------
    h5py.File
        Writable file at a temporary path in the destination directory,
        already carrying format identifiers and source-library groups.

    Raises
    ------
    ValueError
        Source metadata is invalid or conversion inside the context raises
        a ValueError; exceptions from the caller's body propagate unchanged.
    OSError
        Creating, closing, or replacing the HDF5 file fails.

    Notes
    -----
    The destination is untouched while the caller writes the staged file.
    After a normal context exit, the HDF5 handle is closed and ``os.replace``
    publishes the file. Staging in the same directory keeps replacement on
    the same filesystem. If conversion fails before replacement, an existing
    destination remains intact. The staged path is removed in either case.

    Publication covers this file only: other files modified within the
    context, including companion files, are outside this transaction.
    """
    # Stage beside the destination so replacement stays on the same filesystem
    path = Path(path)
    with tempfile.NamedTemporaryFile(
        dir=path.parent, suffix=".h5", delete=False
    ) as temp:
        staged = Path(temp.name)
    try:
        # Yield a stamped output file for the caller to populate
        with h5py.File(staged, "w") as file:
            stamp(file, sources)
            yield file
        # Publish only after writing and closing the staged file successfully
        os.replace(staged, path)
    finally:
        # Clean up the staged path on both success and failure
        staged.unlink(missing_ok=True)
