==================
Serial Performance
==================

The ``performance/serial`` suite measures Python and Numba execution using one
process per job. Each configured case-mode pair is an independent Maestro step,
with that mode's sampling points run sequentially.
Python and Numba jobs have no dependency on each other.
The supplied ``kobayashi-detector`` input defines its multigroup materials
without an external data directory.

Configuring tasks
-----------------

Each case in ``task.yaml`` has a ``numba`` block, a ``python`` block, or both.
Omit a mode to disable it. Each block defines its own logarithmic particle-count
bounds, task count, and walltime factor:

.. code-block:: yaml

   kobayashi-detector:
     numba:
       logN_min: 1
       logN_max: 6
       N_task: 11
       walltime_factor: 1.0
     python:
       logN_min: 1
       logN_max: 4
       N_task: 7
       walltime_factor: 1.0

The example selects eleven Numba and seven Python points on separate grids.
``input.py`` sets ``N_batch``; the processor reads actual history and batch counts
from each output.

``--walltime`` specifies a base duration in hours, multiplied by the mode's
``walltime_factor``, rounded up to scheduler resolution, and capped at the
platform maximum. If omitted, the platform maximum is the base.
Local execution ignores walltime.
An HPC job requests one process on one exclusive node, not every CPU core on that
node. Use :doc:`../campaigns/setup` for shared scheduler and account settings.

Launching and resuming
----------------------

From ``performance/serial`` in MC/DC-VVP, launch locally:

.. code-block:: console

   python launch.py

Or on a configured cluster:

.. code-block:: console

   python launch.py --platform dane --walltime 1.0

Full tally results, standard metadata, runtime components, and the ``performance/``
group are written to platform-specific output directories:

.. code-block:: text

   cases/kobayashi-detector/outputs/<platform>/output_numba_10000.h5
   cases/kobayashi-detector/outputs/<platform>/output_python_10000.h5
   maestro_run_<platform>_<timestamp>/
   study_<platform>.yaml

Completed mode/particle outputs for the selected platform are skipped on relaunch.
A fully complete case-mode pair is skipped independently of the other mode.
Outputs without full tally results are rejected; preserve them separately before
relaunching to regenerate measurements with the required content.

Processing and records
----------------------

From the suite directory, process the latest saved study for each platform:

.. code-block:: console

   python process.py

Alternatively, select one or more saved studies:

.. code-block:: console

   python process.py maestro_run_dane_<timestamp> maestro_run_local_<timestamp>

Processing uses saved ``task.yaml`` files rather than the currently edited one.
Older task files using shared bounds or ``python_mode`` must be converted to the
nested execution-mode schema before processing.
Results are written under ``results/<platform>/`` with a ``records.csv`` and
five log-log plots per case:

- ``runtime.png``: measured total runtime, including compilation.
- ``tracking_rate.png``: histories per second.
- ``precision.png``: inverse maximum relative variance in percent-squared units.
- ``precision_rate.png``: precision per history.
- ``fom.png``: precision per second.

The :doc:`metric definitions <index>` specify the fixed reference and handling
of invalid variances.
Numba uses a solid blue line with hollow circles; Python uses a dashed vermilion
line with crosses.

Modes are processed independently: a missing output in one mode does not discard
the other mode's point. CSV rows span the union of available particle counts,
with blank fields for missing mode measurements.
A mode with no available outputs has no curve or CSV columns.
At least one configured sampling point must be available.
Every output must have one MPI rank, and overlapping particle counts must have
matching total histories and batch counts.

The CSV records derived metrics, rank counts, fractional maximum relative
variance, reference mode and particle/history counts, and the number of
nonzero-reference bins. Raw ``performance/effective_variance`` is retained for
comparison but is not used for the plotted precision metrics.
All runtime-dependent metrics use the measured total runtime.

After jobs finish, ``python cleanup.py`` removes generated outputs, Maestro
records, processed results, and platform-specific study files.
The launcher regenerates the untracked study file for the selected platform.
