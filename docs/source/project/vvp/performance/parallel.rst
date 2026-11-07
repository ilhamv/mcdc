====================
Parallel Performance
====================

The ``performance/parallel`` suite compares full-node Numba runs across platforms,
node counts, workloads, and transport methods.
Each node-count/workload point is an independent job.
The resulting performance envelope provides context for the performance-level
and scaling-evidence summary.
See :doc:`index` for metric definitions and the distinction between timing,
statistical efficiency, and scaling efficiency.

Baseline and resource matrix
----------------------------

Each ``task.yaml`` entry identifies ``cases/<case>/input.py`` and baseline particles
per batch. For example:

.. code-block:: yaml

   kobayashi-detector:
     method: Analog
     N_particle_base: 30_000_000
     walltime: {1: 1.5, 2: 3.0, 4: 6.0, 8: 12.0, 16: 24.0}

Calibrate the baseline for approximately one hour on one full node, including
compilation and full tally output, while holding the batch count and physical
problem fixed.
The scalar baseline applies to every platform; it does not constitute an
independent GPU calibration.
For separate calibrations, use a mapping of positive integers:

.. code-block:: yaml

   kobayashi-detector:
     method: Analog
     N_particle_base:
       dane: 30_000_000
       tuolumne: 120_000_000
     walltime: {1: 1.5, 2: 3.0, 4: 6.0, 8: 12.0, 16: 24.0}

The mapping is an illustrative configuration, not a measured calibration.
Each selected platform must have an entry; missing baselines raise an error.

Node counts are powers of two from one through ``N_node_max`` (up to the configured
platform limit). Workload multipliers are 1, 2, 4, 8, and 16.
For node count :math:`P` and multiplier :math:`m`, particles per batch are
:math:`N_\mathrm{particle,base}Pm`.
Each case requires a ``walltime`` mapping with positive hours for every workload
multiplier. The selected value applies at every node count.
Document the calibration basis and chosen margin in task comments.
A launch ``--walltime`` overrides all matrix jobs uniformly, but every task value
must still be present and valid. Requests use exclusive nodes, are rounded up to
scheduler resolution, and are capped at the platform maximum.

CPU-only platforms use all CPU cores, one MPI rank per core.
GPU platforms use all GPUs, one rank per GPU, ``--target=gpu``, and
``--gpu_state_storage=united``; GPU jobs do not attempt to occupy every CPU core.
Compilation remains in total runtime and Numba caching is disabled.

Launching and output isolation
------------------------------

The repository-wide configuration uses a platforms list and maximum node count:

.. code-block:: python

   "performance/parallel": {
       "enabled": True,
       "platforms": ["dane", "tuolumne"],
       "N_node_max": 128,
       "walltime": None,
   },

On each target machine, run the root launcher with that machine's ``--platform``.
The list does not submit to another host remotely.
For a suite-level preview without job submission:

.. code-block:: console

   python launch.py --platform dane --N_node_max 128 --dry_run
   python launch.py --platform tuolumne --N_node_max 128 --dry_run

These are platform-specific alternatives. Omit ``--dry_run`` to submit on the
selected host. Add ``--walltime HOURS`` to override every job with the same
duration. Load the platform's scheduler/GPU environment first;
see :doc:`../campaigns/setup`.

Platform names distinguish output directories, filenames, task names, study
files, Maestro directories, and workspaces:

.. code-block:: text

   cases/<case>/outputs/tuolumne/output_tuolumne_n004_m16.h5
   workspaces/tuolumne/<task>_*/
   maestro_run_tuolumne_<timestamp>/
   study_tuolumne.yaml

The example output is a four-node, multiplier-16 task.
Each invocation, including a retry, receives a unique workspace and
``NUMBA_CACHE_DIR`` accessible to every rank and allocated node.
The input remains in its case directory; output is written directly to its
absolute platform-specific path. Workspaces remain available until cleanup.
For future inputs that need external data, resolve paths relative to ``__file__``
or use absolute paths because execution occurs in isolated workspaces.

Matching completed outputs are skipped.
An existing output with an incompatible baseline or rank count raises an error
and must be preserved separately before starting the changed campaign.
Filenames encode platform, nodes, and multiplier, not arbitrary campaign revisions.

Comparing transport methods
---------------------------

Inputs can share a comparison when their task entries have the same
``comparison`` label and distinct ``method`` labels.
For example, analog and weight-window inputs can use
``comparison: kobayashi-detector`` with ``method: Analog`` and
``method: Weight windows``.
Without labels, comparison defaults to the case name and method to ``Numba``.
Only group equivalent observables, tally grids, normalization conventions, and
batch counts.
A new performance case requires an input and task entry, not a case-specific
launcher or processor.

Processing platform studies
---------------------------

Bring outputs and saved Maestro configuration directories from different
filesystems together under one suite directory.
Then process the latest saved study for every available platform:

.. code-block:: console

   python process.py

Or select studies explicitly:

.. code-block:: console

   python process.py maestro_run_dane_<timestamp> maestro_run_tuolumne_<timestamp>

Selection uses saved configurations, not the current ``task.yaml``.
Missing points are reported and skipped; duplicate selected outputs are rejected.
Tally scores, shapes, and grids must match the comparison reference, batch counts
must agree, and rank counts must match the saved CPU/GPU resource layout.
Outputs without full tallies are rejected.
Legacy flat outputs such as ``output_n004_m16.h5`` can be read only through an
explicitly identifying legacy saved study; new launches do not reuse them.

Saved configurations are copied under ``results/studies/``.
Each ``results/<comparison>/`` contains:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - Content
   * - ``records.csv``
     - Measured metrics, platform/method, source paths, and reference provenance
   * - ``runtime.png``
     - Total runtime versus total histories
   * - ``tracking_rate.png`` and ``fom.png``
     - Tracking rate and FOM per node versus total histories
   * - ``precision.png`` and ``precision_rate.png``
     - Precision and precision per history
   * - ``performance_envelope.png``
     - Tracking rate and FOM per node, logarithmic histories and linear vertical axes
   * - ``scaling_summary.png``
     - Performance level with weak-scaling and workload-sensitivity wings
   * - ``scaling_summary.csv``
     - Summary targets, values, interpolation brackets, source paths, and reference

The measured CSV also retains whole-job rates, raw ``effective_variance``,
fractional maximum relative variance, contributing-bin counts, and runtime
components. The summary replaces separate strong- and weak-scaling efficiency
plots. Platform/method combinations use distinct colors, and detailed curves use
distinct line styles for node counts.

Main bar and scaling wings
--------------------------

For each platform/method group, the main bar M and right wing use the largest
available node count.
M is evaluated at the largest common histories-per-node workload covered by
**every included curve**, including the smaller-node curves in the left wing.
It can therefore be below a group's largest measured workload and need not be
the tallest bar.

- Moving left away from M reduces nodes at fixed histories per node, providing
  a direct weak-scaling comparison.
- Moving right reduces histories per node at fixed nodes, showing workload
  sensitivity and indirect strong-scaling evidence. It does not measure
  additional node-dependent communication costs.

All groups use the same rightward workload targets and reduction factor.
By default, the summary attempts four right-wing steps with a factor of four,
stopping earlier for every group if any point would require extrapolation.
For a compatible five-point factor-two matrix, use:

.. code-block:: console

   python process.py --workload_factor 2 --workload_steps 4

The detailed envelope retains all measured points, including workloads beyond
the shared summary range.
If main/left workloads have no common coverage, the summary figure explains the
missing overlap and no summary CSV is generated. Collect overlap points or select
compatible studies rather than extrapolating.

Exact measurements are preferred. Otherwise, total runtime and precision are
interpolated separately and linearly between adjacent measured history counts,
then used to derive rates and FOM.
Hatching marks interpolated bars; the summary CSV records endpoints and files.
Invalid precision endpoints are not bridged, so the corresponding FOM bar remains
unavailable.
Sparse interpolation is an approximation: inspect its brackets and collect more
measurements when a conclusion depends on a wide interval.
The envelope, rather than the summary alone, indicates whether M is near a plateau.

Cleanup
-------

After all jobs have stopped, ``python cleanup.py`` removes generated outputs for
all platforms, workspaces/caches, Maestro records, processed results, and generated
study files. Inputs and ``task.yaml`` are preserved.
