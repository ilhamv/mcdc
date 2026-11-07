.. _extending_vvp:

=============================
Contributing Cases and Suites
=============================

Contributions to MC/DC-VVP include benchmark models, reference data, processing
tools, and new suites. Together, these connect scientific questions to repeatable
transport calculations and processed evidence.
Use the :doc:`campaign operating guide <campaigns/index>` before
extending the workflow, and the :doc:`VVP overview <index>` for
the distinction between verification, validation, and performance.
The same extension principles apply to all three kinds of assessment.

Where content belongs
---------------------

MC/DC-VVP READMEs provide short operational entry points for setup, launch,
processing, output locations, and suite-specific requirements.
Sphinx provides the complete guide: those operations, the scientific rationale,
organization, mathematical definitions, interpretation, and extension procedures.
Keep operational summaries consistent with the scripts and link to the relevant
Sphinx page for the full explanation.

A suite README should make its local commands and input requirements easy to find.
A Sphinx case page should explain the problem, source and material assumptions,
quantity of interest, reference provenance, statistical treatment, limitations,
and published results.
Do not imply that a model produces a validated observable when it only generates
an intermediate spectrum or lacks a required response function.

Adding a verification case
--------------------------

Choose the suite by transport formulation and calculation type: analytical
multigroup fixed-source, analytical multigroup eigenvalue, or multigroup
code-to-code. Classify cases by the formulation being assessed, and document
the particle species and interaction data in each case definition.
Start from a neighboring case in the appropriate suite. Keep its file interfaces so the common runner and processor can find it:

.. code-block:: text

   cases/<case>/
     input.py        MC/DC model and simulation execution
     reference.py    Analytical solution or external comparison-data retrieval
     process.py      Convergence calculation across configured sampling levels
     plot.py         A selected-output comparison or animation

#. Define the model and scored observable in ``input.py``.
   Use the existing simulation command-line overrides for sampling and output
   naming rather than embedding a private launch workflow.
   Fixed-source runners pass ``--N_particle``; eigenvalue runners pass
   ``--N_active``. Keep eigenvalue particles per cycle and inactive cycles fixed
   after a preliminary source-convergence study.
#. Put reusable suite data in ``data/`` where appropriate.
   Resolve input paths relative to the source file, not an assumed working
   directory. Record data provenance and any required temperature, isotope,
   normalization, or unit conventions.
#. Implement the reference and processing interfaces of the chosen suite.
   Analytical processing compares against the reference solution; code-to-code
   processing aligns the participating grids and normalizations and constructs
   the fixed largest-sample comparison reference.
   Follow the :doc:`standardized error metrics <verification/error_metrics>`
   for bin normalization, absolute maxima, and zero-reference handling, and
   identify any existing processor conventions that require updating.
   For external files, record their archive and checksums. Use reference data
   for the actual observable, not a similarly named baseline problem.
#. Register the case in ``task.yaml``.
   Fixed-source/code-to-code entries use ``logN_min``, ``logN_max``, ``N_task``,
   and ``walltime``; eigenvalue entries use ``N_active_min``,
   ``N_active_max``, ``N_task``, and ``walltime``.
   The :doc:`verification workflows <campaigns/verification>`
   show both schemas and task-generation rules.
   Require positive walltime hours in each task and document its calibration in
   task-file comments; launch overrides do not replace this requirement.
#. Check a small task locally, inspect its tally grids and normalization, then
   process a small sequence before allocating a full campaign.
   Exercise partial restart and an already complete case; a repeated unchanged
   launch should skip existing outputs.
#. Add a case page under the matching Sphinx verification directory and register
   it in the corresponding suite landing page and its toctree.
   Add its short entry to the suite README and link to the full explanation.

For example, adding a fixed-source ``my_case`` requires a matching
``cases/my_case/`` directory and a task entry such as:

.. code-block:: yaml

   my_case:
     logN_min: 3
     logN_max: 5
     N_task: 3
     walltime: 1.0  # Illustrative hours; calibrate for the case job.

Those sampling levels are illustrative; choose production levels based on the
observable's variance and required statistical resolution.

Adding a validation case
------------------------

The SINBAD suite uses ``<benchmark>/<case>/input.py`` and ``plot.py`` rather than
the verification suite's four-file layout.
Register cases in the nested benchmark/case mapping shown in
:doc:`validation/sinbad`.
Its common runner supplies ``N_particle`` per batch; the input defines batch count.
The shared package helpers in ``validation/sinbad/util.py`` locate the locally
available official package during processing.

Document the experimental observable, its normalization and defensible comparison
range, source law, material representation, geometry approximations, and reported
uncertainty interpretation.
Separate statistical error from source, detector, model, and nuclear-data
uncertainties. State missing covariance information and avoid assigning a
one-standard-deviation meaning to an unspecified combined uncertainty.

Do not redistribute licensed experimental packages.
Keep package lookup and experimental-table parsing in the processor/plotter so
transport inputs can run with their required nuclear-data library independently
of the experimental-table collection.
Models that lack the capability to compute the benchmark observable should be
identified as such and left outside the executable validation task configuration.
The :doc:`existing benchmark pages <validation/index>` illustrate
this distinction.

Adding a performance case or method
-----------------------------------

Performance suites need an ``input.py`` plus an entry in their ``task.yaml``.
They use common runners and processors instead of a per-case convergence script.
The input must save full tally results as well as standard performance metadata.
The :doc:`serial workflow <performance/serial>` uses separate
Numba/Python task blocks; the
:doc:`parallel workflow <performance/parallel>` uses a baseline
and a node-count/workload matrix.

For a method comparison, preserve the physical observable, tally scores and grids,
normalization, and batch count. In parallel tasks, assign the same ``comparison``
label and distinct ``method`` labels to the inputs being compared.
Calibrate platform-specific baselines explicitly when appropriate.
A larger tracking rate alone does not establish a larger FOM; retain the
statistical-efficiency measurements and common fixed-reference calculation.

Parallel inputs execute in isolated working directories.
Use source-relative or absolute paths for data, preserve the absolute
platform-specific output path, and keep each job's workspace accessible to all
its ranks and allocated nodes.
Test restart validation and ensure incompatible old outputs are rejected.

Adding a suite
--------------

Choose a directory reflecting the assessment and physics, following the existing
``verification/``, ``validation/``, or ``performance/`` organization.
A suite is independently runnable and provides at least:

.. code-block:: text

   README.md        Operational entry point and links to full documentation
   launch.py        Study generation and platform-specific launch
   process.py       Processing entry point, writing suite results/
   cleanup.py       Removal of generated suite products

Provide ``task.yaml``, model directories, a runner, and helpers as appropriate.
Use neighboring suite implementations as interface examples rather than assuming
all suites have the same task schema or job granularity.

The top-level integration contract is:

- Register the relative suite path in ``configs/launch_config.py.template`` and
  in the operator's local ``configs/launch_config.py`` when testing it.
- Provide ``enabled`` and a nonempty ``platforms`` list.
  Accept ``--platform`` and only the additional options configured for that suite:
  ``--N_node``, ``--walltime``, or ``--N_node_max``.
- Use shared platform and user configurations for resource and interpreter
  settings. Preserve the active environment unless an interpreter override is
  configured.
- Save effective launch and task configurations with the Maestro run.
  Have processing use those saved configurations rather than a later edit of
  the suite's task file.
- Write processed products under suite ``results/`` so the root processor can
  collect them under the same relative suite path.
  The root processor visits all registered suites, independently of launch
  enablement, and replaces previously collected results for a processed suite.
- Make cleanup safe to scope to generated products and keep reference data and
  model definitions. The root cleaner visits all registered suites.
- Define task-completion checks and collision behavior so resuming work does not
  silently mix distinct studies. Preserve platform identity when combining
  hardware measurements.

If a new suite needs additional top-level launch or processing arguments, update
those dispatchers deliberately. Adding a configuration key does not automatically
make the root launcher forward it.
Add operating guidance and a linked scientific section to Sphinx, then test the
suite both independently and through the root workflow with a small workload.

Publishing new figures
----------------------

After reviewing the numerical evidence, collect the suite results and use
:doc:`campaigns/processing` to prepare flat release assets.
For a new documentation figure, use the actual generated filename, a descriptive
caption, and alt text identifying the quantity and comparison.
Do not embed an assumed asset before it exists or use a notional illustration
as measured evidence.

Keep development evidence associated with ``vvp-dev`` and released evidence with
the corresponding version tag. Archive the MC/DC and MC/DC-VVP revisions and the
saved launch/task metadata needed to interpret the campaign.
Follow :doc:`../../developer_guide/release_policy` when publishing evidence with a software release.
