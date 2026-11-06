======================
Launching and Resuming
======================

Repository-wide launch
----------------------

From the MC/DC-VVP repository root, launch enabled suites assigned to local
execution:

.. code-block:: console

   python launch.py

On a configured HPC machine, select its platform name:

.. code-block:: console

   python launch.py --platform tuolumne

Run separately on each target machine.
The top-level launcher forwards suite-specific node and walltime settings from
``configs/launch_config.py`` and invokes each suite with the active interpreter.
Submission is not completion: wait for the Maestro studies and scheduled jobs
to finish before processing or cleaning their outputs.

Individual suites
-----------------

An integrated suite can also run independently of the top-level launch selection.
For example, from the MC/DC-VVP repository root:

.. code-block:: console

   cd verification/analytical/multigroup/fixed_source
   python launch.py
   python launch.py --platform dane --N_node 1 --walltime 1.0

These are alternative local and HPC launches, not consecutive steps required
for one campaign.
The suite reads its ``task.yaml`` and shared platform/user settings.
Verification and SINBAD launchers support ``--N_node`` and ``--walltime``;
performance launchers have different resource interfaces described in their
:doc:`serial <../performance/serial>` and :doc:`parallel <../performance/parallel>`
workflows.
Only parallel performance provides the documented ``--dry_run`` study preview;
that option is not a general top-level launcher flag.

Resuming work
-------------

Verification runners skip a task when its expected ``output_<sampling-level>.h5``
exists. A case with all expected outputs is omitted from the generated study,
and a fully complete suite does not launch another study.
Partially complete cases execute only their missing sampling levels.
SINBAD similarly skips an existing ``output_<N_particle>.h5``.

These existence checks do not establish that an old output matches a changed
physical model or task configuration.
For a changed campaign, preserve the previous outputs and saved configurations
in a separate checkout or campaign copy before generating new results.

Performance launchers additionally validate output metadata and full tally
content. Outputs incompatible with the requested baseline, execution resources,
or required tallies are rejected rather than silently overwritten.
Performance paths distinguish platforms, but do not distinguish arbitrary
revisions of a campaign on the same platform.

Cleaning a campaign
-------------------

After all jobs have stopped, run from the repository root:

.. code-block:: console

   python cleanup.py

The root cleaner visits **every suite registered in** ``LAUNCH_CONFIG``, including
entries whose ``enabled`` value is false, then removes the collected root
``results/`` directory.
A suite's ``python cleanup.py`` limits cleanup to that suite.
Cleaners remove generated outputs, processed figures, Maestro run directories,
and generated study files; parallel performance also removes retained workspaces
and caches. Reference data and model inputs are retained.

Use cleanup when a fresh run is intended. Relaunching without cleanup is the
normal way to fill missing tasks in an otherwise unchanged campaign.
