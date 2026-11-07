====================
Running SINBAD Cases
====================

The suite lives in ``validation/sinbad`` in MC/DC-VVP.
Its :doc:`benchmark pages <index>` describe the experimental configurations,
source approximations, quantities of interest, and uncertainties.
This page covers the executable campaign and its processing requirements.

Inputs and experimental data
----------------------------

The repository supplies MC/DC inputs and processing scripts but does not
redistribute the licensed SINBAD packages.
Obtain the official packages locally before processing comparisons that read
their experimental tables.
Launching transport requires MC/DC's native nuclear-data library; the experimental
packages are required for the comparison processing, rather than the launch.

Runnable cases contain ``input.py`` and ``plot.py`` inside their benchmark and
case directories. For example:

.. code-block:: text

   validation/sinbad/
     oktavian-si-60/neutron_leakage/
     fng-tud-sic/p1/
     fng-tud-sic/p2/
     fng-tud-sic/p3/
     fng-tud-sic/p4/
     task.yaml
     launch.py
     run_case.py
     process.py
     cleanup.py
     util.py

Task configuration
------------------

Unlike verification suites, ``task.yaml`` is nested by benchmark directory and
case directory. Each entry specifies particles per batch and required walltime hours:

.. code-block:: yaml

   oktavian-si-60:
     neutron_leakage:
       N_particle: 1_000_000_000
       walltime: 8.16
   fng-tud-sic:
     p1:
       N_particle: 1_000_000_000
       walltime: 0.24
     p2:
       N_particle: 1_000_000_000
       walltime: 0.21
     p3:
       N_particle: 1_000_000_000
       walltime: 0.21
     p4:
       N_particle: 1_000_000_000
       walltime: 0.21

The batch count remains in each input.
Enable only cases whose benchmark observable can be calculated directly.
A case is skipped if its expected ``output_<N_particle>.h5`` already exists.

Launching
---------

From the suite directory, launch locally:

.. code-block:: console

   python launch.py

Or launch on a configured cluster:

.. code-block:: console

   python launch.py --platform dane --N_node 10

Each case uses the specified node count and all configured CPU cores per node.
Each task requires positive ``walltime`` hours. The example requests 8.16 hours
for OKTAVIAN, 0.24 hours for P1, and 0.21 hours for P2--P4.
Document measured runtime, resources, workload, and margin in task comments.
A launch ``--walltime`` overrides every case job uniformly; requests are rounded
up to scheduler resolution and capped at the platform maximum.
Local execution uses one node and ignores time limits, but still requires valid
task values even when an override is supplied.

The saved ``launch_config.yaml`` records platform, nodes, processes, Python
interpreter, optional walltime override, and effective case walltimes.
For repository-wide operation, enable ``validation/sinbad`` in
``configs/launch_config.py`` and set ``sinbad_root`` to the local SINBAD root.
See :doc:`../campaigns/setup` for platform and user settings.

Processing
----------

After completion, run from the suite directory:

.. code-block:: console

   python process.py /path/to/sinbad

To select a particular study instead of the latest:

.. code-block:: console

   python process.py /path/to/sinbad maestro_run_<timestamp>

The processor resolves each official package below the supplied root, runs its
case plotter, and collects figures under ``results/<benchmark>/<case>/``.
It preserves effective launch and task configurations in ``results/``.
The repository-wide processor collects this tree under
``results/validation/sinbad/``.

Run ``python cleanup.py`` from this suite only after jobs have stopped.
It removes outputs and generated figures for configured cases, processed results,
Maestro directories, and ``study.yaml``.
Incomplete model inputs outside ``task.yaml`` are left untouched.
