=====================
Setting Up a Campaign
=====================

Environment and inputs
----------------------

Obtain a checkout of `MC/DC-VVP <https://github.com/mcdc-project/mcdc-vvp>`_
and activate the Python environment containing the MC/DC version to assess.
Install MC/DC using the :doc:`installation guide <../../../user_guide/getting_started/installation>`.
Maestro must be available to the launcher; the workflow scripts also use NumPy
and PyYAML, while processing uses h5py, Matplotlib, SciPy, and Pillow for animated
figures.
For an environment that already contains MC/DC, the additional tools can be
installed with:

.. code-block:: console

   python -m pip install maestrowf numpy pyyaml h5py matplotlib scipy pillow

Suite inputs may need transport libraries or external comparison data:

- Analytical SHEM-361 cases use the checked-in multigroup dataset in
  ``verification/analytical/multigroup/fixed_source/data/SHEM-361.npz``.
  The eigenvalue suite uses that same dataset.
- C5G7 verification uses its suite's ``data/MGXS-C5G7.h5``.
  Code-to-code processing also needs the participating-code reference files;
  see :doc:`verification`.
- SINBAD calculations require the native nuclear-data library accessible to
  MC/DC. Experimental-table processing requires a licensed local SINBAD
  collection; see :doc:`../validation/sinbad`.
- The Kobayashi detector performance inputs define their multigroup materials
  directly and require no external data directory.

Use :doc:`../../../user_guide/materials_and_data` for transport-library setup
and material-data conventions.

Campaign selection
------------------

From the MC/DC-VVP repository root, create the local configuration:

.. code-block:: console

   cp configs/launch_config.py.template configs/launch_config.py

Edit ``LAUNCH_CONFIG`` to enable the desired suites and choose their platforms.
For example, a local fixed-source campaign can use:

.. code-block:: python

   LAUNCH_CONFIG = {
       "verification/analytical/multigroup/fixed_source": {
           "enabled": True,
           "platforms": ["local"],
           "N_node": 1,
           "walltime": None,
       },
   }

The keys are suite paths relative to the repository root.
Each entry must have a nonempty ``platforms`` list, containing ``"local"`` or
names defined in ``configs/platform_config.py``.
Only enabled suites that list the selected platform are launched.
The list selects compatible hosts; it does not submit jobs to another machine.

``N_node`` sets the nodes assigned to a verification or SINBAD case.
Serial performance uses one process on one exclusive node.
Parallel performance instead takes ``N_node_max`` and generates power-of-two
node counts through that limit; see :doc:`../performance/parallel`.
Set ``sinbad_root`` on the SINBAD entry before repository-wide processing.
A relative SINBAD root is resolved against the MC/DC-VVP repository root.

Platform and user settings
--------------------------

``configs/platform_config.py`` defines scheduler commands, process launchers,
CPU/GPU counts, maximum nodes, maximum walltime, and scheduler time resolution.
The supplied configurations include Dane (Slurm), Lassen (LSF), and Tuolumne
(Flux). Use the checked-in configuration for the checkout being executed.

For HPC execution, create the user configuration:

.. code-block:: console

   cp configs/user_config.py.template configs/user_config.py

Replace the template entry with the target platform and your account settings:

.. code-block:: python

   USER_CONFIG = {
       "dane": {
           "account": "your_account",
           "queue": "your_queue",
           "reservation": None,
           "maestro_python": None,
           "mcdc_python": None,
       },
   }

``maestro_python`` and ``mcdc_python`` can select separate interpreters for
workflow orchestration and transport execution.
With no override, the launcher uses its active environment.
The transport interpreter and required data must also be accessible on compute
nodes. Load the required GPU runtime and scheduler environment before launching;
on Tuolumne this includes ROCm and, when required by the Python/Flux installation,
running the launcher through ``flux env``.

Walltime and workload
---------------------

Every task definition must include a finite, positive ``walltime`` in hours.
For verification, this covers the entire case job and its sequential sampling
levels; validation uses one value per case, and serial performance uses one per
execution-mode block. Parallel performance requires a mapping from every workload
multiplier to its requested hours; see :doc:`../performance/parallel`.

A non-``None`` ``walltime`` in ``configs/launch_config.py``, or a suite-level
``--walltime HOURS``, overrides all selected job durations uniformly.
Without an override, each job uses its task value. There is no platform-maximum
fallback: a missing or invalid task value is an error even when an override is
provided. Validation also applies to local launches and already completed tasks.
Local execution does not enforce a time limit.
Scheduler requests are rounded up to the scheduler's time resolution and capped
at the platform maximum.

Document calibration in ``task.yaml`` comments: platform and resources, workload,
measured runtime, and the chosen margin.
Use ``None`` for the launch-level walltime to select the task-defined requests.

Review the suite's ``task.yaml`` before launching: repository configurations can
represent large campaigns. Choose sampling levels and resource counts appropriate
to the intended study, and keep physical inputs fixed when comparing sampling,
execution modes, or hardware.
