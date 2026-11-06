.. _vvp_campaigns:

===================
Operating Campaigns
===================

A campaign selects suites, generates their tasks, executes them through Maestro,
and processes the resulting transport outputs into figures and tables.
The instructions here cover both repository-wide operation and individual suites.
Commands run in an MC/DC-VVP checkout, rather than the MC/DC source repository.

Start with :doc:`setup`, then follow :doc:`running` and :doc:`processing`.
Use :doc:`verification` for analytical and code-to-code task configuration,
:doc:`../validation/sinbad` for experimental comparisons, and the
:doc:`serial <../performance/serial>` and :doc:`parallel <../performance/parallel>`
performance workflows for timing and scaling studies.
To add a case or suite, see :doc:`../contributing`.

Small local example
-------------------

With the MC/DC environment and dependencies from :doc:`setup` available, run a
small analytical fixed-source demonstration from the MC/DC-VVP repository root.
Use this minimal ``configs/launch_config.py``:

.. code-block:: python

   LAUNCH_CONFIG = {
       "verification/analytical/multigroup/fixed_source": {
           "enabled": True,
           "platforms": ["local"],
           "N_node": 1,
       },
   }

Set ``verification/analytical/multigroup/fixed_source/task.yaml`` to:

.. code-block:: yaml

   slab_absorbium:
     logN_min: 3
     logN_max: 4
     N_task: 3

This selects one case and three particle-count tasks between 1,000 and 10,000.
It demonstrates the workflow; use a suitably resolved sampling campaign when
assessing convergence.
Launch the study:

.. code-block:: console

   python launch.py

After the Maestro study finishes, process and collect the results:

.. code-block:: console

   python process.py

Plots and saved configurations are collected under
``results/verification/analytical/multigroup/fixed_source/``.
Run ``python prepare_release.py`` to copy its PNG/GIF figures to flat assets in
``release/``. See :doc:`processing` for output replacement, archival requirements,
and publication, and :doc:`running` for restarting or cleaning this campaign.

Repository layout
-----------------

.. code-block:: text

   configs/               Platform, user, and launch configurations
   verification/          Analytical and code-to-code suites
   validation/            Experimental benchmark suites
   performance/           Serial and parallel performance suites
   results/               Collected processed results, organized by suite
   release/               Flat PNG and GIF assets prepared for publication
   launch.py              Launch enabled suites for the selected platform
   process.py             Process registered suites and collect results
   cleanup.py             Clean registered suites and collected results
   prepare_release.py     Prepare release assets from collected figures

Within a suite, ``task.yaml`` defines the cases and sampling levels.
``launch.py`` generates a Maestro study, ``run_case.py`` executes its transport
work, and ``process.py`` produces the comparisons or performance measurements.
A saved ``maestro_run_*`` directory records the launch and task configurations.
Keep these records together with the simulation outputs so that later processing
can recover the configuration actually used.

The suite/case/task hierarchy describes scientific organization, not a fixed
one-to-one mapping to scheduler jobs.
Verification normally uses one Maestro step per case with sequential sampling
levels; serial performance uses one step per case and execution mode; parallel
performance uses one step per node-count/workload point.

.. toctree::
   :maxdepth: 1

   setup
   running
   processing
   verification
