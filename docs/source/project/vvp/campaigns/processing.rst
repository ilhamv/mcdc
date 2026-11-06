=============================
Processing and Release Assets
=============================

Process and collect
-------------------

After the studies finish, run from the MC/DC-VVP repository root:

.. code-block:: console

   python process.py

The processor visits every suite registered in ``configs/launch_config.py``,
regardless of its launch ``enabled`` flag.
When a saved Maestro run is available, it invokes the suite processor.
It then **moves** the suite's generated ``results/`` directory to
``results/<suite-path>/`` at the repository root, replacing any previous collected
results for that suite.
Existing suite results can be collected without a Maestro run; suites with neither
processable runs nor existing results are skipped.

Set ``validation/sinbad``'s ``sinbad_root`` before processing that suite.
Other external inputs, such as OpenMC statepoints, must also be available before
the relevant processor runs.

Selecting a study
-----------------

From a verification suite, process its latest study or an explicitly selected one:

.. code-block:: console

   python process.py
   python process.py maestro_run_<timestamp>

These are alternative study selections.
The processors use saved launch/task metadata, so editing the current
``task.yaml`` does not redefine the study being processed.
Keep the selected outputs and metadata together.
The resulting ``results/launch_config.yaml`` and ``results/task.yaml`` preserve
that configuration alongside the figures.

Performance processors accept multiple saved platform studies and default to
the latest study for each platform. SINBAD requires its local experimental-data
root before the optional study path. Follow the
:doc:`verification`, :doc:`SINBAD <../validation/sinbad>`,
:doc:`serial <../performance/serial>`, or :doc:`parallel <../performance/parallel>`
workflow for the applicable interface and output layout.

Preparing the assets
--------------------

After collecting the results, run from the repository root:

.. code-block:: console

   python prepare_release.py

The script recreates ``release/`` and copies every PNG and GIF under the root
``results/`` directory. It preserves the structured source results and flattens
each relative path by replacing directory boundaries with ``--``.
For example:

.. code-block:: text

   results/verification/analytical/multigroup/fixed_source/convergence/azurv1_flux.png
       becomes
   release/verification--analytical--multigroup--fixed_source--convergence--azurv1_flux.png

The script rejects a duplicate flattened name and reports the number and total
size of prepared assets. It does not upload them, change documentation URLs,
or package HDF5 outputs, CSV records, or configuration files.
Retain the numerical results and saved configurations separately when archiving
the campaign or regenerating its figures.

Publishing and inspecting
-------------------------

Publish development figures to the rolling
`VVP evidence for dev <https://github.com/mcdc-project/mcdc/releases/tag/vvp-dev>`_
release. Evidence for a released MC/DC version belongs to that version's release,
including patch releases, and is retained as a fixed record.
Record the MC/DC and MC/DC-VVP revisions used for the campaign in its publication
metadata.

A documentation figure URL combines the MC/DC release tag and flattened filename:

.. code-block:: text

   https://github.com/mcdc-project/mcdc/releases/download/<tag>/<asset-name>

Use ``vvp-dev`` for development evidence or the corresponding version tag for
released evidence. Retain the same filenames when publishing the corresponding
figures to a new release.
Check that the relevant documentation pages use the intended tag and that every
referenced asset exists and renders. Asset preparation alone does not switch the
links to a different release.
See the :doc:`release policy <../../../developer_guide/release_policy>` for the
software release workflow and :doc:`../contributing` for
integrating new case results into the documentation.
