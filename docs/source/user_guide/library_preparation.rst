.. _library_preparation:

==========================
Native Library Preparation
==========================

MC/DC reads native transport data from HDF5 files in the directory named by ``MCDC_LIB``.
Separate neutron, electron, and proton generators convert source libraries into this common output directory.
Prepare the data for the particle types, materials, and temperatures used in your model.
Multigroup models supplied directly through the Python interface do not require these native libraries.

Library Sources
---------------

The current generators expect the following source data:

.. list-table::
   :header-rows: 1
   :widths: 15 40 45

   * - Particle
     - Source
     - Required input
   * - Neutron
     - `ENDF/B-VIII.1 Lib81 ACE downloads <https://www.nndc.bnl.gov/endf-library/B-VIII.1/ACE/>`_
     - Individual continuous-energy neutron ACE files, organized by nuclide and temperature.
   * - Electron
     - EPRDATA14, included with the `MCNP6.2 distribution <https://rsicc.ornl.gov/codes/ccc/ccc8/ccc-850.html>`_
     - The concatenated EPRDATA14 ACE file containing elemental electron, photon, and relaxation data.
   * - Proton reactions
     - `TENDL-2021 downloads <https://tendl.imperial.ac.uk/tendl_2021/tar.html>`_
     - The proton ACE archive, ``TENDL-ACE-p.tgz``; extract its individual ``.ace`` files.
   * - Proton stopping power
     - `NIST PSTAR <https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html>`_
     - Plain-text tables containing energy and total stopping power for each required element.

The `LANL EPRDATA format report <https://mcnp.lanl.gov/pdf_files/TechReport_2024_LANL_LA-UR-24-30590_LivelyHaeck.pdf>`_ describes the electron library format.
The electron generator needs the ACE data file itself, not a preconverted HDF5 library from another transport code.
For neutron inputs, select the continuous-energy neutron library rather than the separate thermal-scattering archive.
The current neutron filename decoder recognizes the Lib81 temperature suffixes.

Prerequisites and Output Directory
----------------------------------

Use a source checkout of MC/DC containing ``tools/data_library_generator/``.
Install `ACEtk <https://github.com/njoy/ACEtk>`_ with its Python bindings and the generator dependencies in your Python environment:

.. code-block:: sh

   python -m pip install h5py numpy tqdm
   python -c "import ACEtk"
   export MCDC_LIB=/path/to/mcdc/library

Run the commands below from the MC/DC repository root.
All three generators use ``MCDC_LIB`` for output, and simulations use the same variable to locate their data.
The generators create the directory if necessary.

Neutron Data
------------

Point ``MCDC_ACELIB`` to a directory containing the individual neutron ACE files:

.. code-block:: sh

   export MCDC_ACELIB=/path/to/ace/Lib81
   python tools/data_library_generator/neutron/generate.py

The generator converts nuclear cross sections, outgoing distributions, and fission production data into nuclide-temperature files such as ``U235-293.6K.h5``.
Keep input directories limited to ACE tables supported by this generator.

Electron Data
-------------

Point ``MCDC_ACELIB_ELECTRON`` to the concatenated data file, not its containing directory:

.. code-block:: sh

   export MCDC_ACELIB_ELECTRON=/path/to/eprdata14
   python tools/data_library_generator/electron/generate.py

The generator writes elemental files such as ``Al.h5`` with electron reaction and atomic relaxation data.
These coexist with neutron and proton nuclide files in the same output directory.

Proton Data
-----------

Provide a directory of proton ACE files and a directory of PSTAR tables:

.. code-block:: sh

   export MCDC_ACELIB_PROTON=/path/to/tendl2021/proton/ace
   export MCDC_PSTAR_LIB=/path/to/pstar
   python tools/data_library_generator/proton/generate.py

Save each elemental PSTAR table as ``<Symbol>.txt``, for example ``H.txt`` or ``Al.txt``.
Each numeric row must contain exactly two columns: energy in MeV and total stopping power in MeV cm²/g.
Select only total stopping power when exporting the table; do not include additional numeric columns for ranges or separate stopping-power components.

ACE files supply nuclear reactions and secondary-product distributions; PSTAR supplies stopping power for condensed interactions.
The generator also supports a predefined set of stopping-power-only isotopes when no ACE reaction data is available.
Radiation lengths come from the generator's embedded table and require no additional source file.

Updating and Combining Contributions
------------------------------------

All generators accept ``--verbose`` for detailed output and ``--rewrite`` to replace existing generated data:

.. code-block:: sh

   python tools/data_library_generator/neutron/generate.py --rewrite --verbose

Neutron and proton generators skip existing contributions by reaction group, rather than skipping every existing nuclide file.
They preserve the other particle's data when adding or replacing their own contribution.
PSTAR-only updates preserve reaction groups, and absent new stopping-power input does not erase an existing table.
The electron generator skips existing elemental files unless ``--rewrite`` is supplied, in which case it regenerates those files.
Run generators sequentially when they target the same output directory.

Neutron and proton contributions combine only when their nuclide and temperature filenames match.
The proton generator currently writes ``0.0K`` files; it does not automatically attach them to neutron files at other temperatures.
Coupled neutron-proton transport requires both reaction contributions in the nuclide-temperature file selected by the simulation.
Do not rename files across temperatures to bypass this requirement.

Shared nuclide identity must agree when contributions are combined.
Floating-point metadata is checked with a relative tolerance of ``1e-5``; conflicts leave the existing file unchanged.
Per-particle provenance is recorded under ``provenance/neutron`` and ``provenance/proton``, with PSTAR source information attached to its stopping-power data.
The root ``fissionable`` flag describes neutron-induced fission.

Using the Prepared Library
--------------------------

Keep ``MCDC_LIB`` set when running MC/DC:

.. code-block:: sh

   export MCDC_LIB=/path/to/mcdc/library
   python input.py

See :doc:`materials_and_multigroup` for specifying native material compositions and temperatures.
The generator READMEs under ``tools/data_library_generator/`` describe the output schemas and particle-specific details.
