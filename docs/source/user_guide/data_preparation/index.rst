.. _data-preparation:

================
Data Preparation
================

The generators in ``tools/data_library_generator/`` convert source libraries into the :doc:`MC/DC data format <../../reference/data_format/index>`.
Each projectile guide describes its source libraries, requirements, commands, and conversion rules.
The :doc:`format conventions <../../reference/data_format/conventions>` and :doc:`YAML illustrations <../../reference/data_format/illustrations>` describe the shared output structure.

.. list-table:: Source libraries
   :header-rows: 1
   :widths: 20 60 20

   * - Projectile
     - Source libraries
     - Source format
   * - :doc:`Neutron <neutron>`
     - Processed continuous-energy neutron libraries, such as `LANL Lib81 <https://nucleardata.lanl.gov/ace/lib81/>`_.
     - ACE
   * - :doc:`Photon <photon>`
     - `EPICS2025 <https://nuclear.llnl.gov/EPICS/index.html>`_: EPDL for photon interactions and optional EADL for shared atomic-relaxation data.
     - ENDF
   * - Electron
     - To be documented.
     - —
   * - :doc:`Proton <proton>`
     - `TENDL2021 <https://tendl.imperial.ac.uk/tendl_2021/tendl2021.html>`: TENDL2021 ACE data
     - `NIST PSTAR <https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html>`: stopping power tables

.. toctree::
   :maxdepth: 1

   neutron
   photon
   proton
