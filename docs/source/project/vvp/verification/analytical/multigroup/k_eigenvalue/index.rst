.. _vvp_multigroup_eigenvalue:

==================================
Analytical Multigroup k-Eigenvalue
==================================

This suite compares multiplication factors and fundamental-mode flux shapes with
analytical and semi-analytical solutions of multigroup transport eigenproblems.
Finite slabs test leakage-dependent modes; infinite-medium cases isolate the
energy-dependent eigenproblem.

The suite assesses the shared, particle-agnostic multigroup transport formulation.
Each case specifies its particle species, interaction data, and scored observables.

Evidence and interpretation
---------------------------

Tasks increase active cycles while keeping particles per cycle and inactive
cycles fixed within each case. A preliminary source-convergence study establishes
those fixed settings before interpreting active-cycle convergence.
Error plots test :math:`N_\mathrm{active}^{-1/2}` behavior; multiplication-factor
uncertainty plots compare estimates and reported uncertainties with the reference.
Flux-shape and cycle-history plots show the highest-statistics result.

See :ref:`reading_verification_results` and :ref:`verification_error_metrics`
for convergence interpretation, normalization, and existing metric conventions.

Finite-slab criticality problems
--------------------------------

Homogeneous and heterogeneous one-group slabs test leakage-dependent multiplication factors and fundamental modes.

- :doc:`Homogeneous one-group slab criticality <one_group_slab>`
- :doc:`Kornreich-Parsons slab <kornreich>`

Infinite-medium SHEM-361 variants
---------------------------------

These cases test dominant 361-group eigenpairs below and above criticality without spatial leakage.

- :doc:`Subcritical <inf_shem361_subcritical>`
- :doc:`Supercritical <inf_shem361_supercritical>`

Operating and contributing
--------------------------

See :ref:`vvp_run_eigenvalue` for task configuration, launch commands,
reference preparation, and processing.
Use :doc:`../../../../contributing` to add a case or extend a suite.

.. toctree::
   :hidden:

   one_group_slab
   kornreich
   inf_shem361_subcritical
   inf_shem361_supercritical
