.. _vvp_multigroup_fixed_source:

==================================
Analytical Multigroup Fixed-Source
==================================

This suite compares steady-state and transient fixed-source calculations with
analytical, manufactured, and semi-analytical reference solutions.
It exercises source sampling, geometry tracking, material lookup, particle
interactions, and tallying across one-group and multigroup problems.

The suite assesses the shared, particle-agnostic multigroup transport formulation.
Each case specifies its particle species, interaction data, and scored observables.

Evidence and interpretation
---------------------------

Tasks increase the number of source particles while preserving the problem and
tally definition. Convergence plots compare error trends with the expected
:math:`N^{-1/2}` behavior when statistical uncertainty dominates.
Largest-sampling-level comparisons show the solution's spatial, angular, or
energy dependence; transient comparisons use animations.

See :ref:`reading_verification_results` and :ref:`verification_error_metrics`
for interpretation and the distinction between standardized metrics and existing
figure conventions.

Slab problems
-------------

One-dimensional cases exercise geometry, material interfaces, source distributions, energy groups, and time dependence.

- :doc:`Absorbium slab <slab_absorbium>`
- :doc:`Manufactured two-group slab <mms_two_group_slab>`
- :doc:`Time-dependent isotropic-beam slab <slab_isobeam_td>`
- :doc:`Reed's slab problem <reed>`

AZURV1 variants
---------------

Variants of one transient problem test multiplication, time censuses, tally recombination, and variance reduction.

- :doc:`Base critical problem <azurv1>`
- :doc:`Subcritical medium <azurv1_sub>`
- :doc:`Supercritical medium <azurv1_super>`
- :doc:`Time censuses <azurv1-census>`
- :doc:`Census-based tallying <azurv1-census-tally>`
- :doc:`Basic variance reduction <azurv1-basic_techniques>`
- :doc:`Weight windows <azurv1-weight_windows>`

Infinite-medium SHEM-361 variants
---------------------------------

The 361-group cases test steady-state and transient transport, energy-dependent weight windows, and censuses without spatial leakage.

- :doc:`Steady-state <inf_shem361>`
- :doc:`Steady-state with weight windows <inf_shem361-weight_windows>`
- :doc:`Transient <inf_shem361_td>`
- :doc:`Transient with time censuses <inf_shem361_td-census>`

Operating and contributing
--------------------------

See :ref:`vvp_run_fixed_source` for task configuration, launch commands,
reference preparation, and processing.
Use :doc:`../../../../contributing` to add a case or extend a suite.

.. toctree::
   :hidden:

   slab_absorbium
   mms_two_group_slab
   slab_isobeam_td
   reed
   azurv1
   azurv1_sub
   azurv1_super
   azurv1-census
   azurv1-census-tally
   azurv1-basic_techniques
   azurv1-weight_windows
   inf_shem361
   inf_shem361-weight_windows
   inf_shem361_td
   inf_shem361_td-census
