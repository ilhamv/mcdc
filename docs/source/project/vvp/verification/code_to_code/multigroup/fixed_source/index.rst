.. _vvp_multigroup_code_to_code:

====================================
Code-to-Code Multigroup Fixed-Source
====================================

This suite compares independently implemented codes on multigroup fixed-source problems
without tractable analytical solutions. The supplied comparisons use MC/DC and
OpenMC, with the formulation allowing additional participating codes.
The cases exercise moving geometry, delayed emission, void streaming, multiplying
transport, and space-time tallies.

The suite assesses the shared, particle-agnostic multigroup transport formulation.
Each case specifies its particle species, interaction data, and scored observables.

Evidence and interpretation
---------------------------

Tasks increase the participating calculations' sampling efforts together.
The arithmetic mean of their largest-common-sample estimates defines a fixed
comparison reference for every sampling level.
Convergence plots assess whether pairwise relative differences decrease at the
expected :math:`N^{-1/2}` rate, to the extent statistical uncertainty dominates.
Reference plots, solution comparisons, and difference plots show the transient
behavior and where discrepancies occur.

Agreement increases confidence but does not establish correctness because the
codes may share assumptions or biases. Comparisons can still expose differences
in modeling and implementation or identify consistent behavior.
See :ref:`reading_verification_results` and :ref:`verification_error_metrics`
for the common-reference definitions and metric conventions.

The detector variant requires detector-specific participating-code outputs;
its reference downloader currently identifies the baseline archive as a
placeholder. Its case page distinguishes the model from established baseline
comparison evidence.

Cases
-----

The cases compare multidimensional reactor and shielding transients.

- :doc:`C5G7 four-phase transient <c5g7-4phase>`
- :doc:`Kobayashi dog-leg transient <kobayashi>`
- :doc:`Kobayashi detector response <kobayashi_detector>`

Operating and contributing
--------------------------

See :ref:`vvp_run_code_to_code` for task configuration, launch commands,
reference preparation, and processing.
Use :doc:`../../../../contributing` to add a case or extend a suite.

.. toctree::
   :hidden:

   c5g7-4phase
   kobayashi
   kobayashi_detector
