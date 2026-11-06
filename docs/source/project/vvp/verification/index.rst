.. _project_vvp_verification:

============
Verification
============

Verification assesses whether MC/DC correctly solves its stated mathematical and numerical problems.
It uses well-defined problems with a definite target solution.
In analytical verification, that solution is known explicitly or supplied by a
manufactured or accurately converged semi-analytical reference.
In code-to-code verification, the exact solution need not be known explicitly,
but equivalent formulations of the same well-posed problem have a common solution
to which all participating codes should converge as statistical and numerical
errors are reduced.

Establishing that a transport code solves its intended problems correctly
requires verification of both its individual components and their behavior
together in a complete simulation.
:doc:`Unit tests <../../../developer_guide/contributing/unit_testing>` isolate
components, interfaces, and edge cases, making it easier to identify defects in
specific parts of the implementation.
MC/DC-VVP complements these tests with full-transport problems that exercise
model preparation, source sampling, geometry tracking, particle interactions,
and tallying together.
Analytical and code-to-code comparisons then assess whether the resulting
simulation estimates approach the intended solution.

Verification suites and cases
-----------------------------

The verification campaign is organized by verification method, transport
formulation, suite, and case.
The multigroup suites assess the shared, particle-agnostic transport formulation:
interaction cross sections and associated distributions specify the problem,
with cross sections set to zero for interaction channels that are not relevant.
Cases are classified by the transport formulation being verified.
Each case specifies its particle species and interaction data as part of the
problem definition.
See :ref:`reading_verification_results` for guidance on interpreting convergence
curves, uncertainty plots, and solution comparisons.
Use :doc:`../campaigns/verification` for task configuration, launch commands,
reference-data preparation, and processing.

**Analytical verification** uses known reference solutions:

- :doc:`Multigroup fixed-source <analytical/multigroup/fixed_source/index>` tests
  steady-state and transient solutions as source-particle counts increase.
- :doc:`Multigroup k-eigenvalue <analytical/multigroup/k_eigenvalue/index>` tests
  multiplication factors and fundamental modes as active-cycle counts increase.

**Code-to-code verification** compares independent implementations, checking
convergence toward a common solution:

- :doc:`Multigroup fixed-source <code_to_code/multigroup/fixed_source/index>` compares
  steady-state and transient fixed-source solutions across participating codes.
- **Multigroup k-eigenvalue** (planned) will compare multiplication factors and
  fundamental modes across participating codes.

Each suite page describes its scope, evidence, case catalog, and operating workflow.

.. _reading_verification_results:

Reading verification results
----------------------------

Monte Carlo estimates fluctuate statistically, so a single close comparison is not sufficient evidence of convergence.
The verification suites repeat each case over a sequence of sampling levels and compare the observed error or difference with the expected inverse-square-root behavior.
For a source-particle count :math:`N`, ordinary Monte Carlo sampling predicts errors proportional to :math:`N^{-1/2}` when statistical uncertainty dominates.
The fixed-source suites therefore vary the number of source particles and
summarize discrepancies using a relative 2-norm and a maximum relative error or
difference. :doc:`error_metrics` defines the standardized bin-normalized 2-norm
and absolute maximum, including reference selection, bin weighting, and
zero-reference handling. It also records the different conventions used by
existing processors and published figures.

The :math:`k`-eigenvalue suite instead varies the number of active cycles while holding the particles per cycle and inactive cycles fixed within each case.
Its convergence plots test the corresponding :math:`N_\mathrm{active}^{-1/2}` behavior, while its uncertainty plots show whether the reference multiplication factor is consistent with the reported MC/DC uncertainty.

Solution plots use the task with the largest sampling effort to show where MC/DC agrees with or departs from the reference.
Time-dependent solution comparisons animate spatial or energy-dependent behavior throughout the transient.

For code-to-code cases, the arithmetic mean of the participating-code estimates at the largest common sampling level defines one fixed comparison reference for the entire convergence study.
The pairwise relative differences at every sampling level are normalized by that reference.
Decay proportional to :math:`N^{-1/2}` indicates that the codes are approaching the same solution at the expected statistical rate, whereas a persistent plateau may indicate a modeling discrepancy, implementation bias, or insufficiently resolved reference.

.. toctree::
   :hidden:

   analytical/multigroup/fixed_source/index
   analytical/multigroup/k_eigenvalue/index
   code_to_code/multigroup/fixed_source/index
   error_metrics
