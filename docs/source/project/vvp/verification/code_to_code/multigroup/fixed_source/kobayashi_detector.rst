.. _project_vvp_kobayashi_detector:

===================================
Kobayashi Dog-Leg Detector Response
===================================

**Case files:** `Code-to-code model <https://github.com/mcdc-project/mcdc-vvp/tree/dev/verification/code_to_code/multigroup/fixed_source/cases/kobayashi-detector>`_,
`serial performance input <https://github.com/mcdc-project/mcdc-vvp/tree/dev/performance/serial/cases/kobayashi-detector>`_,
and `parallel performance input <https://github.com/mcdc-project/mcdc-vvp/tree/dev/performance/parallel/cases/kobayashi-detector>`_.

Problem and observable
----------------------

This variant extends the :doc:`non-fissioning Kobayashi transient <kobayashi>`
with localized prompt-fission multiplication and a detector at the channel exit.
It retains the outer dimensions, dog-leg channel, source pulse, unit particle
speed, and reflective/vacuum boundary conditions of the baseline.
The material, geometry, and source definitions match the
:ref:`pulsed Kobayashi tutorial <example_kobayashi_td>`.

A 10 cm fuel cube occupies :math:`x\in[30,40]`, :math:`y\in[50,60]`,
:math:`z\in[0,10]` cm at the second channel turn.
Its scattering and fission cross sections are both 0.05 cm\ :sup:`-1`, with zero
capture and prompt-fission multiplicity 2.5.
A 10 cm detector cube at :math:`x\in[30,40]`, :math:`y\in[90,100]`,
:math:`z\in[30,40]` cm uses the shield material.

The tally scores capture integrated over the detector volume and each time bin,
normalized per source particle.
It replaces the baseline flux mesh and particle-density tallies with 500 unit-width
time bins over :math:`[0,500]`.
This longer, finer observation window resolves the outlet response after the
source pulse and subsequent multiplying transport.

Assessment uses
---------------

The code-to-code task configuration varies particles per batch and uses the
same fixed-reference pairwise convergence method as the other multigroup cases.
Detector-specific participating-code outputs are required; the case's reference
downloader identifies its baseline Kobayashi archive as a placeholder.
The published baseline flux comparisons on :doc:`kobayashi` do not verify the
fuel cube or detector observable.
See :doc:`../../../../campaigns/verification` before processing the detector case.

The same problem supplies the :doc:`serial <../../../../performance/serial>` and
:doc:`parallel <../../../../performance/parallel>` performance inputs.
Its compact time-history tally allows the performance suites to preserve full
tally output while measuring tracking rate, precision, and FOM across workloads.
