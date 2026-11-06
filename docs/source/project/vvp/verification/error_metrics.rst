.. _verification_error_metrics:

==========================
Standardized Error Metrics
==========================

Use two complementary summaries of each scored quantity: a bin-normalized
relative 2-norm (root-mean-square relative error) and a maximum absolute relative
error (infinity norm).
The first describes discrepancies across the solution; the second identifies the
worst included bin. Compute them separately for each observable, rather than
combining scores with different units or meanings.

The definitions below specify the common convention for processing campaigns.
The :ref:`published-metric-conventions` section distinguishes existing processor
formulas and figures that have not been regenerated with this convention.

Reference and comparison bins
-----------------------------

Let :math:`u_i(N)` be the calculated estimate in tally bin :math:`i` at sampling
level :math:`N`, and :math:`u_i^*` the fixed reference.
For analytical verification, the reference is the analytical, manufactured, or
accurately converged semi-analytical solution evaluated for the same observable
and bin averaging as the tally.
For code-to-code verification, use the arithmetic mean of all participating-code
estimates at the largest common sampling level.
Keep that reference fixed throughout the study.

Match units, source normalization, score definitions, tally grids, and any
case-defined eigenmode normalization before forming differences.
Use the fixed bin set

.. math::

   \mathcal I=\{i:u_i^*\ne 0\}, \qquad n=|\mathcal I|.

The default gives every included bin equal weight, :math:`1/n`, including on
nonuniform space, time, or energy grids. Thus it summarizes the scored bins,
rather than approximating a volume-, time-, or energy-weighted integral.
Any physically weighted variant must be identified separately with its weights
and normalization; it must not silently replace the default.

Analytical error
----------------

Define the pointwise relative error and its two summaries by

.. math::

   r_i(N)=\frac{u_i(N)-u_i^*}{u_i^*},
   \qquad
   E_2(N)=\left[\frac{1}{n}\sum_{i\in\mathcal I}r_i(N)^2\right]^{1/2},
   \qquad
   E_\infty(N)=\max_{i\in\mathcal I}|r_i(N)|.

The factor :math:`1/n` is inside the square root: this is
:math:`\|r\|_2/\sqrt n`, not :math:`\|r\|_2/n`.
It removes the direct square-root dependence on the number of included bins
that occurs in an unscaled discrete 2-norm.
It does not make different tally grids or observables physically equivalent.
Both metrics are dimensionless fractions; multiply by 100 only when displaying
percentages, and identify that choice on the axis.
For a scalar such as the multiplication factor, both reduce to
:math:`|k-k^*|/|k^*|`.

Pairwise code-to-code difference
--------------------------------

For :math:`M` participating codes and :math:`P=M(M-1)/2` distinct pairs, define

.. math::

   u_i^*=\frac{1}{M}\sum_{m=1}^{M}u_i^{(m)}(N_{\max}),
   \qquad
   d_i^{(a,b)}(N)=\frac{u_i^{(a)}(N)-u_i^{(b)}(N)}{u_i^*}.

Apply the same bin normalization to every pair, then average squared 2-norms
across pairs:

.. math::

   D_2(N)=\left[\frac{1}{Pn}\sum_{a<b}
                      \sum_{i\in\mathcal I}\left(d_i^{(a,b)}(N)\right)^2
                \right]^{1/2},
   \qquad
   D_\infty(N)=\max_{a<b}\max_{i\in\mathcal I}|d_i^{(a,b)}(N)|.

For two codes, these are simply the RMS and maximum absolute values of their
relative-difference field. The maximum uses absolute values so that positive and
negative discrepancies are treated equally.

Zero references and invalid values
----------------------------------

Exactly zero reference bins have undefined pointwise relative errors and are
excluded from both summaries using the same fixed mask at every sampling level.
Report the included and excluded bin counts.
Inspect zero-reference bins in the solution/difference plots: exclusion from a
relative metric does not establish agreement there, and unexpected nonzero
responses must not be hidden by the mask.

Small nonzero reference values remain included; do not introduce a denominator
floor or a sampling-dependent cutoff. Their large relative fluctuations can
make the maximum metric noisy, especially in weakly populated transient bins.
Interpret both summaries together with the solution plots and tally uncertainties.

A nonfinite reference or an included nonfinite estimate makes the comparison
invalid; do not drop such bins silently.
If no nonzero reference bins remain, report the relative metrics as unavailable,
not zero. Any deliberately restricted comparison range must be defined by the
case and kept fixed throughout the study.

Interpreting convergence
------------------------

With fixed bins and reference, ordinary Monte Carlo sampling gives an expected
:math:`N^{-1/2}` trend when statistical uncertainty dominates.
The maximum is especially sensitive to the noisiest bin, so it need not decrease
monotonically in individual realizations.
Reference uncertainty and discretization or modeling differences can limit the
range over which inverse-square-root decay is observed.
Neither norm is itself a confidence interval or a validation pass/fail threshold.

.. _published-metric-conventions:

Existing processor and figure conventions
-----------------------------------------

Existing MC/DC-VVP processors use several different normalizations.
In the formulas below, :math:`n` is the number of entries processed by the
particular helper; only code-to-code processing explicitly masks zero references.
Do not interpret their figures as already using the standardized RMS convention:

- The analytical fixed-source helper ``rerror`` computes :math:`\|r\|_2/n`;
  ``rerror_max`` computes :math:`\max_i|r_i|` without a shared zero-reference mask.
- The analytical helpers ``error`` and ``error_max``, used by several transient
  cases, compute :math:`\|u-u^*\|_2/(\sqrt n\,\|u^*\|_2)` and
  :math:`\max_i(u_i-u_i^*)/\max_i u_i^*`, respectively.
  The latter is signed and is not an infinity norm.
- The eigenvalue helper ``relative_error`` computes the global relative Euclidean
  error :math:`\|u-u^*\|_2/\|u^*\|_2`, rather than an RMS of pointwise relative
  errors.
- The code-to-code processor averages squared **unscaled** pairwise 2-norms
  and reports the maximum absolute pairwise relative difference over nonzero
  reference bins. Its 2-norm is :math:`\sqrt n\,D_2` under the definitions above;
  its all-zero-reference case currently returns zeros.

Changing a fixed normalization changes the ordinate but not the convergence
slope for a fixed bin set. Changing a mask or replacing a signed maximum with an
absolute maximum can also change the observed curve.
Processors must be aligned and affected figures regenerated before they can be
labeled with the standardized definitions. Retain the documented convention for
previously published evidence rather than silently reinterpreting it.
