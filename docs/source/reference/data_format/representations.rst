========================
Reusable Representations
========================

``CrossSection``, ``Data``, ``Distribution``, and ``Interpolation`` describe reusable group layouts.
These type names identify the layouts expanded at individual fields.
``Data`` and ``Distribution`` select a representation through a required scalar string attribute ``type`` on the containing group.

CrossSection
------------

A reaction's ``xs`` group stores only the values needed on its projectile's shared ``xs-energy-grid``.
Let N be the shared grid length and o the offset.

.. list-table:: CrossSection fields
   :header-rows: 1

   * - Field
     - Type / shape / unit
     - Meaning
   * - ``offset``
     - integer
     - Zero-based index o of the first stored cross section on the shared grid.
   * - ``value``
     - float [N - o], barn
     - Cross-section values corresponding to ``xs-energy-grid[o:]``.
   * - ``interpolation``
     - Interpolation
     - Interpolation over the stored values.

The offset is nonnegative and leaves at least two stored grid points.
Cross sections are nonnegative and are zero at shared-grid points below the offset.
Interpolation boundaries index the local ``value`` array, starting at zero, rather than the full shared grid.
CrossSection does not carry its own grid or a ``type`` attribute.

Data
----

Data represents a scalar-valued function of an input coordinate.
The containing field defines the input and output quantities, such as incident energy in eV and neutron multiplicity.
Each representation below lists all of its datasets.

.. list-table:: Data representations
   :header-rows: 1
   :widths: 20 40 40

   * - ``type`` attribute
     - Fields
     - Meaning
   * - ``constant``
     - ``value``: float
     - The same value for every input.
   * - ``tabulated``
     - ``grid``: float [K]; ``value``: float [K]; ``interpolation``: Interpolation
     - Values on a nondecreasing input grid with K >= 2 and at least two distinct coordinates.
   * - ``polynomial``
     - ``coefficient``: float [L]
     - Coefficients in ascending powers of the input, with L >= 1.

A polynomial represents

.. math::

   f(x) = \sum_{l=0}^{L-1} c_l x^l.

The coefficient at index zero is the constant term.
Coefficient units are output units divided by input units raised to the coefficient's index.
For example, a dimensionless multiplicity polynomial of incident energy in eV has coefficient units ``1``, ``1/eV``, ``1/eV^2``, and so on.

Constant and polynomial Data do not have interpolation fields.
Tabulated Data requires an explicit ``interpolation`` group.
No generic extrapolation rule is implied outside a tabulated grid's range.
Repeated coordinates store the left and right limits of a discontinuity in that order.
Interpolation applies between distinct coordinates, not across the zero-width interval.

Interpolation
-------------

An ``interpolation`` group accompanies the tabulated representation whose intervals it describes.
Let K be that representation's number of tabulated points and R its number of interpolation regions.

.. list-table:: Interpolation fields
   :header-rows: 1

   * - Field
     - Type / shape
     - Meaning
   * - ``laws``
     - string [R]
     - One interpolation law per region.
   * - ``boundaries``
     - integer [R]
     - Zero-based upper point index of each region.

The arrays have equal nonzero length.
Boundaries strictly increase, the first boundary is at least one, and the last equals K - 1.
The first region begins at point zero; subsequent regions begin at the preceding region's upper point.
Thus every interval belongs to exactly one region.
For example, boundaries ``[2, 5]`` apply the first law to intervals 0–1 and 1–2, then the second law to intervals 2–3, 3–4, and 4–5.

For neighboring points :math:`(x_0,y_0)` and :math:`(x_1,y_1)`, define :math:`t=(x-x_0)/(x_1-x_0)` and :math:`s=\ln(x/x_0)/\ln(x_1/x_0)`.

.. list-table:: Interpolation laws
   :header-rows: 1
   :widths: 20 45 35

   * - Law
     - Interpolated value
     - Domain condition
   * - ``histogram``
     - :math:`y=y_0` on the interval below :math:`x_1`
     - None beyond an increasing grid.
   * - ``linear``
     - :math:`y=(1-t)y_0+t y_1`
     - None beyond an increasing grid.
   * - ``semilog-x``
     - :math:`y=(1-s)y_0+s y_1`
     - Positive x coordinates.
   * - ``semilog-y``
     - :math:`y=y_0(y_1/y_0)^t`
     - Positive y values.
   * - ``log``
     - :math:`y=y_0(y_1/y_0)^s`
     - Positive x coordinates and y values.

Logarithms operate on dimensionless ratios.
Each tabulated representation specifies its interpolation laws and any restrictions on them.

Distribution
------------

A Distribution describes sampled quantities rather than a deterministic scalar function.
An energy-only distribution samples outgoing energy in eV, an angular distribution samples a dimensionless cosine, and a joint distribution supplies both.
The reaction's reference frame defines the frame of those quantities.
An angular cosine is measured with respect to the incident direction; these layouts do not store an azimuth-dependent law.

Tabulated density
^^^^^^^^^^^^^^^^^

``type = tabulated`` stores a one-dimensional continuous probability density.

.. list-table:: Tabulated distribution fields
   :header-rows: 1

   * - Field
     - Type / shape
     - Meaning
   * - ``value``
     - float [K]
     - Nondecreasing sampled coordinates with positive support width, K >= 2.
   * - ``pdf``
     - float [K]
     - Nonnegative probability density at those coordinates.
   * - ``interpolation``
     - Interpolation
     - Density interpolation between coordinates.

``value`` carries the sampled coordinate's unit, and ``pdf`` carries its inverse.
Both omit ``unit`` for a dimensionless cosine distribution.
Cosines lie between -1 and 1.
The density integrates to one over its tabulated support, within the precision of the source data and the specified interpolation law.
The terminal density of a histogram does not define a further interval.
Repeated coordinates describe a jump in density; their zero-width interval carries no probability mass.
The CDF is obtained by integrating the PDF using those laws.

Isotropic and composed distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Representations
   :header-rows: 1
   :widths: 28 40 32

   * - ``type`` attribute
     - Fields
     - Meaning
   * - ``isotropic``
     - None
     - Angular cosine uniformly distributed on [-1, 1], with uniform azimuth.
   * - ``independent-angle-energy``
     - ``angular-cosine-distribution``: Distribution; ``energy-distribution``: Distribution
     - Independent outgoing angle and energy, conditioned on incident energy where applicable.
   * - ``elastic-angle-energy``
     - ``angular-cosine-distribution``: Distribution
     - Angle distribution with energy determined by elastic two-body kinematics.
   * - ``mixture``
     - ``components/{i}/probability``: Data; ``components/{i}/distribution``: Distribution
     - Select a component by its probability, then sample its distribution.

Mixture probabilities depend on incident energy in eV, are nonnegative, and sum to one at each incident energy in their common domain.
All components sample the same physical quantities in the same frame.

Incident-energy dependence
^^^^^^^^^^^^^^^^^^^^^^^^^^

``type = incident-energy-dependent`` contains the following fields.

.. list-table:: Fields
   :header-rows: 1

   * - Field
     - Type / shape / unit
     - Meaning
   * - ``energy``
     - float [I], eV
     - Nondecreasing incident energies, I >= 2.
   * - ``interpolation``
     - Interpolation
     - Incident-energy dependence; histogram or linear only.
   * - ``distributions/{i}``
     - Distribution
     - One distribution for each incident-energy point.

All child distributions describe the same sampled quantity.
Repeated incident energies describe the distributions on either side of a discontinuity.
Histogram interpolation uses the lower incident-energy distribution within an interval.
Linear interpolation uses the neighboring distributions with weights 1 - t and t.
For outgoing-energy tables with different supports, use unit-base interpolation: linearly interpolate the lower and upper support limits, and map the sampled energy from the selected table's support into that interpolated support.
For angular tables the cosine coordinate is unchanged.

Analytic energy spectra
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table:: Representations
   :header-rows: 1
   :widths: 22 38 40

   * - ``type`` attribute
     - Fields
     - Meaning
   * - ``level-scattering``
     - ``C1``: float, eV; ``C2``: float
     - Outgoing energy :math:`E'=C_2(E-C_1)` for incident energy E.
   * - ``maxwellian``
     - ``temperature``: Data, eV; ``restriction-energy``: float, eV
     - Spectrum proportional to :math:`\sqrt{E'}\exp[-E'/T(E)]`.
   * - ``evaporation``
     - ``temperature``: Data, eV; ``restriction-energy``: float, eV
     - Spectrum proportional to :math:`E'\exp[-E'/T(E)]`.

For Maxwellian and evaporation spectra, T(E) is the positive spectral temperature, evaluated as Data at incident energy E in eV.
The normalized spectrum is restricted to :math:`0\le E'\le E-U`, where U is ``restriction-energy`` and E - U must be positive.
Spectral temperature is expressed in eV, and target temperature in kelvin.

Correlated energy-angle distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``type = kalbach-mann`` and ``type = energy-angle-tabulated`` both contain ``energy`` (float [I], eV) and ``tables/{i}``, one table per strictly increasing incident-energy point.
Their incident-energy dependence is linear with unit-base interpolation of outgoing energy.

Each Kalbach–Mann table contains:

.. list-table:: Kalbach–Mann table fields
   :header-rows: 1

   * - Field
     - Type / unit
     - Meaning
   * - ``outgoing-energy``
     - Tabulated Distribution, eV
     - Marginal outgoing-energy density.
   * - ``precompound-factor``
     - Data
     - Dimensionless precompound fraction r versus outgoing energy in eV.
   * - ``angular-slope``
     - Data
     - Dimensionless angular slope a versus outgoing energy in eV.

The parameters describe the conditional cosine density

.. math::

   p(\mu\mid E,E') = \frac{a}{2\sinh a}
   \left[\cosh(a\mu)+r\sinh(a\mu)\right].

The limit a = 0 is isotropic.
The precompound fraction is between zero and one, and the angular slope is nonnegative.
When tabulated, both parameter grids match the table's outgoing-energy coordinates.

Each energy-angle-tabulated table contains ``outgoing-energy`` (tabulated Distribution, eV) and ``cosine-distributions/{j}`` (tabulated Distribution, dimensionless).
There is one cosine distribution for every outgoing-energy point j.
The outgoing-energy distribution is the marginal density, and each cosine distribution is conditional on that outgoing energy and incident-energy table.

N-body energy distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^

``type = n-body`` contains ``value`` (float [K], eV), ``pdf`` (float [K], 1/eV), and ``interpolation`` (Interpolation).
It stores a tabulated N-body energy density.
Its array and normalization constraints are the same as those of a tabulated energy distribution.
