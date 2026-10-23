========================
Reusable Representations
========================

``CrossSection``, ``Data``, ``Distribution``, and ``Interpolation`` describe reusable group layouts.
These type names identify the layouts expanded at individual fields.
The reusable representations are designed to be as consistent as possible with their counterparts in MC/DC's ``mcdc.object_`` and ``mcdc.transport`` modules.
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
   * - ``table``
     - ``grid``: float [K]; ``value``: float [K]; ``interpolation``: Interpolation
     - Values on a nondecreasing input grid with K >= 2 and at least two distinct coordinates.
   * - ``polynomial``
     - ``coefficients``: float [L]
     - Coefficients in ascending powers of the input, with L >= 1.

A polynomial represents

.. math::

   f(x) = \sum_{l=0}^{L-1} c_l x^l.

The coefficient at index zero is the constant term.
A constant uses ``type = polynomial`` with a single coefficient.
Coefficient units are output units divided by input units raised to the coefficient's index.
For example, a dimensionless multiplicity polynomial of incident energy in eV has coefficient units ``1``, ``1/eV``, ``1/eV^2``, and so on.

Polynomial Data does not have interpolation fields.
Data tables require an explicit ``interpolation`` group.
No generic extrapolation rule is implied outside a table's grid range.
Repeated coordinates store the left and right limits of a discontinuity in that order.
Interpolation applies between distinct coordinates, not across the zero-width interval.

Interpolation
-------------

An ``interpolation`` group accompanies the table representation whose intervals it describes.
Let K be that representation's number of table points and R its number of interpolation regions.

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
     - Exclusive upper-point index of each region.

The arrays have equal nonzero length.
Boundaries strictly increase, the first boundary is at least two, and the last equals K.
The first region begins at point zero; subsequent regions begin at the preceding boundary minus one.
Thus every interval belongs to exactly one region.
For example, boundaries ``[3, 6]`` apply the first law to intervals 0–1 and 1–2, then the second law to intervals 2–3, 3–4, and 4–5.

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
Each table representation specifies its interpolation laws and any restrictions on them.

Distribution
------------

A Distribution describes sampled quantities rather than a deterministic scalar function.
An energy-only distribution samples outgoing energy in eV, an angular distribution samples a dimensionless cosine, and a joint distribution supplies both.
The reaction's reference frame defines the frame of those quantities.
An angular cosine is measured with respect to the incident direction; these layouts do not store an azimuth-dependent law.

Density table
^^^^^^^^^^^^^

``type = table`` stores a one-dimensional continuous probability density.

.. list-table:: Table distribution fields
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
The density integrates to one over its stored support, within the precision of the source data and the specified interpolation law.
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

Multi-table distributions
^^^^^^^^^^^^^^^^^^^^^^^^^

``type = multi-table`` stores a family of continuous distributions represented by tables on another grid.
It contains the following fields.

.. list-table:: Multi-table fields
   :header-rows: 1

   * - Field
     - Type / shape
     - Meaning
   * - ``grid``
     - float [I]
     - Nondecreasing indexing coordinates, I >= 2, with a ``unit`` attribute when dimensional.
   * - ``interpolation``
     - Interpolation
     - Interpolation between neighboring tables; histogram or linear only.
   * - ``tables``
     - Packed PDF tables
     - One conditional distribution per grid point, in grid order.

For incident-energy dependence, ``grid`` has unit ``eV``.
A dimensionless indexing grid omits the unit attribute.
All tables describe the same sampled quantity and share its physical unit.
Repeated grid coordinates describe the distributions on either side of a discontinuity.
Histogram interpolation uses the lower-grid distribution within an interval.
Linear interpolation uses the neighboring distributions with weights 1 - t and t.
For outgoing-energy tables with different supports, use unit-base interpolation: linearly interpolate the lower and upper support limits, and map the sampled energy from the selected table's support into that interpolated support.
For angular tables the cosine coordinate is unchanged.

Packed PDF tables
^^^^^^^^^^^^^^^^^

A collection of PDF tables uses concatenated arrays.
This layout is shared by multi-table and correlated distributions; its containing representation defines how the tables are indexed.
For T tables, P total points, and R total interpolation regions, the collection contains:

.. list-table:: Packed PDF table fields
   :header-rows: 1

   * - Field
     - Type / shape
     - Meaning
   * - ``offsets``
     - integer [T + 1]
     - Start of each table in ``value`` and ``pdf``, followed by the terminal point count P.
   * - ``value``
     - float [P]
     - Concatenated sampled coordinates, carrying their physical unit when dimensional.
   * - ``pdf``
     - float [P]
     - Concatenated densities, carrying the inverse coordinate unit when dimensional.
   * - ``interpolation/offsets``
     - integer [T + 1]
     - Start of each table's regions in ``laws`` and ``boundaries``, followed by R.
   * - ``interpolation/laws``
     - string [R]
     - Concatenated interpolation laws for the individual tables.
   * - ``interpolation/boundaries``
     - integer [R]
     - Exclusive upper-point indices relative to the start of each individual table.

Both offset arrays start at zero and are strictly increasing.
Table i occupies the half-open slice ``offsets[i]:offsets[i + 1]``; its interpolation regions use the corresponding slice of ``interpolation/offsets``.
Each table has at least two points and at least one interpolation region, following the same coordinate, density, normalization, and interpolation rules as ``type = table``.
Coordinates may restart at each table boundary; the concatenated array need not be globally monotonic.
No ``type`` attribute is attached to the collection itself.
CDFs are derived from the stored densities and interpolation rules.

For example, tables with 3, 2, and 4 points have ``offsets = [0, 3, 5, 9]``.
The last table occupies points 5 through 8, and its final interpolation boundary is 4, regardless of its position in the packed array.

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

``type = kalbach-mann`` and ``type = energy-angle-table`` both contain ``energy`` (float [I], eV) and ``tables/outgoing-energy`` (packed PDF tables).
There is one outgoing-energy table per strictly increasing incident-energy point.
Their incident-energy dependence is linear with unit-base interpolation of outgoing energy.

Kalbach–Mann additionally contains:

.. list-table:: Kalbach–Mann parameter fields
   :header-rows: 1

   * - Field
     - Type
     - Meaning
   * - ``tables/precompound-factor/{i}``
     - Data
     - Dimensionless precompound fraction r versus outgoing energy in eV, for incident-energy index i.
   * - ``tables/angular-slope/{i}``
     - Data
     - Dimensionless angular slope a versus outgoing energy in eV, for incident-energy index i.

The parameters describe the conditional cosine density

.. math::

   p(\mu\mid E,E') = \frac{a}{2\sinh a}
   \left[\cosh(a\mu)+r\sinh(a\mu)\right].

The limit a = 0 is isotropic.
The precompound fraction is between zero and one, and the angular slope is nonnegative.
When represented as tables, both parameter grids match the table's outgoing-energy coordinates.

Energy-angle-table distributions additionally contain ``tables/cosine-distributions`` (packed PDF tables, dimensionless).
There is one cosine table for every point in the concatenated ``tables/outgoing-energy/value`` array, in the same order.
The outgoing-energy tables supply marginal densities; each cosine table is conditional on its corresponding outgoing energy and incident-energy table.
The two collections' offsets identify both levels of conditioning.

N-body energy distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^

``type = n-body`` contains ``value`` (float [K], eV), ``pdf`` (float [K], 1/eV), and ``interpolation`` (Interpolation).
It stores an N-body energy-density table.
Its array and normalization constraints are the same as those of an energy distribution table.
The N-body phase-space model assumes isotropic emission in the center-of-mass frame, independently of outgoing energy in that frame.
Emission in the laboratory frame is generally anisotropic after the kinematic transformation.
See the `NJOY manual's N-body phase-space definition <https://github.com/njoy/NJOY2016-manual/blob/master/acex.tex>`_.
