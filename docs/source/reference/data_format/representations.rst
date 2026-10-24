========================
Reusable Representations
========================

``CrossSection``, ``Data``, ``Distribution``, and ``Interpolation`` describe reusable group layouts.
These type names identify the layouts expanded at individual fields.
The reusable representations are designed to be as consistent as possible with their counterparts in MC/DC's ``mcdc.object_`` and ``mcdc.transport`` modules.
``Data`` and ``Distribution`` select a representation through a required scalar string attribute ``type`` on the containing group.
See :doc:`../ace_distributions` for the mapping from ACE laws to these representations and their object and transport implementations.

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

Probability table
^^^^^^^^^^^^^^^^^

``type = table`` stores a continuous PDF, discrete probability masses, or both.

.. list-table:: Table distribution fields
   :header-rows: 1

   * - Field
     - Type / shape
     - Meaning
   * - ``value``
     - float [K]
     - Continuous coordinates, required with ``pdf``; nondecreasing with positive support width and K >= 2.
   * - ``pdf``
     - float [K]
     - Nonnegative continuous density, required with ``value``.
   * - ``interpolation``
     - Interpolation
     - Density interpolation, required only for the continuous component.
   * - ``discrete-value``
     - float [D]
     - Discrete outcomes, required with ``pmf``, D >= 1.
   * - ``pmf``
     - float [D]
     - Nonnegative probability masses, required with ``discrete-value``.

At least one coordinate/probability pair is required.
Coordinates and probabilities are finite, and probability values are nonnegative.
The PMF sum plus the PDF integral differs from one by at most :math:`10^{-6}`.
The integral follows the specified density interpolation laws.
The two components carry their combined probabilities, rather than each being normalized to one.
Discrete outcomes may lie inside or outside the continuous support and are sampled exactly.
``discrete-value`` carries the sampled coordinate's unit, and ``pmf`` is dimensionless.
``value`` carries the sampled coordinate's unit, and ``pdf`` carries its inverse.
Both omit ``unit`` for a dimensionless cosine distribution.
Cosines lie between -1 and 1.
A continuous-only table has a PDF integral of one; a discrete-only table has a PMF sum of one.
The terminal density of a histogram does not define a further interval.
Repeated coordinates describe a jump in density; their zero-width interval carries no probability mass.
CDFs and cumulative probability masses are derived from the stored PDF and PMF.

Isotropic distribution
^^^^^^^^^^^^^^^^^^^^^^

``type = isotropic`` describes an angular cosine uniformly distributed on [-1, 1], with uniform azimuth.
It has no datasets.

Multi-table distributions
^^^^^^^^^^^^^^^^^^^^^^^^^

``type = multi-table`` is a specialized conditional continuous-PDF representation with its own packed arrays.
It does not contain reusable ``table`` Distribution objects and does not support discrete probability masses.
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
All PDFs describe the same sampled quantity, share its physical unit, and individually integrate to one.
Repeated grid coordinates describe the distributions on either side of a discontinuity.
See :ref:`ace-incident-energy-sampling` for table selection and outgoing-energy support scaling.

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
Each table has at least two points and at least one interpolation region, following the coordinate, density, and interpolation rules of a continuous-only probability table.
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
     - Deterministic level-scattering energy parameters.
   * - ``maxwellian``
     - ``temperature``: Data, eV; ``restriction-energy``: float, eV
     - Maxwellian spectral temperature and energy restriction.
   * - ``evaporation``
     - ``temperature``: Data, eV; ``restriction-energy``: float, eV
     - Evaporation spectral temperature and energy restriction.

For Maxwellian and evaporation spectra, T(E) is the positive spectral temperature, evaluated as Data at incident energy E in eV.
Spectral temperature is expressed in eV, and target temperature in kelvin.
See :ref:`ace-law-3`, :ref:`ace-law-7`, and :ref:`ace-law-9` for the energy laws and their sampling.

Correlated energy-angle distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``type = kalbach-mann`` and ``type = energy-angle-table`` both contain ``energy`` (float [I], eV) and ``tables/outgoing-energy`` (packed PDF tables).
There is one outgoing-energy table per strictly increasing incident-energy point.
These representations have inherent linear incident-energy dependence and do not store an outer interpolation group.

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

The precompound fraction is between zero and one, and the angular slope is nonnegative.
Each parameter is evaluated as Data at outgoing energy in eV.
When represented as tables, the parameters may use independent grids and interpolation regions, distinct from one another and from the outgoing-energy density.
See :ref:`ace-law-44` for the conditional angular law and parameter sampling.

Energy-angle-table distributions additionally contain ``tables/cosine-distributions`` (packed PDF tables, dimensionless).
There is one cosine table for every point in the concatenated ``tables/outgoing-energy/value`` array, in the same order.
The outgoing-energy tables supply marginal densities; each cosine table is conditional on its corresponding outgoing energy and incident-energy table.
The two collections' offsets identify both levels of conditioning.

See :ref:`ace-law-61` for selection and sampling of the conditional cosine tables.

N-body phase-space distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``type = n-body`` stores a dimensionless energy-fraction PDF and the final-state mass parameters.

.. list-table:: N-body fields
   :header-rows: 1

   * - Field
     - Type / shape
     - Meaning
   * - ``number-particles``
     - integer
     - Number of bodies N in the final-state phase space, N >= 3.
   * - ``total-mass-ratio``
     - float
     - Total final-state mass divided by the emitted particle mass, R > 1.
   * - ``value``
     - float [K]
     - Normalized outgoing-energy fractions x in [0, 1].
   * - ``pdf``
     - float [K]
     - Dimensionless density of the energy fraction.
   * - ``interpolation``
     - Interpolation
     - Interpolation of the energy-fraction density.

The PDF integrates to one under its interpolation rules, following the continuous probability-table constraints.
The density supports all five Interpolation laws and multiple interpolation regions.
N-body reactions require ``reference-frame = COM`` and omit a separate angular-cosine distribution.
See :ref:`ace-law-66` for the phase-space model, physical energy scale, and angular sampling.
