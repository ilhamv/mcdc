.. _project_vvp_validation:

==========
Validation
==========

Validation assesses whether MC/DC predictions agree with experimental observations
within the relevant measurement, model, and nuclear-data uncertainties.
Each benchmark page identifies the experiment, computational representation,
quantity of interest, comparison range, and limitations needed to interpret its
results. Experimental agreement must be read in the context of those assumptions.

SINBAD benchmarks
-----------------

The :doc:`SINBAD suite <sinbad>` uses experimental shielding benchmarks from the
OECD/NEA Shielding Integral Benchmark Archive and Database.
Its executable campaign contains the OKTAVIAN Si-60 neutron-leakage case and the
four FNG/TUD SiC neutron-spectrum cases.
Other model inputs describe activation, heating, and photon observables whose
complete calculation requires additional response or coupled-particle capabilities.
Those inputs are documented separately and are not enabled validation comparisons.

- :doc:`oktavian_si60`: neutron and gamma leakage configurations, source modeling,
  and the defensible neutron comparison range.
- :doc:`fng_tud_sic`: four detector-depth configurations and their neutron and
  photon spectrum observables.
- :doc:`fng_sic`: activation and TLD-heating configurations, response requirements,
  and intermediate neutron tallies.
- :doc:`rfnc_photon_compounds`: photon leakage from spherical and hemispherical
  compound samples and the required coupled neutron-photon calculation.

Interpreting experimental comparisons
-------------------------------------

A calculated-to-experimental ratio (C/E) describes agreement for a specified
observable, binning, and normalization. Counting uncertainty alone does not capture
source normalization, detector response, geometry, composition, nuclear data, or
other model uncertainties.
The SINBAD packages used here do not supply a complete experimental covariance
matrix and prescribed uncertainties for all modeling assumptions.
Consequently, binwise comparisons and reported uncertainty bands do not establish
a covariance-aware goodness-of-fit test or a formal validation pass/fail criterion.

Use :doc:`sinbad` to configure, run, and process these cases, including access to
the licensed experimental tables. See :doc:`../contributing`
when adding an experimental benchmark or extending its observable.

.. toctree::
   :maxdepth: 1

   sinbad
   oktavian_si60
   fng_tud_sic
   fng_sic
   rfnc_photon_compounds
