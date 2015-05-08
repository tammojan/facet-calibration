# facet-calibration
Facet calibration: scripts updated at the Facet Calibration Workshop in Leiden, 20-24 April 2015.

## Main scripts

Set up the facets:
* doDDE_makefacets_v4.py

Run the direction dependent calibration:
* doDDE_v19_a2256.py - Main script. It requires a configuration file with some variables defined (example a2256.parset).
* doDDE_v21_a2256.py - Wideband version of the main script. *Under heavy development*.

Reimage the field when the full facet calibration is completed:
* doDDE_reimage_v1_bootes_cep3.py - Version that uses the CASA imager.
* doDDE_reimage_wsclean.py - Version that uses wsclean.
