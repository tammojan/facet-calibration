# facet-calibration
Facet calibration: scripts updated at the Facet Calibration Workshop in Leiden, 20-24 April 2015.

## Main scripts

Set up the facets:
* doDDE_makefacets_v4.py - It requires a configuration file with some variables defined (example doDDE_makefacets_parset.py).

Run the direction dependent calibration:
* doDDE_v19_a2256.py - Old main script. It requires a configuration file with some variables defined (example a2256.parset). It is strongly recommended to switch to doDDE_v21_a2256.py.

* doDDE_v21_a2256.py - Main script. Accepts wideband. It requires a configuration file with some variables defined (example a2256.parset).

Reimage the field when the full facet calibration is completed:
* doDDE_reimage_v1_bootes_cep3.py - Version that uses the CASA imager.
* doDDE_reimage_wsclean.py - Version that uses wsclean.

## Helper scripts

* facet_utilities.py - Framework used to run external code
* selfcalv19_ww_cep3.py - Self-cal based on the global solvers
* verify_subtract_v5.py - Check whether the subtraction was OK
* cal_return_slist.py - Return the sources in a facet
* blanker.py - Blank outside the facet mask
* coordinates_mode.py - WCS utilities

[TBD: complete]

PYTHONPATH should be set to the path of the helper scripts if the main
doDDE... script is not located in the same directory
