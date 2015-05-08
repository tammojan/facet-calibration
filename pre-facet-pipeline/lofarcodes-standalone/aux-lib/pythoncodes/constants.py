from math import *


# speed of light
c = 3E8 #(m/s)

#plancks constant
h = 6.626068E-34 #(m^(2)kg/s)

#boltzman constant
k_b = 1.3806503E-23 #(m^(2)kgs^(-2)K^(-1))

#Density params
Omega_m = 0.3
Omega_r = 1E-5
Omega_lambda = 0.7
Omega_k = 0.0

# Hubble constant
H_0 = 72 # (km/s)/Mpc

#Define various angle conversion factors (multiply to undertake operation)
arcsec2deg=1.0/3600
arcmin2deg=1.0/60
deg2rad=pi/180
deg2arcsec = 1.0/arcsec2deg
rad2deg=180.0/pi
arcmin2rad=arcmin2deg*deg2rad
arcsec2rad=arcsec2deg*deg2rad
rad2arcmin=1.0/arcmin2rad
rad2arcsec=1.0/arcsec2rad
steradians2degsquared = (180.0/pi)**2.0
degsquared2steradians = 1.0/steradians2degsquared

# Jansky

jy2si = 1E-26
si2jy = 1/(1*10**(-26))

pc2metres = 3.08E16
metres2pc = 1/pc2metres
boltzman_constant = 1.380648E-23 # J/K
electron_volt = 1.60217657E-19 # J
proton_mass = 1.672E-27
electron_mass = 9.1E-31
magnetic_nu0 =  1.256E-6
