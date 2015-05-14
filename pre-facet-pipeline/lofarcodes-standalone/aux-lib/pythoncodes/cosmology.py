import numpy
import scipy
import scipy.integrate as si
import pylab

# These calculations are with w_0 = -1  (equation of state)

def drange(start, stop, step):
	r = start
	while r < stop:
		yield r
		r += step


def _comoving_integral(z, omega_m, omega_k, omega_lambda, h):
  
    e_z = (omega_m * (1+z)**3. +
           omega_k * (1+z)**2. +
           omega_lambda)**0.5

    Mpc_cm = 3.08568025e24 # cm
    Mpc_km = Mpc_cm * 1.0e-5 # km
    H100_s = 100. / Mpc_km # s^-1
    H_0 = h * H100_s
  
    H_z =  H_0 * e_z

    c_light_cm_s = 29979245800. # cm/s
    c_light_Mpc_s = c_light_cm_s / Mpc_cm # Mpc / s
    
    return c_light_Mpc_s / (H_z)


def comoving_distance(z,omega_m, omega_k, omega_lambda, h):
	
	d_co, err =   si.quad(_comoving_integral, 0.0, z, args=(omega_m, omega_k, omega_lambda, h))
	return d_co

def angular_diameter_distance(z):

    # define cosmology
    h = 0.7
    omega_m = 0.3
    omega_k = 0.0
    omega_lambda = 0.7

    comoving_d = comoving_distance(z,omega_m, omega_k, omega_lambda, h)

    angular_d = comoving_d/(1+z)
	

    return angular_d

def luminosity_distance(z):

    # define cosmology
    h = 0.7
    omega_m = 0.3
    omega_k = 0.0
    omega_lambda = 0.7

    comoving_d = comoving_distance(z,omega_m, omega_k, omega_lambda, h)

    lumin_d = comoving_d*(1+z)
	

    return lumin_d


