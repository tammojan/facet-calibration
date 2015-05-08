import numpy as np
from scipy.optimize import leastsq
import scipy
import pylab


#-----------------------------------------------------------

def prime_factors(n):
    """Find the prime factors of a number
    """
    factors=[]
    d=2
    while(d*d<n):
        while(n>1):            
            while n%d==0:
                factors.append(d)
                n=n/d
            d+=1
    return factors[-1]

#-----------------------------------------------------------
    
def close_higher_bsmooth_number(n,B):
    """Finds the closest higher B-smooth number (good numbers for FFTs) http://en.wikipedia.org/wi$
        """
    fivesmooths = np.array([1])
    i = np.max([n-300,5])
    while np.max(fivesmooths) < n:
        if prime_factors(i) <=B and (i/2.0 == int(i)/2): #Image must have an even number of pixels
            fivesmooths = np.append(fivesmooths,i)
        i+=1
    print fivesmooths
    fivesmooths=np.sort(fivesmooths)
    closest_smooth_index = np.where(abs(fivesmooths-n)==np.min(abs(fivesmooths-n)))
    if len(fivesmooths[closest_smooth_index[0]]) == 2:
	closest_smooth = fivesmooths[closest_smooth_index[0][1]]
        return closest_smooth
    print closest_smooth_index
    if fivesmooths[closest_smooth_index[0]] < n:
        closest_smooth = fivesmooths[closest_smooth_index[0]+1]
    else:
        closest_smooth = fivesmooths[closest_smooth_index[0]]
    return closest_smooth[0]

#-----------------------------------------------------------

def find_chi_squared(values,expected):
    """Given a set of values and a set of expected values this
    will return the chi-squared value
    """

    chisquared = 0.0
    
    for i in range(0,len(values)):
        chisquared += ((values[i]-expected[i])**2.0)/abs(expected[i])

    return chisquared

#-----------------------------------------------------------

def find_chi_squared_werrors(values,expected,error):
    """Given a set of values and a set of expected values this
    will return the chi-squared value
    """

    chisquared = 0.0
    
    for i in range(0,len(values)):
        chisquared += ((values[i]-expected[i])**2.0)/error[i]**2.0

    return chisquared



#-----------------------------------------------------------

def standard_dev(list):
    """ 
    Given a list this returns the standard deviation
    """
    sum = 0.0
    sum2 = 0.0
    for element in list:
        sum = float(line) + sum
        sum2 = float(line)**2 + sum2
        i = i + 1

    mean = sum/i
    stan_dev = (abs((sum2/i)-(mean**2)))**0.5
    rms = sqrt(sum2/i)

    return stan_dev

#-----------------------------------------------------------

def mean(list):
    """ 
    Given a list this routine returns the mean
    """

    i=0.0
    sum = 0.0
    for element in list:
        sum += float(element)
        i += 1.0
    mean = sum/i
    return mean
#-----------------------------------------------------------

def rms(list):
    """
    Given a list this routine returns the rms
    """
    i =0.0
    sum = 0.0
    sum2 = 0.0
    for element in list:
        sum = float(element) + sum
        sum2 = float(element)**2 + sum2
        i = i + 1.0
     
    rms = sqrt(sum2/i)

    return rms

#------------------------------------------------------------

def two_dimensional_gaussian(A,x,x0,sigma_x,y,y0,sigma_y):
    """
    Given the parameters for a two dimensional gaussian this returns the value of the Gaussian
    """
    value = A*exp(-((x-x0)**2.0/(2*sigma_x**2.0) + (y-y0)**2.0/(2*sigma_y**2.0)))
    
    return value

#------------------------------------------------------------

def model_gaussian(t, coeffs):
    return coeffs[0] + coeffs[1] * np.exp( - ((t-coeffs[2])**2.0/(2*coeffs[3]**2)))

#------------------------------------------------------------
def residuals_gaussian(coeffs, y, t):
    return y - model_gaussian(t, coeffs)

#------------------------------------------------------------
def fit_gaussian_histogram(pixelvals,plotting):

    fitnumbers,cellsizes = np.histogram(pixelvals,100)
    sigmaguess = np.std(pixelvals)/(abs(cellsizes[1]-cellsizes[0]))
    x0 = [0.0,max(fitnumbers),np.where(fitnumbers==max(fitnumbers))[0][0],sigmaguess] #Offset amp, amp, x-offset, sigma
    t = np.arange(len(fitnumbers))
    x, flag = scipy.optimize.leastsq(residuals_gaussian, x0, args=(fitnumbers, t))

    if plotting == 'y':
        pylab.plot(fitnumbers)
        pylab.plot(t,fitnumbers,t,model_gaussian(t,x))
        pylab.show()
        pylab.close()
        pylab.cla()

    print 'Sigma is %s'%(x[3]*abs(cellsizes[1]-cellsizes[0]))
    
    return (x[3]*abs(cellsizes[1]-cellsizes[0]))
