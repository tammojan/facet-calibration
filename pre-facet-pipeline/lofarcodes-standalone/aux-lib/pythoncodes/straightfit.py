from scipy.optimize import leastsq
import scipy
import pylab
import numpy as np


def model(t, coeffs):
    return t*coeffs[0] + coeffs[1]
def residuals(coeffs, y, t):
    return y - model(t, coeffs)

x0 = [1, 1]
dataset1 = [1,3,8,10,16,12,18,22,28,40]

t = np.arange(len(dataset1))
x, flag = scipy.optimize.leastsq(residuals, x0, args=(dataset1, t))

resid = residuals(x,dataset1,t)
chisq = sum(resid*resid)

print model(t,x)
print 'Chisquared of fit to the straight line is %s' %chisq
pylab.plot(dataset1)
pylab.plot(t,dataset1,t,model(t,x))
pylab.show()

