import scipy.optimize as opt
import numpy as np
import pylab as plt


#define model function and pass independant variables x and y as a list
def twoD_Gaussian(xy, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    x=xy[:(xy.size)//2]
    y=xy[(xy.size)//2:]
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    #print(g.shape)
    #return g.ravel()
    return g

def Fit2DGauss(x, y, data, initial_guess):
    '''Fit 2D Gaussian distribution and return the parameters'''
    xy=np.append(x,y)
    popt, pcov = opt.curve_fit(twoD_Gaussian, xy, data.ravel(), p0=initial_guess)
    return popt, pcov

def twoD_Gaussian_simple(xy, amplitude, xo, yo, sigma, offset):
    ''''twoD_Gaussian_simple is a simple version of twoD_Gaussian with sigma_x=sigma_y=sigma and theta=0'''
    x=xy[:(xy.size)//2]
    y=xy[(xy.size)//2:]
    xo = float(xo)
    yo = float(yo)    
    a = 0.5/(sigma**2)
    g = offset + amplitude*np.exp( - a*((x-xo)**2 + (y-yo)**2))
    #print(g.shape)
    #return g.ravel()
    return g

def Fit2DGauss_simple(x, y, data, initial_guess, bounds=((-np.inf,-np.inf,-np.inf,0.,-np.inf),(np.inf,np.inf,np.inf,np.inf,np.inf))):
    '''Fit2DGauss_simple is a simple version of Fit2DGauss with sigma_x=sigma_y=sigma and theta=0'''
    xy=np.append(x,y)
    popt, pcov = opt.curve_fit(twoD_Gaussian_simple, xy, data.ravel(), p0=initial_guess, bounds=bounds)
    return popt, pcov

def Gaussian_simple(r, amplitude, sigma):
    '''Gaussian_simple is a simple version of Gaussian function with center=0 and offset=0'''
    a = 0.5/(sigma*sigma)
    g = amplitude*np.exp( - r*r*a)
    return g

def Fit_Gauss_simple(r, data, initial_guess, bounds=((-np.inf,0.),(np.inf,np.inf))):
    '''Fit_Gauss_simple is a simple version of Gaussian fit with center=0 and offset=0'''
    popt, pcov = opt.curve_fit(Gaussian_simple, r, data, p0=initial_guess, bounds=bounds)
    return popt, pcov
