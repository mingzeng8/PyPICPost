import scipy.optimize as opt
import numpy as np

def Lorentzian(x, center, height, half_width):
    f = height/(1+((x-center)/half_width)**2)
    return f

def FitLorentzian(x, y, initial_guess=None):
    '''Fit Lorentzian distribution and return the parameters'''
    if initial_guess is None:
        height_guess = y.max()
        center_guess = x[np.argmax(y)]
        half_width_guess = 0.05*center_guess
        initial_guess = [center_guess, height_guess, half_width_guess]
    popt, pcov = opt.curve_fit(Lorentzian, x, y, p0=initial_guess, bounds=([x.min(), y.max()*0.1, 0.,], [x.max(), y.max()*2., x.max()-x.min()]))
    return popt, pcov
