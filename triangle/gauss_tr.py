# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 17:47:44 2017

@author: taro
"""

import math, numpy, scipy.optimize
import matplotlib.pyplot as plt


OpenfileName='eodam-triangle_30_30_30_50'



def gaussian(x, A, mean, sigma):
    gauss = A/math.sqrt(2.0*math.pi)/sigma * numpy.exp(-((x-mean)/sigma)**2/2)
    return(gauss)
    
def residuals_with_error(param_fit, y, x, yerr):
    A, mean, sigma, base = param_fit
    err = (y - (gaussian(x, A, mean, sigma)+base)) / yerr
    return(err)
    
def main():
        
    count=0
    n=0
    size=140
    cutsize=15
    
    numbers =numpy.zeros(size)  
    fitgauss =numpy.zeros(size)
    
    xmin, xmax, nx = 0.0, size, size
    x = numpy.arange(xmin, xmax, (xmax-xmin)/nx)
    
    for line in open(OpenfileName+'.dat', 'r'):
        items = line.split()
        if( cutsize<count and count<size+cutsize+1 ):
            numbers[n] = abs(float(items[1]) )
            n=n+1
        count=count+1
        
    argmax=numpy.argmax(numbers)
    Vmax=numpy.max(numbers)
        
    base0 = 1.0
    err_frac = Vmax*0.01
    
    y_meas = numbers + err_frac*base0*numpy.random.randn(len(x))
    y_err  = err_frac*base0*numpy.random.randn(len(x))
        
    param0 = [10, argmax, 30, 1.0] # initial guess: χ^2, initial paramator
    param_output = scipy.optimize.leastsq(residuals_with_error, param0,
                                          args=(y_meas, x, y_err), 
                                          full_output=True)
    param_result = param_output[0] # fitted parameters
    covar_result = param_output[1] # covariant matrix
    
    def gauss_result_for_plot(x, param):
        result = gaussian(x, param[0], param[1], param[2]) + param[3]
        return result    

    fitgauss=gauss_result_for_plot(x, param_result)
    
    print numpy.max(fitgauss) 
    
    print("Amplitude: %10.5f +/- %10.5f"
          % (numpy.max(fitgauss), numpy.sqrt(covar_result[0][0])))
          #% (param_result[0], numpy.sqrt(covar_result[0][0])))
    print("Center   : %10.5f +/- %10.5f"
          % (param_result[1], numpy.sqrt(covar_result[1][1])))
    print("Sigma    : %10.5f +/- %10.5f"
          % (param_result[2], numpy.sqrt(covar_result[2][2])))
#    print("Baseline : %10.5f +/- %10.5f"
#          % (param_result[3], numpy.sqrt(covar_result[3][3])))

    
    plt.plot(x, gauss_result_for_plot(x, param_result), # plot result
             x, numbers)  # plot true curve
    plt.title('Least-squares fit to noisy data')
    plt.rcParams['font.family'] = 'IPAGothic'
    plt.title(str(OpenfileName) )
    plt.legend(['Fit', 'Original-data'])
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.ylim(-0.01, 0.05)
    plt.text(20, -.005, r" Ampritude : %10.5f +/- %10.5f"
          % (numpy.max(fitgauss), numpy.sqrt(covar_result[0][0]))+
          "\n Center       : %10.5f +/- %10.5f"
          % (param_result[1], numpy.sqrt(covar_result[1][1]))+
          "\n Sigma       : %10.5f +/- %10.5f"
          % (param_result[2], numpy.sqrt(covar_result[2][2])))
          
    plt.savefig('fit_gaussian '+str(OpenfileName)+'.png', dpi=150)
    plt.show()
    
    
    
main()