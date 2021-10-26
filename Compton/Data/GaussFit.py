# Utility functions
import numpy as np
from scipy.optimize import curve_fit

def gauss(x, A, x0, std_dev, bg):
    return A * np.exp(-(x - x0) ** 2 / (2 * std_dev ** 2)) + bg

def gauss_fit(x, y):
    peak_index = np.argmax(y)
    peak_value = x[peak_index]
    poissonError = np.sqrt(y)
    std_dev = np.sqrt(sum((x - peak_value) ** 2) / (sum(y)-1)) #sample standard deviation
    print(max(y), peak_value, peak_index, std_dev)
    popt, pcov = curve_fit(gauss, x, y, p0=[ max(y), peak_value, std_dev, 0], sigma = poissonError, maxfev=5000)
    return popt, pcov
    
def normal_gauss(x,A,x0,std_dev): 
    return A/np.sqrt(2*np.pi*std_dev**2)*np.exp(-(x-x0)**2/(2*std_dev**2))

def gauss_fit_to_spectrum(input_x,input_y,x1,x2, angle):
    # Select peak to fit
    xdata = input_x[x1:x2]
    ydata = input_y[x1:x2]
    [A, x0, std_dev],pcov = gauss_fit(xdata, ydata)
    # offset, amplitude, center, sigma
    
    [A_err, x0_err, std_dev_err] = np.sqrt(np.diag(pcov))
    print( [A_err, x0_err, std_dev_err])
    # Plot input spectra with correct scaling
    plt.plot(input_x*m+b,input_y, label='Subtracted Spectra')
    
    # Plot selected peak to fit
    plt.plot(xdata*m+b, ydata, label='Selected Peak')
    
    # Plot fitted gaussian to peak
    plt.plot(input_x*m+b, gauss(input_x, A, x0, std_dev), '--r', label='Fit')
    # plt.annotate(f'Mean={round(x0*m+b,2)}keV, $\sigma$={round(sigma,2)}, Amplitude={round(A,2)}counts/s',(0.2,0.5), xycoords='figure fraction')
    plt.legend()
    plt.title(angle + ' deg')
    plt.xlabel('Energy Levels (KeV)')
    plt.ylabel('Counts/s')
    plt.show()
    return np.array([x0*m+b,A,std_dev,x0_err,A_err,std_dev_err])