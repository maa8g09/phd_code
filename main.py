"""
RESOLVENT FORMULATION

Resolvent formulation and projection code. This is the main routine that calls 
all of the functions and files that perform the whole routine of generating 
resolvent modes and then projecting them onto a channel flow solution.


Author details:
    Muhammad Arslan Ahmed
    maa8g09@soton.ac.uk
    
    Aerodynamics and Flight Mechanics Research Group
    Faculty of Engineering and the Environment
    University of Southampton
"""

mac = True

import numpy as np
import utils as ut
import read_construct as rc
import utils_plots as up
import resolvent_formulation as rf
from datetime import datetime


def main(directory, verbose, fourdarray, N, Re, kx, kz, c):
    """
    INPUTS:
     directory:  directory where the flow solution is kept (downloaded from 
                 channelflow.org), this is where the ASCII and GEOM files are 
                 kept.
       verbose:  how verbose the output in the python shell should be, setting 
                 it to:
                   - false: only essentials are output
                   - true: every detail is printed
    fourdarray:  a 4D array used to plot slices from channel flow solutions
             N:  resolution in y axis
            Re:  Reynolds number
            kx:  vector of streamwise wavenumbers
            kz:  vector of spanwise wavenumbers
             c:  phase speed
    
    
    OUTPUTS:
        An ASCII file with th resolvent modes
        
    """
    
    ut.printStart()
    np.set_printoptions(precision = 4)

    # Measure the amount of time it takes to execute this routine.
    startTime = datetime.now()
    
    data = {}
    modesOnly = False
    
    if directory == '':
        ut.printSectionHeader()
        ut.printSectionTitle('We will be generating our own modes!')
        modesOnly = True
        
    else:
        ut.printSectionHeader()
        ut.printSectionTitle('We will be reading in a channelflow solution')
        
        # Read in the channelflow solution and store flowfield
        data = rc.main_read_construct(directory, fourdarray, verbose)
        
        if data['flowField']['is_physical'] == True:
            up.plot2D(data['velslice'])

    # Resolvent Formulation
    generated_flowField = rf.main_resolvent_analysis(N, Re, kx, kz, c, modesOnly, data, fourdarray)
    up.plot2D_modes(generated_flowField, fourdarray)
    
    ut.printSectionHeader()
    ut.printSectionTitle('Calculation Time')
    print '   ', datetime.now() - startTime, '\n'

    return








if mac:
# MAC 
#    direct = '/Users/arslan/Dropbox/python/'
    direct = '/Users/arslan/Documents/phd/code/channelflow_solns/nagata1'
else:    
# DEBIAN
#    direct = '/home/arslan/python/'
#    direct = '/home/arslan/Documents/phd/code/channelflow-1.4.2/solutions/travelling_waves/viswanath2007/tw2'
#    direct = '/home/arslan/Documents/phd/code/channelflow-1.4.2/solutions/travelling_waves/halcrow2008b/tw3'
    direct = '/home/arslan/Documents/phd/code/channelflow-1.4.2/solutions/equilibria/nagata1'
#    direct = '/home/arslan/Documents/channelflow-1.4.2/solutions/equilibria/p47p18'





direct = ''
n = 35
re = 400
kx = np.arange(-15, 17)
kz = np.arange(0, 17)
U_centreline = 20
c = np.linspace(0.12, 1.0, 50) * U_centreline


v = True


fdary = [0, 0, 'all', 'all']


main(direct, v, fdary, n, re, kx, kz, c)

