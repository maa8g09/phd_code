mac = True

import numpy as np
from datetime import datetime
import utils
import read_construct as rc
import utils_plots as up
import resolvent_formulation as rf

"""
Resolvent formulation and projection code. This is the main routine that calls 
all of the functions and files that perform the whole routine of generating 
resolvent modes and then projecting them onto a channel flow solution.

    Author:  Muhammad Arslan Ahmed
     email:  maa8g09@soton.ac.uk

University of Southampton
"""

def main(directory, verbose, fourdarray, N, Re, kx, kz, c):
    """
    INPUTS:
      directory:  directory where the flow solution is kept (downloaded from channelflow.org),
                  this is where the ASCII and GEOM files are kept.
        verbose:  how verbose the output in the python shell should be, setting it to:
                    - false: only essentials are output
                    -  true: every detail is printed
     fourdarray:  a 4D array used to plot slices from channel flow solutions
              N:  resolution in y axis
             Re:  Reynolds number
             kx:  vector of streamwise wavenumbers
             kz:  vector of spanwise wavenumbers
              c:  phase speed
    
    
    OUTPUTS:
        An ASCII file with th resolvent modes
    """
    np.set_printoptions(precision=4)
    
    
    # Measure the amount of time it takes to execute this routine.
    startTime = datetime.now()
    
    
    # Read in the channelflow solution and store flowfield
#    geo_variables, flow_field, vel_slice = rc.main_read_construct(directory, fourdarray, verbose)
    # I am now goint to update the above code so that a dictionary is produced, 
    # this way, only one variable is passed to functions and it can be extended
    # to include whatever I may need ata later time.
    data = rc.main_read_construct(directory, fourdarray, verbose)
    
    up.plot2D(vel_slice)
    
    # Resolvent Formulation
    # If i want ot have random geometrical variables, I can not pass in flow_field.
#    projected_solution = rf.main_resolvent_analysis(flow_field, geo_variables, fourdarray, N, Re, kx, kz, c, True)
    rf.main_resolvent_analysis(data, fourdarray, N, Re, kx, kz, c, True)
#    up.plot2D_modes(projected_solution)
    
    utils.printSectionHeader()
    utils.printSectionTitle('Calculation Time')
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


n = 8
re = 400
kx = np.array([1])
kz = np.array([0])
U_centreline = 10
c = np.linspace(0.01, 1.0, 1) * U_centreline


v = True


fdary = [0, 0, 'all', 'all']


main(direct, v, fdary, n, re, kx, kz, c)

