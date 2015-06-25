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


import numpy as np
import utils as ut
import read_construct as rc
import utils_plots as up
import resolvent_formulation as rf
from datetime import datetime


def main(directory, fourdarray, N, Re, kx, kz, c):
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
    perturbed = True
    
    if directory == '':
        ut.printSectionHeader()
        ut.printSectionTitle('We will be generating our own modes!')
        modesOnly = True
        
    else:
        ut.printSectionHeader()
        ut.printSectionTitle('We will be reading in a channelflow solution')
        
        # Read in the channelflow solution and store flowfield
        data = rc.main_read_construct(directory, fourdarray)
        
        if data['flowField']['is_physical'] == True:
            up.plot2D(data['velslice'])

        if perturbed:
            data['read']=True
            perturbedField = ut.perturbFlowField(data)
            ut.writeASCIIfile_general(perturbedField, directory)
            perturbed_slice = rc.get_vel_slice()##
            up.plot2D(perturbedField)

    # Resolvent Formulation
    generated_flowField = rf.main_resolvent_analysis(N, Re, kx, kz, c, modesOnly, data, fourdarray)
    up.plot2D_modes(generated_flowField, fourdarray, True)
    
    generated_flowField['read']=False
    if perturbed:
        perturbedField = ut.perturbFlowField(generated_flowField)
        up.plot2D_modes(perturbedField, fourdarray, True)
    
    # Write ASCII and geom file for channelflow
    if directory == '':
        directory = '/home/arslan/Desktop'
        directory = ut.makeSolutionDirectory(perturbedField, directory)
        ut.writeASCIIfile(perturbedField, directory)
        ut.writeGEOMfile(perturbedField, directory)
    
    ut.printSectionHeader()
    ut.printSectionTitle('Calculation Time')
    print '   ', datetime.now() - startTime, '\n'

    return





mac = False
linux = True

if mac:
    direct = '/Users/arslan/Documents/phd/code/channelflow_solns/nagata1'
    
elif linux:
#    direct = '/home/arslan/Documents/phd/code/channelflow-1.4.2/solutions/equilibria/nagata1'
    direct ='/home/arslan/Documents/work/channelflow-related/database_solns/HKW/equilibria/eq4'
    direct ='home/arslan/ubest'
#    direct ='/home/arslan/Documents/work/channelflow-related/database_solns/HKW/equilibria/eq4/nonlinear_solver/perturbed'
#    direct ='/home/arslan/Documents/work/channelflow-related/chflow_wavepackets/laminar'
    
else:
    direct = ''
    
    
n = 42
re = 400
kx = np.array([6])
#kz = np.array([-10,10])
kz = np.arange(-1,2) #(-8, to 8)
c = 2.0 / 3.0
fdary = [0, 'all', 'all', 0]

main(direct, fdary, n, re, kx, kz, c)

