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


def main(directory, fourdarray, N, Re, kx, kz, c, A):
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
             A:  amplitude
    
    
    OUTPUTS:
        An ASCII file with th resolvent modes
        
    """
    
    ut.printStart()
    np.set_printoptions(precision = 4)

    # Measure the amount of time it takes to execute this routine.
    startTime = datetime.now()
    

    data = {}
    
    modesOnly = False
    perturbed = False


    
    
    
    if directory == '':
        ut.printSectionHeader()
        ut.printSectionTitle('We will be generating our own modes!')
        modesOnly = True
        
        
        # Resolvent Formulation
        generated_flowField = rf.main_resolvent_analysis(N, Re, kx, kz, c, A, modesOnly, data, fourdarray)
        up.plot2D_modes(generated_flowField, fourdarray, True)
        
        
        generated_flowField['read']=False
        if perturbed:
            perturbedField = ut.perturbFlowField(generated_flowField)
            up.plot2D_modes(perturbedField, fourdarray, True)
        
        directory = '/home/arslan/Documents'
        directory = ut.makeSolutionDirectory(generated_flowField, directory)
        ut.writeASCIIfile(generated_flowField, directory)
        ut.writeGEOMfile(generated_flowField, directory)
        
        
        
        
    else:
        ut.printSectionHeader()
        ut.printSectionTitle('We will be reading in a channelflow solution')
        
        # Read in the channelflow solution and store flowfield
        #
        # The data variable should be a class. See how classes work in python
        # and think about how you could implement them into your code:
        # http://www.jesshamrick.com/2011/05/18/an-introduction-to-classes-and-inheritance-in-python/
        #
        data = rc.main_read_construct(directory, fourdarray)
        
        if data['flowField']['is_physical'] == True:
            up.plot2D(data['velslice'])

        if perturbed:
            data['read']=True
            perturbedField = ut.perturbFlowField(data)
            ut.writeASCIIfile_general(perturbedField, directory)
            perturbed_slice = rc.get_vel_slice()##
            up.plot2D(perturbed_slice)

        # Resolvent Formulation
        generated_flowField = rf.main_resolvent_analysis(N, Re, kx, kz, c, A, modesOnly, data, fourdarray)
        up.plot2D_modes(generated_flowField, fourdarray, True)
        
        directory = '/home/arslan/Documents'
        directory = ut.makeSolutionDirectory(generated_flowField, directory)
        ut.writeASCIIfile(generated_flowField, directory)
        ut.writeGEOMfile(generated_flowField, directory)


    ut.printSectionHeader()
    ut.printSectionTitle('Calculation Time')
    print('   ', datetime.now() - startTime, '\n')

    return





mac = False
linux = False

if mac:
    direct = '/Users/arslan/Documents/phd/code/channelflow_solns/nagata1'
    
elif linux:
    direct = '/home/arslan/Documents/work/channelflow-related/tutorial_invariant_solns/test01'
    direct = '/home/arslan/Documents/work/channelflow-related/database_solns/HKW/equilibria/eq7'
else:
    direct = ''
    
    
n = 32
re = 400
kx = np.arange(1,2)
kz = np.arange(-1,2)
c = 2./3.
fdary = [0, 0, 'all', 'all']
A = 1.0e0

main(direct, fdary, n, re, kx, kz, c, A)

