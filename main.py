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

startTimeLarge = datetime.now()

def main(directory, fourdarray, N, Re, kx, kz, c, A, i):
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
        
        directory = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude_single_mode'
        directory = ut.makeSolutionDirectory(generated_flowField, directory, n, re, kx, kz, c, A, i)
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

#        if perturbed:
#            data['read']=True
#            perturbedField = ut.perturbFlowField(data)
#            ut.writeASCIIfile_general(perturbedField, directory)
#            perturbed_slice = rc.get_vel_slice()##
#            up.plot2D(perturbed_slice)

#        # Resolvent Formulation
#        generated_flowField = rf.main_resolvent_analysis(N, Re, kx, kz, c, A, modesOnly, data, fourdarray)
#        up.plot2D_modes(generated_flowField, fourdarray, True)
#        
#        directory = '/home/arslan/Documents'
#        directory = ut.makeSolutionDirectory(generated_flowField, directory)
#        ut.writeASCIIfile(generated_flowField, directory)
#        ut.writeGEOMfile(generated_flowField, directory)


    ut.printSectionHeader()
    ut.printSectionTitle('Calculation Time')
    print('   ', datetime.now() - startTime, '\n')

    return





mac = False
linux = False

if mac:
    direct = '/Users/arslan/Documents/phd/code/channelflow_solns/nagata1'
    
elif linux:
    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude/triplet_+-1_+-1_0.6667_0.1/unewt8'
    
else:
    direct = ''
    

n = 100
re = 400

# Ideal Packet
kx_a = np.arange(7) #6
kx_b = np.arange(2) #1
kx_c = np.arange(8)
kx = np.array([kx_a, kx_b, kx_c])

kz_a = np.arange(-6,7)
kz_b = np.arange(-6,7)
kz_c = np.arange(-12,13)
kz = np.array([kz_a, kz_b, kz_c])

amplitude_a = 1.0j
amplitude_b = -4.5
amplitude_c = 0.83j
amplitudes = np.array([amplitude_a, amplitude_b, amplitude_c])




kx_a = np.arange(2)
kx = np.array([kx_a])

kz_a = np.arange(-6,7)
kz = np.array([kz_a])



#amplitude_a = 1.0
#amplitude_b = 1.0
#amplitudes = np.array([amplitude_a, amplitude_b])

c = 2.0/3.0
fdary = [0, 'all', 'all', 0]
ampl_weights = np.logspace(-0.1, -2.0, num=5)
#ampl_weights = np.array([1.0])

for i in range(0, len(ampl_weights)):
    main(direct, fdary, n, re, kx, kz, c, amplitudes*ampl_weights[i], i)

#main(direct, fdary, n, re, kx, kz, c, amplitudes, 0)

print('   ', datetime.now() - startTimeLarge, '\n')
