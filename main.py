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

import os
import subprocess as sb


startTimeLarge = datetime.now()


fdary = [0, 'all', 'all', 16] # XY
#fdary = [0, 0, 'all', 'all'] # YZ

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


    # Check to see if there are any 'unewt' files in the directory...
#    files = [fi for fi in os.listdir(direc) if os.path.isfile(os.path.join(direc,fi))]
#    
#    unewt='unewt'
#    
#    # Now we loop through the files in the directory to find if 
#    # there are any unewt files
#    for fi in files:
#        if unewt in str(fi):
#            # now we make a folder with the same name
#            new_directory=directory+'/'+fi[:-3]
#            if not os.path.exists(new_directory):
#                os.makedirs(new_directory)
#                #copy file into the directory
#                
#                sb.call(['cd',new_directory])
    
    
    
    if directory == '':
        ut.printSectionHeader()
        ut.printSectionTitle('We will be generating our own modes!')
        modesOnly = True
        
        
        # Resolvent Formulation
        generated_flowField = rf.main_resolvent_analysis(N, Re, kx, kz, c, A, modesOnly, data, fourdarray)
        up.plot2D_modes(generated_flowField, fourdarray, True)
        
#        
#        generated_flowField['read']=False
#        if perturbed:
#            perturbedField = ut.perturbFlowField(generated_flowField)
#            up.plot2D_modes(perturbedField, fourdarray, True)
#        
#        directory = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude8'
#        directory = ut.makeSolutionDirectory(generated_flowField, directory, n, re, kx, kz, c, A, i)
#        ut.writeASCIIfile(generated_flowField, directory)
#        ut.writeGEOMfile(generated_flowField, directory)
#        
#        
        
        
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
#    direct = '/Users/arslan/Documents/phd/code/channelflow_solns/nagata1'
#    direct = '/Users/arslan/Documents/work/wavenums/wavepacket_008_4modes_(-0.082096307666+0j)/un8'
#    direct = '/Users/arslan/Documents/work/wavenums-copy/wavepacket_1-0.082096307666/un8'
#    direct = '/Users/arslan/Documents/work/wavenums-copy/wavepacket_3-0.0886455539366/un8'
#    direct = '/Users/arslan/Documents/work/wavenums-copy/wavepacket_4-0.0909427472186/un5'
#    direct = '/Users/arslan/Documents/work/wavenums-copy/wavepacket_5-0.091175715568/unewt08'
    direct = '/Users/arslan/Documents/work/wavenums-copy/wavepacket_6-0.0916434441829/unewt02'
    
elif linux:
#    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude/triplet_+-1_+-1_0.6667_0.1/unewt8'
#    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude4/wavepacket_003_2modes_0.00177827941004j/ubest'
#    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude4/wavepacket_000_2modes_0.01j/unewt2'
#    direct='/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude4/wavepacket_000_2modes_0.01j/channelFlow/u47'
#    direct='/home/arslan/Desktop/testing'
#    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude4/wavepacket_011_4modes_(-4.5e-12+0j)/data'
#    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude2/triplet_case_020_6_+-6_0.6667_3.23745754282e-06/data'
#    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude4/wavepacket_000_4modes_(-0.45+0j)/data/u_0010'
#    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude4/wavepacket_002_2modes_0.00695192796178/unewt0'
#    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude8/wavepacket_005_4modes_(-0.0916434441829+0j)/unewt0'
    direct = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude8/wavepacket_005_4modes_(-0.0916434441829+0j)/ubest'
    
else:
    direct = ''
    




Anam's flat heating system:
Ideal heating settings.
input at 4
output at 1


when you want output, put it on 4,

When leaving teh house, put the input on, then upon returning, turn on output nd turn off input. 


simple inut output system, wherre you havea  "tank" of energy you expel at the end of the recycling process. 



input in on the right
output is on the left






n = 37
re = 400



# Ideal response mode
# K1
kx1a = 6.0
kx1b = 6.0
a1a = 1.0j

kz1a = 6.0
kz1b = -6.0
a1b = 1.0j

# K2
kx2a = 1.0
kx2b = 1.0
a2a = -4.5

kz2a = 6.0
kz2b = -6.0
a2b = -4.5


kx3a = 7.0
kx3b = 7.0
a3a = 0.83j

kz3a = 12.0
kz3b = -12.0
a3b = 0.83j












kx1a = 6.0
kx1b = 6.0
a1a = 1.0j

kz1a = 6.0
kz1b = -6.0
a1b = 1.0j

# K2
kx2a = 0.8
kx2b = 0.8
a2a = -4.5

kz2a = 5.0
kz2b = -5.0
a2b = -4.5

# K3
kx3a = 0.3
kx3b = 0.3
a3a = 0.83j

kz3a = 1.0
kz3b = -1.0
a3b = 0.83j





kx = np.array([kx2a, kx2b, kx1a, kx1b, kx3a, kx3b])
kz = np.array([kz2a, kz2b, kz1a, kz1b, kz3a, kz3b])
amplitudes = np.array([a2a, a2b, a1a, a1b, a3a, a3b])


#
### Test K
#kx1_a = 1.0
#kx1_b = 1.0
#a1_a = -4.5
#
#kz1_a = 6.0
#kz1_b = -6.0
#a1_b = -4.5
#
#kx = np.array([kx1_a, kx1_b])
#kz = np.array([kz1_a, kz1_b])
#amplitudes = np.array([a1_a, a1_b])


c = 2.0/3.0



if linux:
    ampl_weights = np.array([1.0])
else:
    ampl_weights = np.logspace(-1.65, -1.75, 10)
    ampl_weights = np.logspace(-1.68, -1.70, 10)
    ##
    ampl_weights = np.array([1e-2])



for i in range(0, len(ampl_weights)):
    ampl = amplitudes*ampl_weights[i]
    main(direct, fdary, n, re, kx, kz, c, ampl, i)


#main(direct, fdary, n, re, kx, kz, c, amplitudes, 0)

print('   ', datetime.now() - startTimeLarge, '\n')
