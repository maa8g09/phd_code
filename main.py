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


fdary = [0, 'all', 'all', 0] # XY # 0 24 36
fdary = [0, 0, 'all', 'all'] # YZ


def main(directory, fourdarray, N, Re, kx, kz, c, A, i, d):
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
#    perturbed = False
    outputPhysical = False

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
#        up.plot2D_modes(generated_flowField, fourdarray, True)
        
#        
#        generated_flowField['read']=False
#        if perturbed:
#            perturbedField = ut.perturbFlowField(generated_flowField)
#            up.plot2D_modes(perturbedField, fourdarray, True)
#        
        #directory = '/home/arslan/Documents/work/channelflow-related/edge_state_varying_amplitude8'
        directory = d
        directory = ut.makeSolutionDirectory(generated_flowField, directory, n, re, kx, kz, c, A, i)
        
        
        if outputPhysical:
            ut.writeASCIIfile(generated_flowField, directory)
            ut.writeGEOMfile(generated_flowField, directory)
        
        else:
            directory = directory + '/spectral_construct'
            if not os.path.exists(directory):
                os.makedirs(directory)
            ut.writeSpectralASCIIfile(generated_flowField, directory)
            ut.writeSpectralGEOMfile(generated_flowField, directory)
        
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
            u = data['flowField']['physical'][0,:,:,:]
            v = data['flowField']['physical'][1,:,:,:]
            w = data['flowField']['physical'][2,:,:,:]
            up.plot2D(data['velslice'])

###        else:
#        # Resolvent Formulation
#        rank = 1
#        
#        generated_flowField = rf.main_gibson_soln(data, rank)
#        
##        up.plot2D_modes(generated_flowField, fourdarray, True)
#        
#        if outputPhysical:
#            ut.writeASCIIfile(generated_flowField, directory)
#            ut.writeGEOMfile(generated_flowField, directory)
#        
#        else:
#            directory = directory + '/reconstructed'
#            
#            if not os.path.exists(directory):
#                os.makedirs(directory)
#                
#            ut.writeSpectralASCIIfile(generated_flowField, directory)
#            ut.writeGEOMfile(generated_flowField, directory)



    ut.printSectionHeader()
    ut.printSectionTitle('Calculation Time')
    print('   ', datetime.now() - startTime, '\n')

    return





mac = False
linux = True

if mac:
    direct = '/Users/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ1'
    
    
elif linux:
    direct = '/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ1'
#    direct = '/home/arslan/Desktop/asciibinarytests/wavepacket_000_4modes_(-0.045+0j)/fromFF2ASC'
    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/from_phys'
    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/from_spec'
    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3'
    #Original Executables
#    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/original_executables/original_ff2ascii'
#    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/original_executables/original_ff2ascii/original_ascii2ff'
#    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/new_executables_spec_only'
#    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/findsoln_rank-1/unewt01'
    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/findsoln_rank-1/pre-unewt0'
#    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/findsoln_original/pre-unewt0'
#    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/findsoln_rank-full/pre-unewt0'
#    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/TEMP'
#    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ3/TEMP/back2ff'
    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ5/fromASCII/fromFF'

    
    
else:
    direct = ''
    










#=================== K1
kx1a = 6.0
kx1b = 6.0

kz1a = 6.0
kz1b = -6.0

a1a = 1.0j
a1b = 1.0j



#=================== K2
kx2a = 1.0
kx2b = 1.0

kz2a = 6.0
kz2b = -6.0

a2a = -4.5
a2b = -4.5



#=================== K3
kx3a = 7.0
kx3b = 7.0

kz3a = 12.0
kz3b = -12.0

a3a = 0.83j
a3b = 0.83j




#=================== K1-STAR
kx1a_star = 6.0
kx1b_star = 6.0

kz1a_star = 8.0
kz1b_star = -8.0

a1a_star = 1.0j
a1b_star = 1.0j



#=================== K3-STAR
kx3a_star = 7.0
kx3b_star = 7.0

kz3a_star = 14.0
kz3b_star = -14.0

a3a_star = 0.83j
a3b_star = 0.83j











#=================== K4
kx4a = 0.3
kx4b = 0.3

kz4a = 3.0
kz4b = -3.0

a4a = 0.3
a4b = 0.3

#=================== K5
kx5a = 1.5
kx5b = 1.5

kz5a = 4.0
kz5b = -4.0

a5a = 1.0
a5b = 1.0

#=================== K6
kx6a = 2.1
kx6b = 2.1

kz6a = 5.0
kz6b = -5.0

a6a = 3.0
a6b = 3.0

#=================== K7
kx7a = 1.0
kx7b = 1.0

kz7a = 6.0
kz7b = -6.0

a7a = 2.0
a7b = 2.0













## KA
#kx = np.array([kx1a, kx1b])
#kz = np.array([kz1a, kz1b])
#amplitudes = np.array([a1a, a1b])
#
#kx = np.array([kx2a, kx2b])
#kz = np.array([kz2a, kz2b])
#amplitudes = np.array([a2a, a2b])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re00400/kA/ampls01'






## KB
#kx = np.array([kx2a, kx2b, kx1a, kx1b])
#kz = np.array([kz2a, kz2b, kz1a, kz1b])
#amplitudes = np.array([a2a, a2b, a1a, a1b])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re01800/kB/ampls01'
##d = '/home/arslan/Desktop/asciibinarytests'





## KC
#kx = np.array([kx2a, kx2b, kx1a, kx1b, kx3a, kx3b])
#kz = np.array([kz2a, kz2b, kz1a, kz1b, kz3a, kz3b])
#amplitudes = np.array([a2a, a2b, a1a, a1b, a3a, a3b])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re01800/kC/ampls01'





## KD
#kx = np.array([kx2a, kx2b, kx1a_star, kx1b_star, kx3a_star, kx3b_star])
#kz = np.array([kz2a, kz2b, kz1a_star, kz1b_star, kz3a_star, kz3b_star])
#amplitudes = np.array([a2a, a2b, a1a_star, a1b_star, a3a_star, a3b_star])
#d = '/home/arslan/Documents/work/channelflow-related/set01/Re01800/kD/ampls01'





# KE
kx = np.array([kx4a, kx4b, kx5a, kx5b, kx6a, kx6b, kx7a, kx7b])
kz = np.array([kz4a, kz4b, kz5a, kz5b, kz6a, kz6b, kz7a, kz7b])
amplitudes = np.array([a4a, a4b, a5a, a5b, a6a, a6b, a7a, a7b])
d = '/home/arslan/Documents/work/channelflow-related/set01/Re01800/kE/ampls01'





n = 37
re = 400
c = 2.0/3.0



if linux:
    ampl_weights = np.array([1.0])
else:
    ## Step 1
    ampl_weights = np.logspace(1.0, -3.0, 10)



    ## Step 2
    ## B C D
#    ampl_weights = np.logspace(-0.7 , -1.2 , 10)
    ## E
#    ampl_weights = np.logspace(-1.2, -1.6, 10)



    ## Step 3
    ## B C D
#    ampl_weights = np.logspace(-0.97, -1.09, 10)
    ## E
#    ampl_weights = np.logspace(-1.47, -1.51, 10)



    ## Step 4
    ## B C D
#    ampl_weights = np.logspace(-0.99, -1.01, 15)
    ## E
#    ampl_weights = np.logspace(-1.497, -1.501, 5)






#    ampl_weights = np.array([1e-2])
    
    



for i in range(0, len(ampl_weights)):
    ampl = amplitudes*ampl_weights[i]
    main(direct, fdary, n, re, kx, kz, c, ampl, i, d)


#main(direct, fdary, n, re, kx, kz, c, amplitudes, 0)

print('   ', datetime.now() - startTimeLarge, '\n')




