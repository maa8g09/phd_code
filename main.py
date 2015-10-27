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

import os
import subprocess as sb
from datetime import datetime

def main(directory, fourdarray, geom, Re, kx, kz, c, A, i, output_directory, approximate_soln, outputPhysicalASC, outputSpectralASC, plotting):
    """
    INPUTS:
     directory:  directory where the flow solution is kept (downloaded from 
                 channelflow.org), this is where the ASCII and GEOM files are 
                 kept.
    fourdarray:  a 4D array used to plot slices from channel flow solutions
             N:  resolution in y axis
            Re:  Reynolds number
            kx:  vector of streamwise wavenumbers
            kz:  vector of spanwise wavenumbers
             c:  wave speed
             A:  vector of prescribed amplitudes
             i:  iterator index (for when you're varying amplitudes)
             d:  when generating your own modes, this is the output directory
    
    OUTPUTS:
        Depends...
        
    """
    
    ut.printStart()
    np.set_printoptions(precision = 4)

    # Measure the amount of time it takes to execute this routine.
    startTime = datetime.now()
    

    data = {}
    
    if directory == '':
        #===========================================================================================
        #===========================================================================================
        # OWN MODES
        #===========================================================================================
        #===========================================================================================
        ut.printSectionHeader()
        ut.printSectionTitle('We will be generating our own modes!')
        
        
        #===========================================================================================
        # Resolvent formulation
        generated_flowField = rf.resolvent_analysis(geom, Re, kx, kz, c, A, data, fourdarray)
        
        
        #===========================================================================================
        # Plotting
        if plotting:
            up.plot2D_modes(generated_flowField, fourdarray, True)
        
        
        #===========================================================================================
        # Output Routine
        directory = output_directory
        directory = ut.makeSolutionDirectory(generated_flowField, directory, geom['m'], Re, kx, kz, c, A, i)
                
        if outputPhysicalASC:
            ut.writeASCIIfile(generated_flowField, directory)
            ut.writeGEOMfile(generated_flowField, directory)

        
        if outputSpectralASC:
            directory = directory + '/spectral_construct'
            
            if not os.path.exists(directory):
                os.makedirs(directory)

            ut.writeSpectralASCIIfile(generated_flowField, directory)
            ut.writeSpectralGEOMfile(generated_flowField, directory)
        
        
    else:
        #===========================================================================================
        #===========================================================================================
        # APPROXIMATION
        #===========================================================================================
        #===========================================================================================
        ut.printSectionHeader()
        ut.printSectionTitle('We will be reading in a channelflow solution')
        
        
        #===========================================================================================
        # Read solution
        data = rc.main_read_construct(directory, fourdarray)
        
        
        #===========================================================================================
        # Plotting
        if data['flowField']['is_physical'] == True:
            up.plot2D(data['velslice'])

        #===========================================================================================
        # Resolvent Formulation
        if approximate_soln:
    
            rank = 90
            
            aprxmtd_flowField = rf.resolvent_approximation(data, rank)
            
            
            #=======================================================================================
            # Output Routine
            if outputPhysicalASC:
                directory = directory + '/approximated_physical'
                
                if not os.path.exists(directory):
                    os.makedirs(directory)
                    
                ut.writeASCIIfile(aprxmtd_flowField, directory)
                ut.writeGEOMfile(aprxmtd_flowField, directory)
            
            
            
            if outputSpectralASC:
                directory = directory + '/approximated_spectral'
                
                if not os.path.exists(directory):
                    os.makedirs(directory)
                    
                ut.writeSpectralASCIIfile(aprxmtd_flowField, directory)
                ut.writeGEOMfile(aprxmtd_flowField, directory)



    ut.printSectionHeader()
    ut.printSectionTitle('Calculation Time')
    print('   ', datetime.now() - startTime, '\n')

    return 0








