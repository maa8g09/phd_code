"""
UTILITIES

This file contains miscellaneous functions, which are used for small tasks such 
as printing errors, messages and checking input validity.


Author details:
    Muhammad Arslan Ahmed
    maa8g09@soton.ac.uk
    
    Aerodynamics and Flight Mechanics Research Group
    Faculty of Engineering and the Environment
    University of Southampton
"""

import os
import sys
import math

def printStart():
    print('##################################################################')
    print('_________ .__                                .__            ')
    print('\_   ___ \|  |__ _____    ____   ____   ____ |  |           ')
    print('/    \  \/|  |  \\__   \  /    \ /    \_/ __ \|  |           ')
    print('\     \___|   Y  \/ __ \|   |  \   |  \  ___/|  |__         ')
    print(' \______  /___|  (____  /___|  /___|  /\___  >____/         ')
    print('        \/     \/     \/     \/     \/     \/               ')
    print('    __________                    .__                      __   ')
    print('    \______   \ ____   __________ |  |___  __ ____   _____/  |_ ')
    print('     |       _// __ \ /  ___/  _ \|  |\  \/ // __ \ /    \   __\.')
    print('     |    |   \  ___/ \___ (  <_> )  |_\   /\  ___/|   |  \  |  ')
    print('     |____|_  /\___  >____  >____/|____/\_/  \___  >___|  /__|  ')
    print('            \/     \/     \/                     \/     \/      ')
    print('Muhammad Arslan Ahmed')
    print('University of Southampton')
    print('##################################################################')
    return


def printSectionHeader():
    print('__________________________________________________________________\n')
    return


def printSectionTitle(str):
    # Print the section headers for each main section of the code output.
    print(' **', str, '\n')
    return


def error(str):
    # Print the error and then exit from the program entirely
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('ERROR:', str)
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('\n\n\n')
    print('  ______    ______   __       __  ________       ')
    print(' /      \  /      \ |  \     /  \|        \      ')
    print('|  $$$$$$\|  $$$$$$\| $$\   /  $$| $$$$$$$$      ')
    print('| $$ __\$$| $$__| $$| $$$\ /  $$$| $$__          ')
    print('| $$|    \| $$    $$| $$$$\  $$$$| $$  \         ')
    print('| $$ \$$$$| $$$$$$$$| $$\$$ $$ $$| $$$$$         ')
    print('| $$__| $$| $$  | $$| $$ \$$$| $$| $$_____       ')
    print(' \$$    $$| $$  | $$| $$  \$ | $$| $$     \      ')
    print('  \$$$$$$  \$$   \$$ \$$      \$$ \$$$$$$$$      ')                          
    print('                                                 ')
    print('  ______   __     __  ________  _______          ')
    print(' /      \ |  \   |  \|        \|       \         ')
    print('|  $$$$$$\| $$   | $$| $$$$$$$$| $$$$$$$\        ')
    print('| $$  | $$| $$   | $$| $$__    | $$__| $$        ')
    print('| $$  | $$ \$$\ /  $$| $$  \   | $$    $$        ')
    print('| $$  | $$  \$$\  $$ | $$$$$   | $$$$$$$\        ')
    print('| $$__/ $$   \$$ $$  | $$_____ | $$  | $$        ')
    print(' \$$    $$    \$$$   | $$     \| $$  | $$        ')
    print('  \$$$$$$      \$     \$$$$$$$$ \$$   \$$        ')
    sys.exit('')
    return
    
    
def message(str):
    print('   ', str, '\n')
    return


def openFile(str):
    """
    INPUTS:
         str:  string of the directory where the file exists.
    OUTPUTS:
           f:  file object
    """
    
    f = open(str, 'r')
    message('Opened the file: ' + str)
    return f
    

def printGeoVars(geom_variables):
    """
    This function prints out the grid dimensions of the channelflow solution.
    """
    print('\n    The geometry variables are:')
    print('      Nx:', geom_variables['Nx'])
    print('      Ny:', geom_variables['Ny'])
    print('      Nz:', geom_variables['Nz'])
    
    return


def checkInputValidity(fourdarray, geom_variables):
    """
    This function checks that the values in the fourdarray are valid to allow
    the plotting of the slice specified. If it is valid the function executes
    as normal and a message is output in the shell to let the user know all's
    well.
    
    If, however, the values in the fourdarray are invalid, the whole routine 
    is stopped and an error message is output which specifies the invalidity
    of the input. 
    
    
    INPUTS:
        fourdarray: the plane co-ordinates of data you want to plot, indexed
                      as (i, nx, ny, nz)
      geom_variables: a dictionary contianing all the geometrical values
      
    OUTPUTS:
               
    """
    printSectionHeader()
    printSectionTitle('Checking validity of 4D array input')
    
    for i in range(1,4):
        if fourdarray[i] != 'all':
            if i == 1:
                if fourdarray[i] >= geom_variables['Nx']:
                    print('\n  ! The 4D array input is invalid')
                    printGeoVars(geom_variables)
                    error('X point given exceeds the maximum grid points available')
            elif i == 2:
                if fourdarray[i] >= geom_variables['Ny']:
                    print('\n  ! The 4D array input is invalid')
                    printGeoVars(geom_variables)
                    error('Y point given exceeds the maximum grid points available')
            else:
                if fourdarray[i] >= geom_variables['Nz']:
                    print('\n  ! The 4D array input is invalid')
                    printGeoVars(geom_variables)
                    error('Z point given exceeds the maximum grid points available')
    
    if fourdarray[0] < 0 or fourdarray[0] >= 3:
        message('Invalid velocity component given, velocity component must be in range 0 to 2.')
        error('Invalid velocity component given!')
    
    message('The 4D array input is valid.')
    return
    
    
def writeASCIIfile(data, directory):
    
    """
    The data variable is a dictionary which has the generated flowField and 
    geometry information. The dictionary can be unpacked and an ASCII file with 
    the flowfield can be written. It should be written with the following indices:
    
    u(nx, ny, nz, i)
    
    The file us written with the title:
    wave_packet_kx_kz_c_A.asc
    
    The wave number triplet followed by the amplitude.
    
    The way that I store the flow field is in the following indexing
    U[i, nx, ny, nz]. Just be weary of this when unpacking the matrix to file.
    
    """

    kx = data['kx']
    kz = data['kz']
    c = data['c']
    A = data['A']
    
    flowField = data['resolvent_flowField']
    
    fileName = "/wave_packet_" + str(kx) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A) + ".asc"
    
    file = open(directory + fileName, "w")
    
    for nx in range(0, data['Nx']):
        for ny in range(0, data['Ny']):
            for nz in range(0, data['Nz']):
                for nd in range(0, data['Nd']):
                    tmp = flowField[nd, nx, ny, nz]
                    tmp = format(tmp, '.16f')
                    file.write(tmp + "\n")

    file.close()
    
    print('Boom, the ASCII is done.')
    
    return 0
    
    
    
def writeGEOMfile(data, directory):
    
    kx = data['kx']
    kz = data['kz']
    c = data['c']
    A = data['A']
    
    fileName = "/wave_packet_" + str(kx) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A) + ".geom"
    
    file = open(directory + fileName, "w")
    
    file.write(str( int(data['Nx']) ) + '\t\t\t\t\t\t%Nx' + "\n")
    file.write(str( int(data['Ny']) ) + '\t\t\t\t\t\t%Ny' + "\n")
    file.write(str( int(data['Nz']) ) + '\t\t\t\t\t\t%Nz' + "\n")
    file.write(str( int(data['Nd']) ) + '\t\t\t\t\t\t%Nd' + "\n")
    
    
    Lx = format( data['Lx'], '.16f')
    Lz = format( data['Lz'], '.16f')
    file.write(Lx + '\t\t%Lx' + "\n")
    file.write(Lz + '\t\t%Lz' + "\n")
    
    lx = data['Lx'] / (2. * math.pi) 
    lz = data['Lz'] / (2. * math.pi) 
    file.write(str( lx ) + '\t\t\t\t\t\t%lx=Lx/(2pi)' + "\n")
    file.write(str( lz ) + '\t\t\t\t\t\t%lz=Lz/(2pi)' + "\n")
    
    alpha = (2.* math.pi) / data['Lx']
    file.write(str( alpha ) + '\t\t\t\t\t\t%alpha=2pi/Lx' + "\n")
    
    file.close()
    
    return 0
    
    
def makeSolutionDirectory(data, directory):
    kx = data['kx']
    kz = data['kz']
    c = data['c']
    A = data['A']
    
    folderName = "/wave_packet_" + str(kx) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A)
    
    directory = directory + folderName
    
    
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    return directory
    
    
    
    
    
def perturbFlowField(data):
    if data['read']:
        u = data['flowField']['physical'][0, :, :, :]
        # perturb u vel
        u[:, 1, :] = 0.01
        u[:, 2, :] = 0.02
        u[:, 3, :] = 0.01
        data['flowField']['physical'] = u
        
        
    else:
        u = data['resolvent_flowField'][0, :, :, :]
        # perturb u vel
        u[:, 1, :] = 0.4
        u[:, 2, :] = 0.7
        u[:, 3, :] = 0.6
        u[:, 4, :] = 0.3
        data['resolvent_flowField'][0, :, :, :] = u
        data['U'] = u
        
    return data



def writeASCIIfile_general(data, directory):
    
    flowField = data['flowField']['physical']
    
    fileName = "/wave_packet_" + str(kx) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A) + ".asc"
    
    file = open(directory + fileName, "w")
    
    for nx in range(0, data['Nx']):
        for ny in range(0, data['Ny']):
            for nz in range(0, data['Nz']):
                for nd in range(0, data['Nd']):
                    tmp = flowField[nd, nx, ny, nz]
                    tmp = format(tmp, '.16f')
                    file.write(tmp + "\n")

    file.close()
    
    print('Boom, the ASCII is done.')
    
    return 0
    
    