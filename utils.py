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
import numpy as np
from scipy.interpolate import interp1d


verbose = False



def printStart():
    
    if verbose:
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
    if verbose:
        print('__________________________________________________________________\n')
    return


def printSectionTitle(str):
    # Print the section headers for each main section of the code output.
    if verbose:
        print(' **', str, '\n')
    return


def error(str):
    # Print the error and then exit from the program entirely
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('ERROR:', str)
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('\n')
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
    print('\n')
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('ERROR:', str)
    print('!!!!====!!!!====!!!!====!!!!====!!!!====')
    print('\n\n')
    sys.exit('')
    return
    
    
def message(str):
    if verbose:
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
    if verbose:
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
    
#    fileName = "/wave_packet_" + str(kx) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A) + ".asc"
    fileName = "/u0.asc"
    
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
    
#    fileName = "/wave_packet_" + str(kx) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A) + ".geom"
    fileName = "/u0.geom"
    
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
    gamma = (2.* math.pi) / data['Lz']
    file.write(str( gamma ) + '\t\t\t\t\t\t%gamma=2pi/Lz' + "\n")
    
    
    file.close()
    
    return 0
    
    
def makeSolutionDirectory(data, directory, n, re, kx, kz, c, amplitudes, i):
    i = i+1
    kxstr = data['kx']
    
    if kx.shape[0] > 1:
        i = str(i).zfill(3)
        folderName = "/wavepacket_" + str(i) + "_" + str(kx.shape[0]) + "modes_" + str(amplitudes[0])
        
    else:
        kxstr = data['kx']
        kz = data['kz']
        c = data['c']
        A = data['A']
        i = str(i).zfill(3)
        
        if kx[0][0] < 0:
            folderName = "/triplet_case_"+ str(i) +"_+-" + str(kxstr) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A)
        else:
            folderName = "/triplet_case_"+ str(i) +"_"   + str(kxstr) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A)
        
    directory = directory + folderName
    
    
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    fileName = '/details.txt'
    file = open(directory + fileName, "w")
    file.write('n (Ny - 2)'+ "\n")
    file.write(str(n)+ "\n"+ "\n")
    file.write('Re'+ "\n")
    file.write(str(re)+ "\n"+ "\n")
    file.write('kx'+ "\n")
    file.write(str(kx)+ "\n"+ "\n")
    file.write('kz'+ "\n")
    file.write(str(kz)+ "\n"+ "\n")
    file.write('c'+ "\n")
    file.write(str(c)+ "\n"+ "\n")
    file.write('amplitudes'+ "\n")
    file.write(str(amplitudes)+ "\n"+ "\n")
    
    file.close()
    
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
    














def writeSpectralASCIIfile(data, directory):
    
    flowField = data['spectral_ff']
    fileName = "/u0_spec_rank-"+str(data['Rank'])+".asc"
    file = open(directory + fileName, "w")
    
    for mx in range(0, data['Mx']):
        for ny in range(0, data['Ny']):
            for mz in range(0, data['Mz']):
                for nd in range(0, data['Nd']):
                    tmpReal = flowField[nd, mx, ny, mz].real
                    tmpImag = flowField[nd, mx, ny, mz].imag
#                    tmpReal = format(tmpReal, '.8f')
#                    tmpImag = format(tmpImag, '.8f')
                    
                    output = '(' + str(tmpReal) + ', ' + str(tmpImag) + ')'
#                    if tmpImag >= 0.0:                        
#                        output = str(tmpReal) +'+'+ str(tmpImag)+'j'
#                    elif tmpImag < 0.0:
#                        output = str(tmpReal) + str(tmpImag)+'j'
                    
                    file.write(output + "\n")
    
    file.close()
    
    print('Boom, the SPECTRAL ASCII is done.')
    
    return 0
    
    
    
    
def writeSpectralGEOMfile(data, directory):
    
    kx = data['kx']
    kz = data['kz']
    c = data['c']
    A = data['A']
    
#    fileName = "/wave_packet_" + str(kx) + "_+-" + str(kz) + "_" + str(c) + "_" + str(A) + ".geom"
    fileName = "/u0_spec_rank-"+str(data['Rank'])+".geom"
    
    file = open(directory + fileName, "w")
    
    
#    file.write(str( int(len(data['Mx'])) ) + '\t\t\t\t\t\t%Mx' + "\n")
#    file.write(str( int(len(data['Mz'])) ) + '\t\t\t\t\t\t%Mz' + "\n")
    
    
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
    file.write(str( lx ) + '\t\t%lx=Lx/(2pi)' + "\n")
    file.write(str( lz ) + '\t\t%lz=Lz/(2pi)' + "\n")
    
    alpha = str(data['fund_alpha'])
    file.write(str( alpha ) + '\t\t%alpha=2pi/Lx' + "\n")
    gamma = str(data['fund_beta'])
    file.write(str( gamma ) + '\t\t%gamma=2pi/Lz' + "\n")
    
    
    file.close()
    
    return 0



def makeOutputDictionary(generated_ff, geom, y_cheb, uniform, string_kx, string_kz, string_c, string_A):
    
    outputDic = {}
    
    U = np.zeros((geom['Nd'], geom['Nx'], 3*geom['m'], geom['Nz']))
    
    U_u = generated_ff.real[:,           0:geom['m']  , :]
    U_v = generated_ff.real[:,   geom['m']:2*geom['m'], :]
    U_w = generated_ff.real[:, 2*geom['m']:3*geom['m'], :]
    
    for i in range(0, geom['Nd']):
        for nx in range(0, geom['Nx']):
            for ny in range(0, geom['m']):
                for nz in range(0, geom['Nz']):
                    if i == 0: # u direction
                        U[i, nx, ny, nz] = U_u[nx, ny, nz]
                    elif i == 1: # v direction
                        U[i, nx, ny, nz] = U_v[nx, ny, nz]
                    elif i == 2: # w direction
                        U[i, nx, ny, nz] = U_w[nx, ny, nz]

    
#    L2Norm = np.linalg.norm(U)
#    print(np.allclose(L2Norm, np.sqrt(np.sum(np.square(U[:,:,:,:])))))
#    magn = 10.0
#    U *= magn / L2Norm
    
    
    # Interpolation to go from y_cheb toy_uniform
    Ny = geom['m']
    y_uniform = np.linspace(1.0, -1.0, Ny*1.0)
    y_cheb = np.asarray(y_cheb)
    y_cheb = np.squeeze(y_cheb)
    
    
    U_u_uniform = np.zeros((geom['Nx'], geom['m'], geom['Nz']))
    U_v_uniform = np.zeros((geom['Nx'], geom['m'], geom['Nz']))
    U_w_uniform = np.zeros((geom['Nx'], geom['m'], geom['Nz']))
    
    for nx in range(0, geom['Nx']):
        for nz in range(0, geom['Nz']):
            uprofile = U_u[nx, :, nz] # 1-d vector
            # fill value is the no-slip boundary condition
            fu = interp1d(y_cheb, uprofile, bounds_error=False, fill_value=0.0, kind='cubic') 
            fu = fu(y_uniform)
            U_u_uniform[nx, :, nz] = fu
            
            
            vprofile=U_v[nx, :, nz] # 1-d vector
            fv = interp1d(y_cheb, vprofile, bounds_error=False, fill_value=0.0, kind='cubic') 
            fv = fv(y_uniform)
            U_v_uniform[nx, :, nz] = fv
            
            wprofile=U_w[nx, :, nz] # 1-d vector
            fw = interp1d(y_cheb, wprofile, bounds_error=False, fill_value=0.0, kind='cubic') 
            fw = fw(y_uniform)
            U_w_uniform[nx, :, nz] = fw


#    plt.plot(y_cheb, uprofile, 'r-', y_uniform, fu, 'g--')
#    plt.legend(['data', 'cubic'], loc='best')
#    plt.grid(True)
#    plt.show()
    

    outputDic['resolvent_flowField'] = U
    

    if uniform:
        outputDic['U'] = U_u_uniform
        outputDic['V'] = U_v_uniform
        outputDic['W'] = U_w_uniform
        
        outputDic['Y'] = y_uniform
    
    else:
        outputDic['U'] = U_u
        outputDic['V'] = U_v
        outputDic['W'] = U_w
        
        outputDic['Y'] = y_cheb
        
    
    outputDic['X'] = geom['x']
    outputDic['Z'] = geom['z']

    outputDic['Nx'] = geom['Nx']
    outputDic['Ny'] = geom['m']
    outputDic['Nz'] = geom['Nz']
    outputDic['Nd'] = geom['Nd']

    outputDic['Lx'] = geom['Lx']
    outputDic['Lz'] = geom['Lz']
    
    outputDic['kx'] = string_kx
    outputDic['kz'] = string_kz
    outputDic['c'] = string_c
    outputDic['A'] = string_A
    
    
    return outputDic
    
