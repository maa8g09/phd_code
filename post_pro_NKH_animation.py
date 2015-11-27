#!/usr/bin/env python
import os
import subprocess as sp
import read_construct as rc
import utils as ut
import numpy as np
import math
####################################################################################################
#### Post-processing the newton search cases...                                                 ####
####################################################################################################

# First we have to go through and find all newton steps and put them into their own folders...

main_direc = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_007_4modes_(-1.51191822883+0j)'

case_direc = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_007_4modes_(-1.51191822883+0j)/nkh-vdt'


# Change to the case directory
os.chdir(case_direc)

# Loop through all files and find newton step flow fields
whatWereLookingForIs = 'unewt'
correctFileType = '.ff'

files = [fi for fi in os.listdir(case_direc) if os.path.isfile(os.path.join(case_direc,fi))]

files = sorted(files)


init_cond = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_007_4modes_(-1.51191822883+0j)/initialCondition/'


initialFF = rc.main_read_construct(init_cond, [0,0,'all', 'all'], False)
meanFile = ut.openFile(main_direc + '/data-skew/umean.asc')
meanVec = np.zeros((initialFF['geometry']['physical']['Nd'],
                    initialFF['geometry']['physical']['Nx'],
                    initialFF['geometry']['physical']['Ny'],
                    initialFF['geometry']['physical']['Nz']))

for u, line in enumerate(meanFile):
    values = line.split()
    nx = int(values[0])
    ny = int(values[1])
    nz = int(values[2])
    nd = int(values[3])
    vel = float(values[4])
    if vel <= 1e-8:
        vel = 0.0
    meanVec[nd, nx, ny, nz] = vel
    
ut.message('Closing the physical ASCII file')
meanFile.close()





data_each_NS = {}





for k in files:
    if whatWereLookingForIs in str(k) and correctFileType in str(k):
        folderName = k[0:-3]
        newtonStep = folderName[5:]
        
        folder_direc = case_direc + '/' + folderName

        if os.path.exists(folder_direc):
            sp.call(['rm', '-rf', folder_direc])
        
        if not os.path.exists(folder_direc):
            os.mkdir(folder_direc)
            ascOutput = folderName+'/'+folderName
            sp.call(['field2ascii', '-p', k, ascOutput]) # channelflow command line arguments
            
            # Now conver the ASC to DAT
            #### READ ASC
            data = rc.main_read_construct(folder_direc, [0, 'all', 'all', 17], True)
            
            
            # Calculate the fluctuations
            data['flowField']['physical'] = data['flowField']['physical'] - meanVec
            data_each_NS[newtonStep] = data['flowField']['physical']
            
            #### CONVERT 2 DAT
#            ut.write_DAT_file_from_ff(data, case_direc + '/', newtonStep, folderName)
            
print('Written separate DAT files')
print('Now writing a large one for animation')

sortedKeys = sorted(data_each_NS)
sortedKeys = map(int, sortedKeys)
sortedKeys = sorted(sortedKeys)


# Write a '.dat' file for the geometry info for Paraview
fileName = 'uAllSteps.dat'
file = open(case_direc + '/' + fileName, "w")


Nx = data['geometry']['physical']['Nx']
Ny = data['geometry']['physical']['Ny']
Nz = data['geometry']['physical']['Nz']

Lx = data['geometry']['physical']['Lx']
Lz = data['geometry']['physical']['Lz']

x = np.linspace(0.0, Lx, Nx)
z = np.linspace(-Lz/2.0, Lz/2.0, Nz)
y = np.linspace(-1.0,1.0,Ny)
for ny in range(0, Ny):
    y[ny] = math.cos(ny*math.pi/(Ny-1))


title   = 'TITLE= "Flow field through all the Newton-Krylov-hookstep steps"\n' 
columns = 'VARIABLES = "X", "Y", "Z", "U", "V", "W", "t"\n'
zones   = 'ZONE I=' + str( int(Nx)) + ', J=' + str( int(Ny) ) + ', K=' + str( int(Nz) ) + ', F=POINT\n'


file.write(title)
file.write(columns)
file.write(zones)


for k in sortedKeys:
    flowField = data_each_NS[str(k)]
    for nx in range(0, Nx):
        for ny in range(0, Ny):
            for nz in range(0, Nz):
                string = format(x[nx], '.8f') + ' '
                string += format(y[ny], '.8f') + ' '
                string += format(z[nz], '.8f') + ' '
                string += format(flowField[0, nx, ny, nz], '.16f') + ' '
                string += format(flowField[1, nx, ny, nz], '.16f') + ' '
                string += format(flowField[2, nx, ny, nz], '.16f') + ' '
                string += format(k, '.2f')
                file.write(string+'\n')
                
file.close()


print('\nDone!\n')










