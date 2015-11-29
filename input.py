#!/usr/bin/env python


import numpy as np
import main_resolvent as m
import wavepackets as wp
import dns_routine as dns
from datetime import datetime
import utils as ut
import time
import os

pwd=os.getcwd()
print(pwd)

startTimeLarge = datetime.now()



# Case Booleans_____________________________________________________________________________________

generateInitialFF   = True

# Initial Condition
outputPhysASC       = True
outputSpecASC       = False
plotting            = False

# Convert to Binary for DNS
convert2ff          = False

# Write symmetry file
writeSymmFile       = False

# Run DNS
DNS                 = False
if DNS:
    convert2ff = False

# Check convergence
checkConvergence    = False
runSeriesDist       = False
if checkConvergence:
    runSeriesDist = False


"""
====================================================================================================
Initialize Variables 
====================================================================================================
"""


####################################################################################################
#### INITIALISE VARIABLES - Using existing flow field
####################################################################################################
mac = False
linux = False
approximate_soln  = linux

if mac:
    direct = '/Users/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ1'
elif linux:
    direct='/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25/wavepacket_017_4modes_(-0.0118649290443+0j)/data-skew/u200.00'
    ampl_weights = np.array([1.0])
    step = ''
else:
    direct = ''
    ampl_weights = np.logspace(1.0, -3.0, 5)
    #ampl_weights = np.logspace(-0.26, -0.47, 5)
    #ampl_weights = np.array([1e-1])
    date = time.strftime("%Y_%m_%d")
    step = '-DNS-' + str(date)

####################################################################################################
#### INITIALISE VARIABLES - Creating your own flow field
####################################################################################################
# Resolution in wall-normal direction________________________________
geom={}
geom['Nd']=3; geom['Nx']=36; geom['Ny']=37; geom['Nz']=36

# Reynolds number ___________________________________________________
# based on centreline/wall velocity for channel/coutte flow
Re = 3000

# Wave number triplet________________________________________________
# The mode combinations are taken from Sharma & McKeon (2013) paper
wavepacket = 'KB'

xlabel = wavepacket+'_x'
zlabel = wavepacket+'_z'
kx = wp.wavepackets[xlabel]
kz = wp.wavepackets[zlabel]

# Wave speed_________________________________________________________
c = 2.0 / 3.0

# Amplitudes_________________________________________________________
alabel = wavepacket+'_a'
a = wp.wavepackets[alabel]

# Plotting information for a slice___________________________________
fdary = [0 , 'all', 'all', 17] # XY 
fdary = [0, 0, 'all', 'all'] # YZ

# Output directory___________________________________________________
main_direc = '/home/arslan/Documents/work/channelflow-related/set01/Re' + str(Re) + '/' + wavepacket + '/ampls' + step

# Geometry___________________________________________________________
geom['Lx'] = 2.0*np.pi / kx[0]
geom['Lz'] = 2.0*np.pi / kz[0]

# Stationary nodes along each axis___________________________________
# X axis
geom['Mx_full'] = np.arange((-geom['Nx']/2.0), (geom['Nx']/2.0)+1) # Full range of Fourier modes
geom['Mx'] = np.arange(-1.0, 2.0) # 1 harmonic

# Z axis_____________________________________________________________
geom['Mz_full'] = np.arange(0, 0.5*(geom['Nz']) + 1.0) # full range of Fourier modes
geom['Mz'] = np.arange(2.0) # 1 harmonic

# X & Z axes_________________________________________________________
geom['x'] = np.linspace(0.0, geom['Lx'], geom['Nx'])
geom['z'] = np.linspace(-geom['Lz']/2.0, geom['Lz']/2.0, geom['Nz'])


geom['t'] = 0
# Total number of modes____________________________________________________
geom['m'] = geom['Ny'] - 2


####################################################################################################
#### ROUTINE: Generate Initial Flow Field & Conver to BINARY
####################################################################################################
if generateInitialFF:
    for i in range(0, len(ampl_weights)):
        ampl = a*ampl_weights[i]
        m.main(direct, fdary, geom, Re, kx, kz, c, ampl, i, main_direc, approximate_soln, outputPhysASC, outputSpecASC, plotting)


os.chdir(main_direc)


# Now we are going to loop through the directory and list directories
dirs = [d for d in os.listdir(main_direc) if os.path.isdir(d)]
dirs = sorted(dirs)
# now we loop thorugh all directories in main_direc and run the conver2ff method
for i in range(0, len(dirs)):
    case_direc = dirs[i]

    ####################################################################################################
    #### ROUTINE: Convert to Binary
    ####################################################################################################
    if convert2ff:
        ut.convert2ff(case_direc)

    ###################################################################################################
    #### ROUTINE: Write symmetry file
    ####################################################################################################
    symData = {}
    symData['fileName'] = 'sigma.asc'
    symData['writtenSymmsFile'] = writeSymmFile
    if writeSymmFile:
    # symmetry file: s sx sy sz ax az
        ut.writeSymmsFile(case_direc, symmsFileName, 1, ['1 1 1 1 0.5 0.0'])

    ####################################################################################################
    #### ROUTINE: DNS
    ####################################################################################################        
    if DNS:
        dns.main(case_direc, Re, symData)

print('\n    All donne!\n')
print('   ', datetime.now() - startTimeLarge, '\n')
