import numpy as np
import main as m
import wavepackets as wp
from datetime import datetime
import time
import os
pwd=os.getcwd()
print(pwd)
startTimeLarge = datetime.now()

approximate_soln  = False
outputPhysicalASC = False
outputSpectralASC = False
plotting          = True

####################################################################################################
#### INPUT PARAMETERS - Using existing flow field
####################################################################################################
mac = False
linux = True

if mac:
    direct = '/Users/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ1'
elif linux:
    direct='/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25/wavepacket_017_4modes_(-0.0118649290443+0j)/data-skew/u200.00'
    ampl_weights = np.array([1.0])
    step = ''
else:
    direct = ''
    ampl_weights = np.logspace(1.0, -3.0, 20)
#    ampl_weights = np.array([1e-1])
    date = time.strftime("%Y_%m_%d")
    step = '-DNS-' + str(date)
    
####################################################################################################
#### INPUT PARAMETERS - Creating your own flow field
####################################################################################################
# Resolution in wall-normal direction________________________________
geom={}
geom['Nd']=3; geom['Nx']=36; geom['Ny']=37; geom['Nz']=36

# Reynolds number ___________________________________________________
# based on centreline/wall velocity for channel/coutte flow
Re = 1200

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
#fdary = [0, 0, 'all', 'all'] # YZ

# Output directory___________________________________________________
d = '/home/arslan/Documents/work/channelflow-related/set01/Re' + str(Re) + '/' + wavepacket + '/ampls' + step

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

# run-time___________________________________________________________
geom['t'] = 0            
# number of modes____________________________________________________
geom['m'] = geom['Ny'] - 2  
    
####################################################################################################
#### MAIN ROUTINE
####################################################################################################
for i in range(0, len(ampl_weights)):
    ampl = a*ampl_weights[i]
    m.main(direct, fdary, geom, Re, kx, kz, c, ampl, i, d, approximate_soln, outputPhysicalASC, outputSpectralASC, plotting)

print('   ', datetime.now() - startTimeLarge, '\n')
