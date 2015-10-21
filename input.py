import numpy as np
import main as m
import wavepackets as wp
from datetime import datetime
startTimeLarge = datetime.now()

####################################################################################################
#### INPUT PARAMETERS - Using existing flow field
####################################################################################################
mac = False
linux = False

if mac:
    direct = '/Users/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ1'
elif linux:
    direct='/home/arslan/Documents/work/channelflow-related/database_solns/W03/equilibria/EQ1/findsoln_stuff/switch-OFF/rank-090'
    ampl_weights = np.array([1.0])
else:
    direct = ''
    ampl_weights = np.logspace(1.0, -4.0, 10)
#    ampl_weights = np.array([1.0])
    step = '01'
    
####################################################################################################
#### INPUT PARAMETERS - Creating your own flow field
####################################################################################################
# Resolution in wall-normal direction
geom={}
geom['Nd']=3; geom['Nx']=40; geom['Ny']=52; geom['Nz']=40

# Reynolds number based on centreline/wall velocity for channel/coutte flow
Re = 1200

# Wave number triplet
# The mode combinations are taken from Sharma & McKeon (2013) paper
wavepacket = 'KB'

xlabel = wavepacket+'_x'
zlabel = wavepacket+'_z'
kx = wp.wavepackets[xlabel]
kz = wp.wavepackets[zlabel]

# Wave speed
c = 2.0 / 3.0

# Amplitudes
alabel = wavepacket+'_a'
a = wp.wavepackets[alabel]

# Plotting information for a slice
fdary = [0, 'all', 'all', 0] # XY 
fdary = [0, 0, 'all', 'all'] # YZ

# Output directory
d = '/home/arslan/Documents/work/channelflow-related/set01/Re' + str(Re) + '/' + wavepacket + '/ampls' + step

# Geometry
geom['Lx'] = 2.0*np.pi / kx[0]
geom['Lz'] = 2.0*np.pi / kz[0]

# Stationary nodes along each axis
# X axis
geom['Mx_full'] = np.arange((-geom['Nx']/2.0), (geom['Nx']/2.0)+1) # Full range of Fourier modes
geom['Mx'] = np.arange(-1.0, 2.0) # 1 harmonic

# Z axis
geom['Mz_full'] = np.arange(0, 0.5*(geom['Nz']) + 1.0) # full range of Fourier modes
geom['Mz'] = np.arange(2.0) # 1 harmonic



# X & Z axes
geom['x'] = np.linspace(0.0, geom['Lx'], geom['Nx'])
geom['z'] = np.linspace(-geom['Lz']/2.0, geom['Lz']/2.0, geom['Nz'])


geom['t'] = 0               # run-time
geom['m'] = geom['Ny'] - 2  # number of modes
    
    
for i in range(0, len(ampl_weights)):
    ampl = a*ampl_weights[i]
    m.main(direct, fdary, geom, Re, kx, kz, c, ampl, i, d)

print('   ', datetime.now() - startTimeLarge, '\n')
