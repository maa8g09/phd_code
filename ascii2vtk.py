import os
import numpy as np
import subprocess as sp
from colorama import Fore, Back, Style
import read_construct as rc
from pyevtk.hl import gridToVTK



data_dir = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_007_4modes_(-1.51191822883+0j)/data-skew/'


# Name change all variables
# Take all the decimals out. 
# so a range of 500 - 999


t_start = 0.0
t_end = 10.0
steps = 101
t_range = np.linspace(t_start, t_end, steps)



def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


os.chdir(data_dir)

files = [fi for fi in os.listdir(data_dir) if os.path.isfile(os.path.join(data_dir,fi))]

for k in files:
    if str(k) != 'couette.args' and str(k)[0] == 'u' and str(k) != 'umean.asc':
        
        k = str(k)
        folderName = k[1:-3]
        
        bool_start_pt = isclose(t_start, float(folderName))
        bool_end_pt   = isclose(t_end,   float(folderName))
        bool_in_range = False
        
        if float(folderName) > t_start and float(folderName) < t_end:
            bool_in_range = True
        
        
        if bool_in_range or bool_start_pt or bool_end_pt:
            iteration_time = folderName
            folderName = folderName.zfill(8)
            folderName = folderName.replace(".", "")
            
            
            
            print('Made\n\n',Fore.BLACK + Back.CYAN + folderName + Style.RESET_ALL)
            
            asc_dir = data_dir + folderName
            #if directory doesn't exist:
            if os.path.exists(asc_dir):
                # Delete the folder witht eh ASCII file and geom file.
                sp.call(['rm', '-rf', asc_dir])
                
                
            if not os.path.exists(asc_dir): # per flow field
                os.mkdir(asc_dir)
                ascOutput = folderName+'/'+folderName
                sp.call(['field2ascii', '-p', k, ascOutput])
    
    
                data = rc.main_read_construct(asc_dir, [0, 'all', 'all', 17])
                
                
                
                #### CONVERSION
                # Now we convert the data from numpy to vtk.
                u_comp = data['flowField']['physical'][0, :, :, :]
                x = np.linspace(0, data['geometry']['physical']['Lx'], data['geometry']['physical']['Nx'])
                y = 0
                
                gridToVTK("./julia", x, y, z, cellData = {'julia': u_comp})

                                
            # Delete the folder witht eh ASCII file and geom file.
            sp.call(['rm', '-rf', asc_dir])

            print('Deleted\n\n', folderName)


print('\n   All done!   \n')

