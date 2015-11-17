import os
import numpy as np
import subprocess as sp

dirc = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_007_4modes_(-1.51191822883+0j)/data-skew/'


#wavepacket_008_4modes_(-0.931112136502+0j)



# Name change all variables
# Take all the decimals out. 
# so a range of 500 - 999


t_start = 0
t_end = 999
steps = 1000
t_range = np.linspace(t_start, t_end, steps)




new_dir = dirc + 'ints'
if not os.path.exists(new_dir):
    os.mkdir(new_dir)        
os.chdir(new_dir)




files = [fi for fi in os.listdir(dirc) if os.path.isfile(os.path.join(dirc,fi))]

for k in files:
    if str(k) != 'couette.args' and str(k)[0] == 'u' and str(k) != 'umean.asc':
        oldFile = k
        k = k[1:-3]
        k = float(k)
        
        if k in t_range:
            k = int(k)
            
            newFile = 'u' + str(k) + '.ff'
            oldFile = '../' + oldFile
            command = ['cp', oldFile, newFile]
            
            sp.call(command)
            
    
                


print('\n   All done!   \n')

