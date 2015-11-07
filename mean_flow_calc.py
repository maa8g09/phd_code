import utils as ut
import read_construct as rc
import utils_plots as up
import matplotlib.pyplot as plt
import numpy as np
import subprocess as sp
import os
from colorama import Fore, Back, Style

from images2gif import writeGif
from PIL import Image
from time import time 
import sys

# In order to get the mean flow from the flow fields I have generated, I need to first:
#  - list the range I want
#  - convert to ASCII
#  - read ASCII and store the data for the ff in a dictionary:
#             flowFields[time] = U[nd, nx, ny, nz]
#  - all ffs have now been read in and stored.
#  - loop through dictionary
#       - sum all ffs
#  - divide answer by range length (number of ffs observed)
# This is the mean flow. 

# Now feed to get the fluctuations ASCII
#  - delete the mean ff from each ff in the massive dictionary
#  - output the fluctuations matrix in a format that Paraview can interpret

# Load into Paraview and viola.



def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)



main_directory = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_007_4modes_(-1.51191822883+0j)'
case_directory = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_007_4modes_(-1.51191822883+0j)/data-skew'
t_start = 0.0
t_end   = 999.0
steps = 1000
t_range = np.linspace(t_start, t_end, steps)



# Loop through this directory and make ascii files
os.chdir(case_directory)

directories = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
directories = sorted(directories)

files = [fi for fi in os.listdir(case_directory) if os.path.isfile(os.path.join(case_directory,fi))]
files = sorted(os.listdir(case_directory), key=os.path.getctime)


dict_inst = {}
dict_fluc = {}


cumulative = 0

final_data = {}

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
            
            if float(folderName) in t_range:
                print(folderName)
                
                ff_dir = case_directory + '/' + folderName
                
                if not os.path.exists(ff_dir):
                    os.mkdir(ff_dir)
                    ascOutput = folderName+'/'+folderName
                    sp.call(['field2ascii', '-p', k, ascOutput]) # channelflow command line arguments
                
                    
                    # Read solution
                    data = rc.main_read_construct(ff_dir, [0, 'all', 'all', 17])
                    final_data = data
                    cumulative = np.zeros((data['geometry']['physical']['Nd'],
                                           data['geometry']['physical']['Nx'],
                                           data['geometry']['physical']['Ny'],
                                           data['geometry']['physical']['Nz']))
                    mean = cumulative
                                           
                    dict_inst[folderName] = data['flowField']['physical']
                
                sp.call(['rm', '-rf', ff_dir])
                
                
print('Now i\'ve collated all the data from', t_start, 'till', t_end)





# Looping through the dictioanry
for key, value in dict_inst.items():
    cumulative += dict_inst[key]


magn = 1.0 / steps

# Mean
mean = cumulative * magn

ut.writeMeanASCIIfile(mean, case_directory)


post_image_dir = main_directory + '/images'
if not os.path.exists(post_image_dir):
    os.mkdir(post_image_dir)


mean_slice = rc.get_data_slice(mean, [0, 'all', 'all', 17], final_data['geometry']['physical'])
if os.path.exists(post_image_dir):
    fileName = 'mean'
    up.plot2D(mean_slice, post_image_dir, fileName, '')
        
        
        
        
        
        

for key, value in dict_inst.items():
    dict_fluc[key] = dict_inst[key] - mean
    


    
# Now we can plot the images
# First we need a dictionary of vel_slices.

dict_fluc_velSlices = {}
for key, value in dict_inst.items():
    ff = dict_fluc[key]
    
    fluc_slice = rc.get_data_slice(ff, [0, 0, 'all', 'all'], final_data['geometry']['physical'])
    
    if os.path.exists(post_image_dir):
        key = key.zfill(8)
        fileName = 'fluc_t_'+str(key)
        up.plot2D(fluc_slice, post_image_dir, fileName, str(key))
    
    print('\nkey done:', str(key), '\n')
    
        
        


