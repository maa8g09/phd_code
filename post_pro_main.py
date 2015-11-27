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




convPlots = True
makingMovie = False
outputPlots = False
convertDAT = False


calculatedMean = False

'''

INPUTS:
    String main_directory:      directory where all cases are kept.
    String  dns_directory:      DNS-data folder name
    String     dns_output:      DNS-output file name

Change the directory of the execution to main_directory.

List all folders and loop through
change directory to each folder and in each folder:

    read the DNS-output .txt file

    loop through each line:
        store t, norm, ubulk, dissipation etc.
        
    plot any variable against t. 
    save the plot in the current working directory...
    
    go into the data directory and depending on the end time, make intervalled folders
    i.e. for 200time units, make 50 folders, separated by 4 flowfields.
    In each folder:
        execute the field2ascii command
        run main.py to plot the flow field and save the image in post-pro/xy-plane/flowfieldName.png



Make an automated input file, and keep the manual input separate. This way, the post-pro remains entirely blackboxed.

Just set main.py to output nothing and only plot the flow field. 





Questions:

What slice do you want the movies for? it's not steady state, so you either go through 
and output all images at each XY plane at each time step?


I think you should specify a slice you wanna look at (a la 4D array in main)

then look at that slice for however many time units you want to see. 

Put each 4D slice in its own folder.

'''
from datetime import datetime

pwd=os.getcwd()
print(pwd)
startTimeLarge = datetime.now()




def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)







main_direc = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/'
dns_output = 'data_skew.txt'
t_start = 0.0
t_end   = 500.0
steps = 501
t_range = np.linspace(t_start, t_end, steps)




# Changing directories and looping through listing all subdirectories.
os.chdir(main_direc)


directories = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
directories = sorted(directories)

legnd = []
var = {}


data_per_case = {}






















for i in range(0, len(directories)):
    startTimePerDir = datetime.now()
    
    
    if directories[i] != 'further':
        
        print('\n\n',directories[i], '\n')
        os.chdir(main_direc + '/' + directories[i])   
        pwd = os.getcwd() 
        print(pwd)
        
        
#        initialFF = rc.main_read_construct(pwd, [0,0,'all', 'all'], False)
#        
#        data_per_case = initialFF
#        
#        meanFile = ut.openFile(pwd + '/data-skew/umean.asc')
#        meanVec = np.zeros((initialFF['geometry']['physical']['Nd'],
#                            initialFF['geometry']['physical']['Nx'],
#                            initialFF['geometry']['physical']['Ny'],
#                            initialFF['geometry']['physical']['Nz']))
#        
#        for u, line in enumerate(meanFile):
#            values = line.split()
#            nx = int(values[0])
#            ny = int(values[1])
#            nz = int(values[2])
#            nd = int(values[3])
#            vel = float(values[4])
#            if vel <= 1e-8:
#                vel = 0.0
#            meanVec[nd, nx, ny, nz] = vel
#            
#        ut.message('Closing the physical ASCII file')
#        meanFile.close()
        
        
        
        
        files = [fi for fi in os.listdir(pwd) if os.path.isfile(os.path.join(pwd,fi))]
        
        tmp_dir = str(directories[i])
        tmp = '_0'
        tmp_i = tmp_dir.find(tmp)
        tmp_dir = tmp_dir[tmp_i+1: tmp_i+4]
        tmp_dir = str(int(tmp_dir))
        legnd.append('$\chi_{' + tmp_dir+'}$')
        
        
        if convPlots:
        
            for j in files:
                if str(j) == dns_output:
                    ut.printSectionHeader()
                    ut.printSectionTitle('Reading the DNS output file')
                            
                    
                    f = ut.openFile(pwd + '/' + dns_output)
                    
                    dict_per_case = {}
                    
                    
                    ut.message('Constructing dictionary')
                    
                    '''
                    Sample output, you can decide which variables you want ot plot from here.
                    
                               t == 42.8
                       L2Norm(u) == 0.704647
                    chebyNorm(u) == 1.0705
                     dissip(u+U) == 2.59421
                      input(u+U) == 1.00164
                      divNorm(u) == 0.0407906
                           Ubulk == 1.33333
                           ubulk == 0.666667
                            dPdx == -0.0040663
                      Ubulk*h/nu == 1600
                     Uparab*h/nu == 2400
                             CFL == 0.0220104
                    '''
                    
                    dict_per_case['t'] = []
                    dict_per_case['L2Norm(u)'] = []
                    dict_per_case['dissip(u+U)'] = []
                    dict_per_case['CFL'] = []
                    
                    for k, line in enumerate(f):
                        values = line.split()
                        
                        if len(values) == 0:
                            print('Empty line')
                            print('')
                            
                        else:
                            
                            if values[0] == 't':
                                if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                                    dict_per_case['t'].append(float(values[2]))
                    
                            elif values[0] == 'L2Norm(u)':
                                if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                                    dict_per_case['L2Norm(u)'].append(float(values[2]))
                        
                            elif values[0] == 'dissip(u+U)':
                                if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                                    dict_per_case['dissip(u+U)'].append(float(values[2]))
                    
                            elif values[0] == 'CFL':
                                if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                                    dict_per_case['CFL'].append(float(values[2]))
                    
                    ut.message('Closing file')
                    f.close()
                    
                    dict_per_case['t']=np.asarray(dict_per_case['t'])
                    dict_per_case['L2Norm(u)']=np.asarray(dict_per_case['L2Norm(u)'])
                    dict_per_case['dissip(u+U)']=np.asarray(dict_per_case['dissip(u+U)'])
                    dict_per_case['CFL']=np.asarray(dict_per_case['CFL'])
                    var[tmp_dir] = dict_per_case
                    
           
        
        if makingMovie:
            
            
            # Now we are gonna go through and make subdirectories for each time step.
            subdirecs = [d for d in os.listdir(pwd) if os.path.isdir(d)]
            subdirecs = sorted(subdirecs)
            
            for j in range(0, len(subdirecs)):
                
                folder_string='data'
                
                if folder_string in str(subdirecs[j]):
                    
                    os.chdir(pwd + '/' + subdirecs[j])   
                    subpwd = os.getcwd()
        #            print(subpwd)
                    files = [fi for fi in os.listdir(subpwd) if os.path.isfile(os.path.join(subpwd,fi))]
                    files = sorted(os.listdir(subpwd), key=os.path.getctime)
            
                    
                    post_image_dir = pwd + '/images'
                    if not os.path.exists(post_image_dir):
                        os.mkdir(post_image_dir)
                        
                    
                    
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
                                
                                    iteration_time = folderName
                                    folderName = folderName.zfill(8)
                                    folderName = folderName.replace(".", "")
                                    folderName = 'img'+folderName[:-2]
                                    
                                    
                                    print('Made\n',Fore.BLACK + Back.CYAN + folderName + Style.RESET_ALL)
                                    
                                    ff_dir = subpwd + '/' + folderName
                                    #if director does exist:
                                    if os.path.exists(ff_dir):
                                        sp.call(['rm', '-rf', ff_dir])
                                    
                                    #if directory doesn't exist:
                                    if not os.path.exists(ff_dir):
                                        os.mkdir(ff_dir)
                                        ascOutput = folderName+'/'+folderName
                                        
                                        #### CONVERT 2 ASC
                                        sp.call(['field2ascii', '-p', k, ascOutput]) # channelflow command line arguments
                                        
                                        
                                        #### READ ASC
                                        data = rc.main_read_construct(ff_dir, [0, 'all', 'all', 17], True)
                                        data_per_case['flowField'] = data['flowField']
                                        
                                        #### CALCULATE FLUCTUATIONS
                                        data_per_case['flowField']['physical'] = data['flowField']['physical'] - meanVec
                                        
                                        
                                        #### CONVERT 2 DAT
                                        if convertDAT:
                                            ut.write_DAT_file_from_ff(data_per_case, subpwd, iteration_time, folderName[3:]+'_flucs')
                                        
                                        # Plotting
                                        if outputPlots:
                                            if data['flowField']['is_physical'] == True:
                                                if os.path.exists(post_image_dir):
                                                    up.plot2D(data['velslice'], post_image_dir, folderName, iteration_time)
                                        
                                    
                                    # Delete the folder witht eh ASCII file and geom file.
                                    sp.call(['rm', '-rf', ff_dir])
                        
                                    print('Deleted\n\n', folderName)
                                
                    text04='Done with the flow fields'
                    print('\n',Fore.WHITE + Back.MAGENTA + text04 + Style.RESET_ALL,'\n')
        
                else: # We're now in images directory
                    print(str(subdirecs[j]))
    

    
    
    
    text01 = str(datetime.now() - startTimePerDir)
    print('\n',Fore.BLACK + Back.GREEN + text01 + Style.RESET_ALL,'\n')











if convPlots:
    # Now plot everything contained in the dictionary
    t=[]
    L2norm=[]
    dissip=[]
    cfl=[]
         
    for key, value in var.items():
        # Loop through each of the keys (cases)
        t.append(var[key]['t'])
        L2norm.append(var[key]['L2Norm(u)'])
        dissip.append(var[key]['dissip(u+U)'])
        cfl.append(var[key]['CFL'])
    
    num_plots = len(t)
    
    colormap = plt.cm.gist_ncar
    plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, num_plots)])
    #plt.figure()
    
    for i in range(len(t)):
        plt.plot(t[i], L2norm[i])
        plt.xlabel('Time units, t')
        plt.ylabel('L2Norm(u)')
        plt.grid(True)
        plt.xlim([0,1000])
        plt.ylim([0.68,0.74])
    
    plt.legend(legnd, loc='best') 
    plt.savefig(main_direc + '/convergence_L2norm.png')   # save the figure to file   
    
#    
#    for i in range(len(t)):
#        plt.plot(t[i], dissip[i])
#        plt.ylabel('dissip(u+U)')
#        plt.grid(True)
#        plt.xlim([0,520])
#        plt.ylim([0.0,20])
#    
#    plt.legend(legnd, loc='best')
#    plt.savefig(main_direc + '/convergence_dissip.png')   # save the figure to file
    
    #plt.legend(legnd, ncol=4, loc='upper center', 
    #           bbox_to_anchor=[0.5, 1.1], 
    #           columnspacing=1.0, labelspacing=0.0,
    #           handletextpad=0.0, handlelength=1.5,
    #           fancybox=True, shadow=True)
    #           
    #    
    #
    #
    #    if key != 't':
    #        color = tuple(np.random.rand(4))
    #        fig = plt.figure()
    #        plt.plot(var['t'], var[key], c=color)
    #        
    #        plt.grid(True)
    #        plt.xlabel('Time units, t')
    #        plt.ylabel(key)
    #        plt.legend(legnd, loc='best')
    #        plt.show()
    ##        
    ##        plt.savefig(main_direc + '/convergence_' + key + '.png')   # save the figure to file
    ##        
    ##        plt.close(fig)

    





text02 = str(datetime.now() - startTimeLarge)
print('   ', Fore.YELLOW + Back.BLUE + text02 + Style.RESET_ALL, '\n')










#
## A CRUDE VERSION (NEEDS TO BE TIDIED):
#
## Convergence plot
#main_direc = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25/'
#amp = 'wavepacket_018_4modes_(-0.00730699532635+0j)'
#main_direc = main_direc + amp
#dns_output = 'data_skew.txt'
#
#
#
#ut.printSectionHeader()
#ut.printSectionTitle('Reading the DNS output file')
#        
#
#f = ut.openFile(main_direc + '/' + dns_output)
#
#ut.message('Constructing dictionary')
#var = {}
#'''
#           t == 42.8
#   L2Norm(u) == 0.704647
#chebyNorm(u) == 1.0705
# dissip(u+U) == 2.59421
#  input(u+U) == 1.00164
#  divNorm(u) == 0.0407906
#       Ubulk == 1.33333
#       ubulk == 0.666667
#        dPdx == -0.0040663
#  Ubulk*h/nu == 1600
# Uparab*h/nu == 2400
#         CFL == 0.0220104
#'''
#var['t'] = []
#var['L2Norm(u)'] = []
#var['dissip(u+U)'] = []
#
#for i, line in enumerate(f):
#    values = line.split()
#    
#    if len(values) == 0:
#        print('Empty line')
#        
#    else:
#        
#        if values[0] == 't':
#            if values[2] != '-nan' or values[2] != 'nan':
#                var['t'].append(float(values[2]))
#
#        elif values[0] == 'L2Norm(u)':
#            if values[2] != '-nan' or values[2] != 'nan':
#                var['L2Norm(u)'].append(float(values[2]))
#    
#        elif values[0] == 'dissip(u+U)':
#            if values[2] != '-nan' or values[2] != 'nan':
#                var['dissip(u+U)'].append(float(values[2]))
#
#
#ut.message('Closing file')
#f.close()
#        
#
#var['t']=np.asarray(var['t'])
#var['L2Norm(u)']=np.asarray(var['L2Norm(u)'])
#var['dissip(u+U)']=np.asarray(var['dissip(u+U)'])
#
#
#plt.plot(var['t'], var['L2Norm(u)'], 'r-')
#plt.grid(True)
#
#plt.xlabel('Time units, t')
#plt.ylabel('L2Norm(u)')
#
#plt.show()
