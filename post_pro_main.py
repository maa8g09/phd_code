import utils as ut
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess as sp
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











main_direc = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/postprotests/'
dns_output = 'data_skew.txt'


# Changing directories and looping through listing all subdirectories.
os.chdir(main_direc)


directories = [d for d in os.listdir(os.getcwd()) if os.path.isdir(d)]
directories = sorted(directories)

legnd = []

for i in range(0, len(directories)):
    print(directories[i])    
    os.chdir(main_direc + directories[i])   
    pwd = os.getcwd() 
    print(pwd)
    files = [fi for fi in os.listdir(pwd) if os.path.isfile(os.path.join(pwd,fi))]
    
    tmp_dir = str(directories[i])
    tmp = '_0'
    tmp_i = tmp_dir.find(tmp)
    tmp_dir = tmp_dir[tmp_i+1: tmp_i+4]
    tmp_dir = str(int(tmp_dir)+1)
    tmp_dir = '$\chi_{' + tmp_dir+'}$'
    legnd.append(tmp_dir)
    
    for j in files:
        if str(j) == dns_output:
            ut.printSectionHeader()
            ut.printSectionTitle('Reading the DNS output file')
                    
            
            f = ut.openFile(pwd + '/' + dns_output)
            
            ut.message('Constructing dictionary')
            var = {}
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
            
            var['t'] = []
            var['L2Norm(u)'] = []
            var['dissip(u+U)'] = []
            var['CFL'] = []
            
            for k, line in enumerate(f):
                values = line.split()
                
                if len(values) == 0:
                    print('Empty line')
                    
                else:
                    
                    if values[0] == 't':
                        if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                            var['t'].append(float(values[2]))
            
                    elif values[0] == 'L2Norm(u)':
                        if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                            var['L2Norm(u)'].append(float(values[2]))
                
                    elif values[0] == 'dissip(u+U)':
                        if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                            var['dissip(u+U)'].append(float(values[2]))
            
                    elif values[0] == 'CFL':
                        if values[2] != '-nan' or values[2] != 'nan' or values[2] != 'done!':
                            var['CFL'].append(float(values[2]))
            
            ut.message('Closing file')
            f.close()
            
            var['t']=np.asarray(var['t'])
            var['L2Norm(u)']=np.asarray(var['L2Norm(u)'])
            var['dissip(u+U)']=np.asarray(var['dissip(u+U)'])
            var['CFL']=np.asarray(var['CFL'])
            
            
            
            # Now plot everything contained in the dictionary:
            for key, value in var.items():
                
                if key != 't':
                    color = tuple(np.random.rand(4))
                    plt.plot(var['t'], var[key], c=color)
                    
                    plt.grid(True)
                    plt.xlabel('Time units, t')
                    plt.ylabel(key)
                    plt.legend(legnd, loc='best')
                    
                    plt.savefig(main_direc + 'convergence_' + key + '.png')   # save the figure to file
        
            
    
    # Now we are gonna go through and make subdirectories for each time step.
    subdirecs = [d for d in os.listdir(pwd) if os.path.isdir(d)]
    for j in range(0, len(subdirecs)):
        os.chdir(pwd + '/' + subdirecs[j])   
        subpwd = os.getcwd()
        print(subpwd)
        files = [fi for fi in os.listdir(subpwd) if os.path.isfile(os.path.join(subpwd,fi))]
    
        for k in files:
            k = str(k)
            folderName = k[:-3]
            ff_dir = subpwd + '/' + folderName
            #if directory doesn't exist:
            if not os.path.exists(ff_dir):
                os.mkdir(ff_dir)
                ascOutput = folderName+'/'+folderName
                sp.call(['field2ascii', '-p', k, ascOutput]) # channelflow command line arguments
            
                
                subsubdir = os.chdir(ff_dir)
                
                plotSlice(inputDir, outputDir)
                
                # Now take the asc files and plot slice. 
                
                # NOTE: change field2ascii to only output the physical information, there's no need for the spectral 
                # information.
                
                # isolate the looping through directories and executing commands.
            






















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
