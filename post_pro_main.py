import utils as ut
import numpy as np
import matplotlib.pyplot as plt
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
'''


# A CRUDE VERSION (NEEDS TO BE TIDIED):

# Convergence plot
main_direc = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25/'
amp = 'wavepacket_018_4modes_(-0.00730699532635+0j)'
main_direc = main_direc + amp
dns_output = 'data_skew.txt'



ut.printSectionHeader()
ut.printSectionTitle('Reading the DNS output file')
        

f = ut.openFile(main_direc + '/' + dns_output)

ut.message('Constructing dictionary')
var = {}
'''
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

for i, line in enumerate(f):
    values = line.split()
    
    if len(values) == 0:
        print('Empty line')
        
    else:
        
        if values[0] == 't':
            if values[2] != '-nan' or values[2] != 'nan':
                var['t'].append(float(values[2]))

        elif values[0] == 'L2Norm(u)':
            if values[2] != '-nan' or values[2] != 'nan':
                var['L2Norm(u)'].append(float(values[2]))
    
        elif values[0] == 'dissip(u+U)':
            if values[2] != '-nan' or values[2] != 'nan':
                var['dissip(u+U)'].append(float(values[2]))


ut.message('Closing file')
f.close()
        

var['t']=np.asarray(var['t'])
var['L2Norm(u)']=np.asarray(var['L2Norm(u)'])
var['dissip(u+U)']=np.asarray(var['dissip(u+U)'])


plt.plot(var['t'], var['L2Norm(u)'], 'r-')
plt.grid(True)

plt.xlabel('Time units, t')
plt.ylabel('L2Norm(u)')

plt.show()
