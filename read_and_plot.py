#!/usr/bin/env python

import os
import sys
import math
import numpy as np
from pylab import *
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
#from matplotlib import animation
from matplotlib import cm as cm
from mpl_toolkits import mplot3d as axes3D

"""

Make a 2D array based on how many lines are in the files and how many variables
are on each line.

"""

def openFile(str):
    """
    INPUTS:
         str:  string of the directory where the file exists.
    OUTPUTS:
           f:  file object
    """

    f = open(str, 'r')
#    message('Opened the file: ' + str)
    return f


# Read in the file
distance_file='/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_008_4modes_(-0.931112136502+0j)/seriesdist-output-500-999-8-8.asc'

direct='/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_008_4modes_(-0.931112136502+0j)/'

f=openFile(distance_file)

data=[]

perLine = 0

for j, line in enumerate(f):
    
    
    if j==0:
        # do nothing
        print('nothing to do here...')
        
    else:
        # do something
        values = line.split()
        vl_len=len(values)
        
        perLine = vl_len
        
        for i in range(0,vl_len):
            data.append(values[i])
        


        

print('storing finished')

data = np.asarray(data, dtype = np.float128)

lines = len(data) / float(perLine)
lines = int(lines)
data = data.reshape((perLine, lines))





# window data:
data_plot = data





# Now we print.
y=np.arange(data_plot.shape[0])
x=np.arange(data_plot.shape[1])
x=x+500
x, y = np.meshgrid(x, y)
origin = 'lower'




v = 7
v = np.linspace(0.0, 0.0015, 7, endpoint=True)

CS = plt.contourf(x, y, data_plot, v, cmap=cm.jet, origin=origin)
cbar = plt.colorbar(CS)
plt.xlabel('t')
plt.ylabel('T')
plt.title('L2norm( u(t)  -  u(t+T) )')
plt.grid(True)

plt.xlim([500,900])
plt.ylim([0,60])

plt.savefig(direct+'recurrancePlot.png')









#
#
## plot period ... 
#fig = plt.figure(figsize=plt.figaspect(1.0))
#plot1 = fig.add_subplot(1, 1, 1)
#
#
#first = data[:, 10]
#y=np.arange(data.shape[0])
#x=np.arange(len(first))
#
#plot1.grid(True)
#plot1.plot(y, first, 'r-')
#xticks = np.linspace(0, 100, 25)
##yticks = np.linspace(y[0], y[-1], 101)
#plot1.set_xticks(xticks)
##plot1.set_yticks(yticks)
#plt.show
#
#
#




