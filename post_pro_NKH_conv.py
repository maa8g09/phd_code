#!/usr/bin/env python
import os
import subprocess as sp
import read_construct as rc
import utils as ut
import numpy as np
from pylab import *
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
#from matplotlib import animation
from matplotlib import cm as cm
from mpl_toolkits import mplot3d as axes3D
####################################################################################################
#### Post-processing the newton search cases...                                                 ####
####################################################################################################

# First we have to go through and find all newton steps and put them into their own folders...
case_direc = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further/wavepacket_007_4modes_(-1.51191822883+0j)/nkh-vdt/'

conv_file = 'convergence.asc'

conv_file = ut.openFile(case_direc + conv_file)
conv_dict = {}
conv_dict['L2Norm(G)'] = []
conv_dict['r'] = []
conv_dict['delta'] = []
conv_dict['dT'] = []
conv_dict['L2Norm(du)'] = []
conv_dict['L2Norm(u)'] = []
conv_dict['L2Norm(dxN)'] = []
conv_dict['dxNalign'] = []
conv_dict['L2Norm(dxH)'] = []
conv_dict['GMRESerr'] = []
conv_dict['ftotal'] = []
conv_dict['fnewt'] = []
conv_dict['fhook'] = []
conv_dict['CFL'] = []
conv_dict['ns'] = []

for i, line in enumerate(conv_file):
    values = line.split()
    if i != 0:
        conv_dict['L2Norm(G)'].append(values[0])
        conv_dict['r'].append(values[1])
        conv_dict['delta'].append(values[2])
        conv_dict['dT'].append(values[3])
        conv_dict['L2Norm(du)'].append(values[4])
        conv_dict['L2Norm(u)'].append(values[5])
        conv_dict['L2Norm(dxN)'].append(values[6])
        conv_dict['dxNalign'].append(values[7])
        conv_dict['L2Norm(dxH)'].append(values[8])
        conv_dict['GMRESerr'].append(values[9])
        conv_dict['ftotal'].append(values[10])
        conv_dict['fnewt'].append(values[11])
        conv_dict['fhook'].append(values[12])
        conv_dict['CFL'].append(values[13])
        conv_dict['ns'].append(i-1)


plt.plot(conv_dict['ns'], 
         conv_dict['L2Norm(u)'])
         
plt.xlabel('Newton steps')
plt.ylabel('||u(t)||')

plt.grid(True)

plt.show()
