#!/usr/bin/env python
###############################################################################
##                                                                           ##
## Total process for finding TWs/RPOs                                        ##
## ========================================================================= ##
## The code presented provides a complete method of finding edge states      ##
## using velocity fields generated using the rsolvent formulation as         ##
## an initial condition.                                                     ##
##                                                                           ##
## The program works by simply providing a configuration file (there is a    ##
## configuration file maker provided) and the output directory where the     ##
## process is output.                                                        ##
##                                                                           ##
## The resolvent formulation was developed by Sharma & McKeon (2010).        ##
##                                                                           ##
## Author:       Muhammad Arslan Ahmed                                       ##
## Emal:         maa8g09@soton.ac.uk                                         ##
## Institute:    University of Southampton                                   ##
## Department:   Aerodynamics and Flight Mechanics                           ##
##                                                                           ##
###############################################################################

import sys
import os
import argparse
pwd = os.getcwd()

sys.path.insert(0, pwd+'/modes')
sys.path.insert(1, pwd+'/dns')
sys.path.insert(2, pwd+'/post_process')
sys.path.insert(3, pwd+'/nkh')
sys.path.insert(4, pwd+'/misc')

import utils as ut
import generate_modes as gm
import dns
import postpro_standard as ps
#import nkh
from datetime import datetime
import time

parser = argparse.ArgumentParser(description='Total process to find edge states. \nPlease provide a directory to carry out the process in.')
parser.add_argument('-d', help='directory to output process in.', action='store')

args = parser.parse_args()
ut.printEdgeStatesStart()
print('\nThe process will be carried out in:')
print(args.directory, '\n')


startTime_1 = datetime.now()
print('=================================================')
print('\n#1 Generating the initial condition (using resolvent formulation)\n')
gm.generate_initial_condition(args.directory)
endTime_1 = datetime.now() - startTime_1

startTime_2 = datetime.now()
print('=================================================')
print('\n#2 Running DNS\n')
dns.run(args.directory)
endTime_2 = datetime.now() - startTime_2

startTime_3 = datetime.now()
print('=================================================')
print('\n#3 Post-processing the DNS simulation\n')
ps.run_postpro_standard(args.directory)
endTime_3 = datetime.now() - startTime_3

print('\nFinished!\n')
