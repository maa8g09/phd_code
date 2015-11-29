#!/usr/bin/env python
###############################################################################
##                                                                           ##
## Total process for finding TWs/RPOs                                        ##
## ========================================================================= ##
## The code presented provides a complete method of finding edge states      ##
## using velocity fields generated using the rsolvent formulation as         ##
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

import utils as ut
import generate_modes as gm
import dns as d
import postpro_standard as ps
import nkh as n
import argparse
import sys

parser = argparse.ArgumentParser(description='Total process to find edge states. \nPlease provide a directory to carry out the process in.')
parser.add_argument('-d', help='directory to output process in.', action='store')

args = parser.parse_args()
ut.printEdgeStatesStart()
print('\nThe process will be carried out in:')
print(args.directory, '\n')

print('=================================================')
print('\n#1 Generating the initial condition (using resolvent formulation)\n')
gm.generate_initial_condition(args.directory)

print('=================================================')
print('\n#2 Running DNS\n')
d.run_dns(args.directory)

print('=================================================')
print('\n#3 Post-processing the DNS simulation\n')
ps.run_postpro_standard(args.directory)

print('\nFinished!\n')
