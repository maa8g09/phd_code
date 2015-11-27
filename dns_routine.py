#!/usr/bin/env python
import os
import subprocess as sp
from colorama import Fore, Back, Style
import utils as ut


def main(case_direc, Re, symData):
    os.chdir(case_direc)
    pwd = os.getcwd()




    nonlinearity = 'skew'
    output_dir  = '/data-' + nonlinearity
    output_file = 'data_'+nonlinearity+'.txt'




    t_start   = 0
    if t_start == 0:
        initialFF = 'u0.ff'
    else:
        initialFF = output_dir + '/u'+str(t_start) + '.000.ff'
    t_end     = 800
    t_dt      = 0.001
    t_dtmin   = 0.0001
    t_dtmax   = 0.01
    t_dtsave  = 0.1

    cfl_min   = 0.01
    cfl_max   = 1.0

    reynolds = Re

    Ubulk = 1.33333333

               
    
    print(Style.BRIGHT + '                T0:' + Style.RESET_ALL, t_start)
    print(Style.BRIGHT + '                T1:' + Style.RESET_ALL, t_end)
    print(Style.BRIGHT + '                dt:' + Style.RESET_ALL, t_dt)
    print(Style.BRIGHT + '            dt min:' + Style.RESET_ALL, t_dtmin)
    print(Style.BRIGHT + '            dt max:' + Style.RESET_ALL, t_dtmax)
    print(Style.BRIGHT + '                dT:' + Style.RESET_ALL, t_dtsave)
    print('')
    print(Style.BRIGHT + '           CFL min:' + Style.RESET_ALL, cfl_min)
    print(Style.BRIGHT + '           CFL max:' + Style.RESET_ALL, cfl_max)
    print('')
    print(Style.BRIGHT + '      Nonlinearity:' + Style.RESET_ALL, nonlinearity) # Method of calculating nonlinearity
    print('')
    print(Style.BRIGHT + '  Output directory:' + Style.RESET_ALL, output_dir)
    print(Style.BRIGHT + '       Output file:' + Style.RESET_ALL, output_file)
    print('')
    print(Style.BRIGHT + '          Reynolds:' + Style.RESET_ALL, reynolds)
    print(Style.BRIGHT + '                 b:' + Style.RESET_ALL,'Holding bulk velocity fixed')
    print(Style.BRIGHT + '             Ubulk:' + Style.RESET_ALL, Ubulk)
    print('')
    
    
    print('Executing the following command:')
    command  = 'couette --channel'
    command += ' -T0 ' + t_start
    command += ' -T1 ' + t_end
    command += ' -dt ' + t_dt
    command += ' -dtmin ' + t_dtmin
    command += ' -dtmax ' + t_dtmax
    command += ' -dT ' + t_dtsave
    command += ' -CFLmin ' + cfl_min
    command += ' -CFLmax ' + cfl_max
    command += ' -nl ' + nonlinearity
    if symData['writtenSymmsFile']:
        command += ' -symms ' + symData['fileName']
    command += ' --outdir ' + output_dir
    command += ' -R' + reynolds
    command += ' -b -U ' + Ubulk
    command += ' -cfl -l2 -D -I -dv -u -Up -p'
    command += ' ' + initialFF 
    command += ' >> ' + output_file
    print('')

    
    
