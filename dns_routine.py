#!/usr/bin/env python
import os
import subprocess as sp
from colorama import Fore, Back, Style

main_directory = '/home/arslan/Documents/work/channelflow-related/set01/Re1200/KB/ampls-DNS-2015_10_25-further'
os.chdir(main_directory)
pwd = os.getcwd()
print('Present Working Directory (after changing dir):\n')
print(pwd, '\n')
directories = [d for d in os.listdir(pwd) if os.path.isdir(d)]
directories = sorted(directories)


# Loop through the main directory and in each subfolder set-off a DNS simulation.
for i in range(0, len(directories)):
    print(Fore.BLACK + Back.GREEN + str(directories[i]) + Style.RESET_ALL,'\n')
    
    t_start   = 200;                print(Style.BRIGHT + '                T0:' + Style.RESET_ALL, t_start)
    t_end     = 800;                print(Style.BRIGHT + '                T1:' + Style.RESET_ALL, t_end)
    t_dt      = 0.001;              print(Style.BRIGHT + '                dt:' + Style.RESET_ALL, t_dt)
    t_dtmin   = 0.0001;             print(Style.BRIGHT + '            dt min:' + Style.RESET_ALL, t_dtmin)
    t_dtmax   = 0.01;               print(Style.BRIGHT + '            dt max:' + Style.RESET_ALL, t_dtmax)
    t_dtsave  = 0.1;                print(Style.BRIGHT + '                dT:' + Style.RESET_ALL, t_dtsave)
    print('')
    cfl_min   = 0.01;               print(Style.BRIGHT + '           CFL min:' + Style.RESET_ALL, cfl_min)
    cfl_max   = 1.0;                print(Style.BRIGHT + '           CFL max:' + Style.RESET_ALL, cfl_max)
    print('')
    nonlinearity = 'skew';          print(Style.BRIGHT + '      Nonlinearity:' + Style.RESET_ALL, nonlinearity) # Method of calculating nonlinearity
    print('')
    output_dir  = '/data-skew';     print(Style.BRIGHT + '  Output directory:' + Style.RESET_ALL, output_dir)
    output_file = 'data_skew.txt';  print(Style.BRIGHT + '       Output file:' + Style.RESET_ALL, output_file)
    print('')
    reynolds = 1200;                print(Style.BRIGHT + '          Reynolds:' + Style.RESET_ALL, reynolds)
    print(Style.BRIGHT + '                 b:' + Style.RESET_ALL,'Holding bulk velocity fixed')
    Ubulk = 1.33333333;             print(Style.BRIGHT + '             Ubulk:' + Style.RESET_ALL, Ubulk)
    print('')
    
    print('Executing the following command:')
    print('couette --channel -T0', t_start, '-T1', t_end, '-dt', t_dt, '-dtmin', t_dtmin, '-dtmax', t_dtmax,'-dT', t_dtsave, '-CFLmin', cfl_min, '-CFLmax', cfl_max, '-nl', nonlinearity, '-symms sigma.asc --outdir', output_dir, '-R', reynolds,'-b -U', Ubulk,'-cfl -l2 -D -I -dv -u -Up -p data-skew/u200.000.ff >>', output_file)
    print('')
    
    
    sp.call(['ls -l'])
#    sp.call(['couette', 
#    '--channel', 
#    '-T0', str(t_start),
#    '-T1', str(t_end),
#    '-dt', str(t_dt),
#    '-dtmin', str(t_dtmin),
#    '-dtmax', str(t_dtmax),
#    '-dT', str(t_dtsave),
#    '-CFLmin', str(cfl_min),
#    '-CFLmax', str(cfl_max),
#    '-nl', str(nonlinearity),
#    '-symms', 'sigma.asc',
#    '--outdir', output_dir,
#    '-R', str(reynolds),
#    '-b', 
#    '-U', str(Ubulk),
#    '-cfl', 
#    '-l2', 
#    '-D', 
#    '-I', 
#    '-dv', 
#    '-u',
#    '-Up', 
#    '-p',
#    'data-skew/u200.000.ff',
#    '>>', 
#    output_file
#    ])
    
    