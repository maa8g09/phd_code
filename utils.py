"""
UTILITIES

This file contains miscellaneous functions, which are used for small tasks such 
as printing errors, messages and checking input validity.


Author details:
    Muhammad Arslan Ahmed
    maa8g09@soton.ac.uk
    
    Aerodynamics and Flight Mechanics Research Group
    Faculty of Engineering and the Environment
    University of Southampton
"""

import sys


def printStart():
    print '##################################################################'
    print '_________ .__                                .__            '
    print '\_   ___ \|  |__ _____    ____   ____   ____ |  |           '
    print '/    \  \/|  |  \\__  \  /    \ /    \_/ __ \|  |           '
    print '\     \___|   Y  \/ __ \|   |  \   |  \  ___/|  |__         '
    print ' \______  /___|  (____  /___|  /___|  /\___  >____/         '
    print '        \/     \/     \/     \/     \/     \/               '
    print '    __________                    .__                      __   '
    print '    \______   \ ____   __________ |  |___  __ ____   _____/  |_ '
    print '     |       _// __ \ /  ___/  _ \|  |\  \/ // __ \ /    \   __\.'
    print '     |    |   \  ___/ \___ (  <_> )  |_\   /\  ___/|   |  \  |  '
    print '     |____|_  /\___  >____  >____/|____/\_/  \___  >___|  /__|  '
    print '            \/     \/     \/                     \/     \/      '
    print 'Muhammad Arslan Ahmed'
    print 'University of Southampton'
    print '##################################################################'
    return


def printSectionHeader():
    print '__________________________________________________________________\n'
    return


def printSectionTitle(str):
    # Print the section headers for each main section of the code output.
    print ' **', str, '\n'
    return


def error(str):
    # Print the error and then exit from the program entirely
    print '!!!!====!!!!====!!!!====!!!!====!!!!===='
    print 'ERROR:', str
    print '!!!!====!!!!====!!!!====!!!!====!!!!===='
    print '\n\n\n'
    print '  ______    ______   __       __  ________       '
    print ' /      \  /      \ |  \     /  \|        \      '
    print '|  $$$$$$\|  $$$$$$\| $$\   /  $$| $$$$$$$$      '
    print '| $$ __\$$| $$__| $$| $$$\ /  $$$| $$__          '
    print '| $$|    \| $$    $$| $$$$\  $$$$| $$  \         '
    print '| $$ \$$$$| $$$$$$$$| $$\$$ $$ $$| $$$$$         '
    print '| $$__| $$| $$  | $$| $$ \$$$| $$| $$_____       '
    print ' \$$    $$| $$  | $$| $$  \$ | $$| $$     \      '
    print '  \$$$$$$  \$$   \$$ \$$      \$$ \$$$$$$$$      '                          
    print '                                                 '
    print '  ______   __     __  ________  _______          '
    print ' /      \ |  \   |  \|        \|       \         '
    print '|  $$$$$$\| $$   | $$| $$$$$$$$| $$$$$$$\        '
    print '| $$  | $$| $$   | $$| $$__    | $$__| $$        '
    print '| $$  | $$ \$$\ /  $$| $$  \   | $$    $$        '
    print '| $$  | $$  \$$\  $$ | $$$$$   | $$$$$$$\        '
    print '| $$__/ $$   \$$ $$  | $$_____ | $$  | $$        '
    print ' \$$    $$    \$$$   | $$     \| $$  | $$        '
    print '  \$$$$$$      \$     \$$$$$$$$ \$$   \$$        '
    sys.exit('')
    return
    
    
def message(str):
    print '   ', str, '\n'
    return


def openFile(str):
    """
    INPUTS:
         str:  string of the directory where the file exists.
    OUTPUTS:
           f:  file object
    """
    
    f = open(str, 'r')
    message('Opened the file: ' + str)
    return f
    

def printGeoVars(geom_variables):
    """
    This function prints out the grid dimensions of the channelflow solution.
    """
    print '\n    The geometry variables are:'
    print '      Nx:', geom_variables['Nx']
    print '      Ny:', geom_variables['Ny']
    print '      Nz:', geom_variables['Nz']
    
    return


def checkInputValidity(fourdarray, geom_variables):
    """
    This function checks that the values in the fourdarray are valid to allow
    the plotting of the slice specified. If it is valid the function executes
    as normal and a message is output in the shell to let the user know all's
    well.
    
    If, however, the values in the fourdarray are invalid, the whole routine 
    is stopped and an error message is output which specifies the invalidity
    of the input. 
    
    
    INPUTS:
        fourdarray: the plane co-ordinates of data you want to plot, indexed
                      as (i, nx, ny, nz)
      geom_variables: a dictionary contianing all the geometrical values
      
    OUTPUTS:
               
    """
    printSectionHeader()
    printSectionTitle('Checking validity of 4D array input')
    
    for i in range(1,4):
        if fourdarray[i] != 'all':
            if i == 1:
                if fourdarray[i] >= geom_variables['Nx']:
                    print '\n  ! The 4D array input is invalid'
                    printGeoVars(geom_variables)
                    error('X point given exceeds the maximum grid points available')
            elif i == 2:
                if fourdarray[i] >= geom_variables['Ny']:
                    print '\n  ! The 4D array input is invalid'
                    printGeoVars(geom_variables)
                    error('Y point given exceeds the maximum grid points available')
            else:
                if fourdarray[i] >= geom_variables['Nz']:
                    print '\n  ! The 4D array input is invalid'
                    printGeoVars(geom_variables)
                    error('Z point given exceeds the maximum grid points available')
    
    if fourdarray[0] < 0 or fourdarray[0] >= 3:
        message('Invalid velocity component given, velocity component must be in range 0 to 2.')
        error('Invalid velocity component given!')
    
    message('The 4D array input is valid.')
    return
    