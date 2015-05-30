"""
UTILITIES

This file contains miscellaneous functions, which are used for small
tasks such as printing errors. 

 Author: Muhammad Arslan Ahmed
  email: maa8g09@soton.ac.uk
  
University of Southampton
"""

import sys

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
    sys.exit('')
    return
    
    
def message(str):
    print '    Message:'
    print '   ', str, '\n'


def openFile(str):
    """
    INPUTS:
         str:  string of the directory where the file exists.
    OUTPUTS:
           f:  file object
    """
    
    f = open(str, 'r')
    print '\n    Opened the file:', str
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
          fourdarray: the plane co-ordinates of data you want to plot
      geom_variables: a dictionary contianing all the geometrical values
      
    OUTPUTS:
               -
               
    """
    
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
        print 'Invalid velocity component given, velocity component must be in range 0 to 2.'
        error('Invalid velocity component given!')
    
    print '\n    The 4D array input is valid.'
    return
    

def sparsify(skip, array):
    
    # Sparsify to reduce number of elements in a 1D array.
    #
    # This function controls the resolution of the plot.
    # The way it works is you pass in the array with a skip value.
    #
    # skip  = how many values to skip
    # array = The array to process
    #
    # For example:
    # array = [1,2,3,4,5,6,7,8]
    # skip  = 1
    # it becomes...
    # array = [1,3,5,7]

    tmp_array = []
    for i in range(0, len(array)):
        #print '',i,' i%(skip+1) = ',i%(skip+1)
        if (i % (skip + 1)) is 0:
            tmp_array.append(array[i])
    print '\n\n', tmp_array, '\n\n'
    
    return tmp_array
    