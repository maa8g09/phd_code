"""
READ and CONSTRUCT

Read channelflow solutions and then use the data stored in those files to 
construct the solution.

The functions below check the directory where the files (should) be kept and 
then reads them for contruction (construction of flow field should be a 4D array,
where (nx, ny, nz, i) will give you the i-th conponent of velocity at 
co-ordinate nx, ny, nz).


    Author: Muhammad Arslan Ahmed
     email: maa8g09@soton.ac.uk
  
University of Southampton
"""

import math
import numpy as np
import utils as ut
import os


def main_read_construct(directory, fourdarray, verbosity):
    """
    The main function that handles all of the reading and constructing of the 
    flowfield. 

    
    INPUTS:
     directory:  the directory where the solution files are
    fourdarray:  the plane co-ordinates of data you want to plot
     verbosity:  the verbosity of the output in the python shell
        
        
    OUTPUTS:
           geo:  a dictionary contianing all the values from *.geo file
       flowfield:  a dictionary containing the geometrical co-ordinates as 
                 well as velocity vectors as 3D arrays at each velocity 
                 component
      vecslice:  a dictionary containing information necassary to plot the 
                 plane specified using the fourdarray input. 
    """
    
    # READ & CONSTRUCT
    # Check that the files exist in the directory specified.
    f_geo, f_ascii = check_files(directory, verbosity)
    # Read and construct geometry dictionary or channelflow solution
    geo = read_construct_geometry(directory, f_geo, verbosity)
    # Check that the plane specified in fourdarray is valid in terms of the 
    # geometry of the channelflow solution read in.
    ut.checkInputValidity(fourdarray, geo)
    # We can only continue if the input is valid so to get this far means that
    # the input was valid.
    # Read and construct flow field dictionary or channelflow solution.
    flowfield = read_construct_flow_field(directory, f_ascii, geo, verbosity)
    # We need to plot a plane according to the elements of fourdarray
    velslice = get_data_slice(flowfield, fourdarray, geo)
    
    data = {
            'geometry': geo,
            'flowfield': flowfield,
            'velslice': velslice
            }
    
    return data


def check_files(direc, verbosity):
    """
    INPUTS:
         direc:  the directory where the solution files *should* exist
     verbosity:  boolean that controls how verbose the output of the 
                 function is
     
     
    OUTPUTS:
          a, b:  names of the solution files with extensions included
                   a - geometry file
                   b - ascii file with velocity data
    """
    
    ut.printSectionHeader()
    ut.printSectionTitle('Checking the ASCII and GEOM files')

    if verbosity:
        print '\n        The working directory is:    ', direc
        print '\n        Checking to see if the necassary files exist in directory...'

    # These are strings which are used to check if the necassary files exist in
    # the directory. If at the end of the check these strings remain unchanged
    # the files DO NOT exist.
    a = 'geo'
    b = 'asc'
    c = 'ec.ge'
    d = 'ec.as'
    files = [fi for fi in os.listdir(direc) if os.path.isfile(os.path.join(direc,fi))]
    for fi in files:
        if a in str(fi):
            if c in str(fi):
                # spectral geometry file
                c = str(fi)
            else:
                a = str(fi)
            
        if b in str(fi):
            if d in str(fi):
                #spectral ascii file
                d = str(fi)
            else:
                b = str(fi)

    if a is 'geo':
        ut.error('Could not find physical geometry file')

    if b is 'asc':
        ut.error('Could not find physical ASCII file')
        
    if c is 'ec.ge':
        ut.error('Could not find spectral geometry file')
        
    if d is 'ec.as':
        ut.error('Could not find spectral ASCII file')

    # If at this stage, no errors have occured and the relevant files are 
    # in the directory...
    if verbosity:
        print '\n        Necassary files exist, we are good to go!'



    files_in_direc = {}
    physical = False
    spectral = False
    

    # Check to see if there are BOTH physical
    if a is 'geo' and b is 'asc':
        # both physical files are missing.
        physical = False
        ut.message('MISSING: Physical geometry and ASCII files are not in directory')
    elif a is 'geo' or b is 'asc':
        # if not, then point out which one is missing...DONT RAISE ERRORS
        if a is 'geo':
            ut.message('MISSING: Could not find physical geometry file')
        if b is 'asc':
            ut.message('MISSING: Could not find physical ASCII file')
    else:
        physical = True
        files_in_direc['p_geo'] = a
        files_in_direc['p_asc'] = b
    
    
    if c is 'ec.ge' and d is 'ec.as':
        # both spectral files are missing.
        spectral = False
        ut.message('MISSING: Spectral geometry and ASCII files are not in directory')
    elif c is 'ec.ge' or d is 'ec.as':
        # if not, then point out which one is missin... DONT RAISE ERRORS
        if c is 'ec.ge':
            ut.message('MISSING: Could not find spectral geometry file')
        if d is 'ec.as':
            ut.message('MISSING: Could not find spectral ASCII file')
    else:
        spectral = True
        files_in_direc['s_geo'] = c
        files_in_direc['s_asc'] = d
        

    if spectral == True and physical == True:
        ut.message('We have both physical and spectral files.')
        
    elif spectral == False and physical == True:
        ut.message('We only have PHYSICAL')
        
    elif spectral == True and physical == False:
        ut.message('We only have SPECTRAL')
        
    else:
        ut.error('No solution files found, terminating read...')
            


    return files_in_direc


def read_construct_geometry(direc, f_geo, verbosity):
    """
    INPUTS:
         direc:  the directory where the solution files are
         f_geo:  geometry file
     verbosity:  boolean that controls how verbose the output of the function is 
     
    OUTPUTS:
       var_geo:  a dictionary of the geometrical variables that were contained
                 in the geometry file. The key and value is exactly the same as
                 in the file, i.e. the keys are the same labels as in the file.
    """

    
    # READ
    ut.printSectionHeader()
    ut.printSectionTitle('Reading the geometry file')
    
    f = ut.openFile(direc + '/' + f_geo)
    
    print '\n    Constructing geometry dictionary'
    var_geo = {}
    
    
    # CONSTRUCT
    for i, line in enumerate(f):
        values = line.split()
        
        if values[1] == '%Nx':
            var_geo['Nx'] = int(values[0])

        elif values[1] == '%Ny':
            var_geo['Ny'] = int(values[0])

        elif values[1] == '%Nz':
            var_geo['Nz'] = int(values[0])

        elif values[1] == '%Nd':
            var_geo['Nd'] = int(values[0])

        elif values[1] == '%Lx':
            var_geo['Lx'] = float(values[0])

        elif values[1] == '%Lz':
            var_geo['Lz'] = float(values[0])

        elif values[1] == '%lx=Lx/(2pi)':
            var_geo['lx'] = float(values[0])

        elif values[1] == '%lz=Lz/(2pi)':
            var_geo['lz'] = float(values[0])

        elif values[1] == '%alpha=2pi/Lx':
            var_geo['alpha'] = float(values[0])

        elif values[1] == '%gamma=2pi/Lz':
            var_geo['gamma'] = float(values[0])

    print '\n    Closing the geometry file'
    f.close()
    
    if verbosity:
        print '\n    ', var_geo

    return var_geo    


def read_construct_flow_field(direc, f_ascii, geom, verbosity):
    """
    INPUTS:
         direc:  the directory where the solution files are
       f_ascii:  ascii file containing flow field data
          geom:  a dictionary of the geometrical variables
     verbosity:  boolean that controls how verbose the output of the function 
                 is 
     
     
    OUTPUTS:
            ff:  a dictionary of the flow-field and the grid point variables.
                 The 4D array which contains the flow field variables, the way 
                 the data is indexed is: (i, nx, ny, nz). i is the i-th compo-
                 nent of velocity and then the rest is self-explanatory, i.e.
                 the i-th component of velocity at the nx,ny,nz grid point.
    
    """
    
    
    # READ
    ut.printSectionHeader()
    ut.printSectionTitle('Reading the ASCII file')
    f = ut.openFile(direc + '/' + f_ascii)
    
    print '\n    Making a list out of values contained in the ASCII file'
    var_ff_list = []
    
    for i, line in enumerate(f):
        var_ff_list.append(float(line))

    print '\n    Closing the ASCII file'
    f.close()
    
    
    # CONSTRUCT
    ut.printSectionHeader()
    ut.printSectionTitle('Creating the flow field vector')
    
    # Create an empty 4D array to store the velocity 
    var_ff = np.zeros((geom['Nd'], geom['Nx'], geom['Ny'], geom['Nz']))
    
    allU = []
    allV = []
    allW = []
    
    # The flow_variables[0::3] means create subset collection of elements that 
    # (index % 3 is 0) i.e. the u variables
    
    # The flow_variables[1::3] means create subset collection of elements that 
    # (index % 3 is 1) i.e. the v variables
    
    # The flow_variables[2::3] means create subset collection of elements that 
    # (index % 3 is 2) i.e. the v variables
    
    for u, v, w in zip(var_ff_list[0::3], var_ff_list[1::3], var_ff_list[2::3]):
        if u <= 1.0e-12:
            u = 0.0
        if v <= 1.0e-12:
            v = 0.0
        if w <= 1.0e-12:
            w = 0.0
        allU.append(u)
        allV.append(v)
        allW.append(w)
        
#        
#    
#    var_ff_list_len_third = len(var_ff_list) / 3
#    
#    allU = var_ff_list[0:var_ff_list_len_third]
#    allV = var_ff_list[var_ff_list_len_third:2*var_ff_list_len_third]
#    allW = var_ff_list[2*var_ff_list_len_third:3*var_ff_list_len_third]
#    
    
    allU = np.asarray(allU)
    allV = np.asarray(allV)
    allW = np.asarray(allW)
    
    var_ff[0,:,:,:] = allU.reshape((geom['Nx'], geom['Ny'], geom['Nz']))
    var_ff[1,:,:,:] = allV.reshape((geom['Nx'], geom['Ny'], geom['Nz']))
    var_ff[2,:,:,:] = allW.reshape((geom['Nx'], geom['Ny'], geom['Nz']))
    
    
    x_grid_points = np.zeros((geom['Nx']))
    for nx in range(0, geom['Nx']):
        x_grid_points[nx] = nx * geom['Lx'] / geom['Nx']
        
        
    y_grid_points = np.zeros((geom['Ny']))
    for ny in range(0, geom['Ny']):
        y_grid_points[ny] = math.cos(ny*math.pi/(geom['Ny']-1))
        
        
    z_grid_points = np.zeros((geom['Nz']))
    for nz in range(0, geom['Nz']):
        z_grid_points[nz] = nz * geom['Lz'] / geom['Nz']
    
    
    ff = {'X': x_grid_points,
          'Y': y_grid_points,
          'Z': z_grid_points,
          'flow_field': var_ff}
    
    
    return ff


def get_other_vels(four_d_array_vel):
    """
    The other velocity components are used to plot the quiver plots. So that 
    the velocity vectors of the arrows matches the axes that are being looked
    at. 
    
    INPUTS:
    four_d_array_vel:  the velocity component passed in four_d_array
    
    OUTPUTS:
                  v0:  the velocity of axis_0 (x axis of plot)
                  v1:  the velocity of axis_1 (y axis of plot)
    """
    if four_d_array_vel == 0: # U
        v0, v1 = 2, 1 # W, V
    elif four_d_array_vel == 1: # V
        v0, v1 = 0, 2 # U, W
    elif four_d_array_vel == 2: # W
        v0, v1 = 0, 1 # U, V
        
    return v0, v1
    


   
#def get_data_slice(ff, four_d_array, var_geo):
#    """
#    INPUTS:
#            ff:  a dictionary of the flow-field and the grid point variables.
#  four_d_array:  the plane co-ordinates of data you want to plot
#  
#  
#    OUTPUTS:
#    slice_data:  a dictionary of all the necassary information require to plot
#                 the velocity slice specified in four_d_array
#    
#    """
#    ut.printSectionHeader()
#    ut.printSectionTitle('Get the velocity plane and set-up for plotting')
#    
#    a0title = a1title = title = ''
#    v0 = v1 = a0 = a1 = 0
#    
#
#
##    U_mag = np.zeros(var_geo['Nx'], var_geo['Ny'], var_geo['Nz'])
##    
##    u2 = ff['flow_field'][0,:,:,:] * ff['flow_field'][0,:,:,:]
##    v2 = ff['flow_field'][1,:,:,:] * ff['flow_field'][1,:,:,:]
##    w2 = ff['flow_field'][2,:,:,:] * ff['flow_field'][2,:,:,:]
##    
##    U2 = u2 + v2 + w2
##    U_mag = np.sqrt(U2)
#    
#    
#    
#    
#    
#    if four_d_array[1] != 'all': #ZY plane
#        contour_data = ff['flow_field'][four_d_array[0], four_d_array[1], :, :]
#        v0, v1 = get_other_vels(four_d_array[0])
#        v0 = ff['flow_field'][v0, four_d_array[1], :, :]
#        v1 = ff['flow_field'][v1, four_d_array[1], :, :]
#        a0, a0title = ff['Z'], 'z'
#        a1, a1title = ff['Y'], 'y'
#        title = 'ZY plane at X: ' + str(four_d_array[1])
#
#    elif four_d_array[2] != 'all': # XZ plane
#        # we transpose the matrices because otherwise we would be plotting
#        # X data on the y axis of the plot. 
#        contour_data = ff['flow_field'][four_d_array[0], :, four_d_array[2], :].T
#        v0, v1 = get_other_vels(four_d_array[0])
#        v0 = ff['flow_field'][v0, :, four_d_array[2], :].T
#        v1 = ff['flow_field'][v1, :, four_d_array[2], :].T
#        a0, a0title = ff['X'], 'x'
#        a1, a1title = ff['Z'], 'z'
#        title = 'XZ plane at Y: ' + str(four_d_array[2])
#        
#    elif four_d_array[3] != 'all': # XY plane
#        contour_data = ff['flow_field'][four_d_array[0], :, :, four_d_array[3]]
#        v0, v1 = get_other_vels(four_d_array[0])
#        v0 = ff['flow_field'][v0, :, :, four_d_array[3]]
#        v1 = ff['flow_field'][v1, :, :, four_d_array[3]]
#        a0, a0title = ff['X'], 'x'
#        a1, a1title = ff['Y'], 'y'
#        title = 'XY plane at Z: ' + str(four_d_array[3])
#
# 
#    data_slice = {'contourData': contour_data,
#                  'axis_0': a0,
#                  'axis_1': a1,
#                  'axis_0_title': a0title,
#                  'axis_1_title': a1title,
#                  'vel_0': v0,
#                  'vel_1': v1,
#                  'plotTitle': title}
#    
#    return data_slice
    
def get_data_slice(ff, four_d_array, var_geo):
    """
    INPUTS:
            ff:  a dictionary of the flow-field and the grid point variables.
  four_d_array:  the plane co-ordinates of data you want to plot
  
  
    OUTPUTS:
    slice_data:  a dictionary of all the necassary information require to plot
                 the velocity slice specified in four_d_array
    
    """
    ut.printSectionHeader()
    ut.printSectionTitle('Get the velocity plane and set-up for plotting')
    
    a0title = a1title = title = ''
    v0 = v1 = a0 = a1 = 0
    


    U_mag = np.zeros((var_geo['Nx'], var_geo['Ny'], var_geo['Nz']))
    U_mag = (ff['flow_field'][0,:,:,:] * ff['flow_field'][0,:,:,:]) + (ff['flow_field'][1,:,:,:] * ff['flow_field'][1,:,:,:]) + ff['flow_field'][2,:,:,:] * ff['flow_field'][2,:,:,:]
    U_mag = np.sqrt(U_mag)
    
    
    
    
    
    if four_d_array[1] != 'all': #ZY plane
        contour_data = U_mag[four_d_array[1], :, :]
        
        v0, v1 = get_other_vels(four_d_array[0])
        v0 = ff['flow_field'][v0, four_d_array[1], :, :]
        v1 = ff['flow_field'][v1, four_d_array[1], :, :]
        a0, a0title = ff['Z'], 'z'
        a1, a1title = ff['Y'], 'y'
        title = 'ZY plane at X: ' + str(four_d_array[1])

    elif four_d_array[2] != 'all': # XZ plane
        # we transpose the matrices because otherwise we would be plotting
        # X data on the y axis of the plot. 
        contour_data = U_mag[:, four_d_array[2], :]
        
        v0, v1 = get_other_vels(four_d_array[0])
        v0 = ff['flow_field'][v0, :, four_d_array[2], :].T
        v1 = ff['flow_field'][v1, :, four_d_array[2], :].T
        a0, a0title = ff['X'], 'x'
        a1, a1title = ff['Z'], 'z'
        title = 'XZ plane at Y: ' + str(four_d_array[2])
        
    elif four_d_array[3] != 'all': # XY plane
        contour_data = U_mag[:, :, four_d_array[3]]
        
        v0, v1 = get_other_vels(four_d_array[0])
        v0 = ff['flow_field'][v0, :, :, four_d_array[3]]
        v1 = ff['flow_field'][v1, :, :, four_d_array[3]]
        a0, a0title = ff['X'], 'x'
        a1, a1title = ff['Y'], 'y'
        title = 'XY plane at Z: ' + str(four_d_array[3])

 
 
    data_slice = {'contourData': contour_data,
                  'axis_0': a0,
                  'axis_1': a1,
                  'axis_0_title': a0title,
                  'axis_1_title': a1title,
                  'vel_0': v0,
                  'vel_1': v1,
                  'plotTitle': title}
    
    return data_slice
    