"""
READ and CONSTRUCT

Read channelflow solutions and then use the data stored in those files to 
construct the solution as a 4D matrix.

The functions below check the directory where the files are kept and 
then reads them for contruction (construction of flow field should be a 4D array,
where (i, nx, ny, nz) will give you the i-th component of velocity at 
co-ordinate nx, ny, nz).


Author details:
    Muhammad Arslan Ahmed
    maa8g09@soton.ac.uk
    
    Aerodynamics and Flight Mechanics Research Group
    Faculty of Engineering and the Environment
    University of Southampton
"""

import math
import numpy as np
import utils as ut
import os


def main_read_construct(directory, fourdarray):
    """
    The main function that controls all reading and storing of flowfield. 

    
    INPUTS:
     directory:  the directory where the solution files are
    fourdarray:  the plane co-ordinates of data you want to plot
        
        
    OUTPUTS:
          data:  a dictionary containing the flowfield, geometry and (if physical
                 solution given) velocity slice data. This is then used later
                 for use in projection of resolvent modes.
    """

    data = {}
    
    
    # READ & CONSTRUCT
    # Check that the files exist in the directory specified.
    dict_files = check_files(directory)
    
    # Read and construct geometry dictionary
    dict_geometry = read_construct_geometry(directory, dict_files)
    
    # Check that the plane specified in fourdarray is valid in terms of the 
    # geometry of the channelflow solution read in.
    if dict_geometry['is_physical'] == True:
        ut.checkInputValidity(fourdarray, dict_geometry['physical'])

    # Read and construct flow field dictionary.
    dict_flowField = read_construct_flow_field(directory, dict_files, dict_geometry)

    # We need to plot a plane according to the fourdarray, if physical solution
    # is provided.
    if dict_geometry['is_physical'] == True:
        velslice = get_data_slice(dict_flowField['physical'], fourdarray, dict_geometry['physical'])
        data['velslice'] = velslice
    
    # Store the geometry and flow field information in a master dictionary which 
    # can be passed to functions easily.
    data['geometry'] = dict_geometry
    data['flowField'] = dict_flowField
    
    return data


def check_files(direc):
    """
     Check to see if the ASCII and geometry files exist. We are also going to
     see if the files are for physical or spectral solutions.
     
     If both exist, we can use both. Otherwise we will only use the one present,
     this presents different code paths to take later on in the code when we have
     to project the calculated resolvent modes onto the solution.
     
     If niether exist the function and the whole process will terminate.
     
    
    INPUTS:
         direc:  the directory where the solution files *should* exist
     
     
    OUTPUTS:
          a, b:  names of the solution files with extensions included
                   a - geometry file
                   b - ascii file with velocity data
    """
    
    ut.printSectionHeader()
    ut.printSectionTitle('Checking the ASCII and GEOM files')

    ut.message('The working directory is: ' + direc)
    ut.message('Checking to see if the necassary files exist in directory...')

    # These are strings which are used to check if the necassary files exist in
    # the directory. If at the end of the check these strings remain unchanged
    # the files DO NOT exist.
    a = 'geo'
    b = 'asc'
    c = 'spec.ge'
    d = 'spec.as'
    
    files = [fi for fi in os.listdir(direc) if os.path.isfile(os.path.join(direc,fi))]

    # Now we loop through the files in the directory to find which ones exist.
    for fi in files:
        if a in str(fi):
            a = str(fi) # physical geometry file
        if c in str(fi):
            c = str(fi) # spectral geometry file
        if b in str(fi):
            b = str(fi) # physical ascii file
        if d in str(fi):
            d = str(fi) # spectral ascii file



    files_in_direc = {}
    physical = False
    spectral = False
    missing_p_geo, missing_p_asc,  missing_s_geo, missing_s_asc = '', '', '', ''
    
    # Check to see if there are BOTH physical files (ASCII and GEO)
    if a is 'geo' and b is 'asc':
        ut.message('MISSING: Physical ASCII and geometry files.')
        
    elif a is 'geo' or b is 'asc':
        # if not, then point out which one is missing...DONT RAISE ERRORS...yet
        if a is 'geo':
            ut.message('MISSING: Could not find physical geometry file')
            missing_p_geo = '\nphysical geometry file'
            
        if b is 'asc':
            ut.message('MISSING: Could not find physical ASCII file')
            missing_p_asc = '\nphysical ascii file'
    else:
        physical = True
        files_in_direc['physical'] = physical
        files_in_direc['phy_geo'] = a
        files_in_direc['phy_asc'] = b
    
    
    # We carry out the same check for the spectral files.    
    if c is 'spec.ge' and d is 'spec.as':
        ut.message('MISSING: Spectral geometry and ASCII files are not in directory')
        
    elif c is 'spec.ge' or d is 'spec.as':
        # if not, then point out which one is missing...DONT RAISE ERRORS...yet
        if c is 'spec.ge':
            ut.message('MISSING: Could not find spectral geometry file')
            missing_s_geo = '\nspectral geometry file'
            
        if d is 'spec.as':
            ut.message('MISSING: Could not find spectral ASCII file')
            missing_s_asc = '\nspectral ascii file'
            
    else:
        spectral = True
        files_in_direc['spectral'] = spectral
        files_in_direc['spc_geo'] = c
        files_in_direc['spc_asc'] = d
        
        

    if spectral == True and physical == True:
        ut.message('We have both physical and spectral files.')
        
    elif spectral == False and physical == True:
        ut.message('We only have PHYSICAL')
        
    elif spectral == True and physical == False:
        ut.message('We only have SPECTRAL')
        
    else:
        # We are missing some files...
        errmsg = 'Missing the following files:' + missing_p_asc + missing_p_geo + missing_s_asc + missing_s_geo + '\n\nTerminating...'
        ut.error(errmsg)
        

    return files_in_direc


def read_construct_geometry(direc, dict_files):
    """
    Construct a dictionary with all the geometry variables inside it. The top-
    level keys will be physical/spectral, then each value will contain a further
    dictionary which will store all the key value pairs found in the geo file.
    
    
    INPUTS:
         direc:  the directory where the solution files are
    dict_files:  a dictionary with all the solution files in it.
     
     
    OUTPUTS:
       var_geo:  a dictionary of the geometrical variables that were contained
                 in the geometry file. The key and value is exactly the same as
                 in the file, i.e. the keys are the same labels as in the file.
    """

    dict_geometry = {}
    
    if dict_files['physical'] == True:
        # READ
        ut.printSectionHeader()
        ut.printSectionTitle('Reading the physical geometry file')
        
        f = ut.openFile(direc + '/' + dict_files['phy_geo'])
        
        ut.message('Constructing geometry dictionary')
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
    
        ut.message('Closing physical geometry file')
        f.close()
        
        
        var_geo['x'] = np.zeros((var_geo['Nx']))
        var_geo['y'] = np.zeros((var_geo['Ny']))
        var_geo['z'] = np.zeros((var_geo['Nz']))
        
        for nx in range(0, var_geo['Nx']):
            var_geo['x'][nx] = nx * var_geo['Lx'] / var_geo['Nx']
        for ny in range(0, var_geo['Ny']):
            var_geo['y'][ny] = math.cos(ny*math.pi/(var_geo['Ny']-1))
        for nz in range(0, var_geo['Nz']):
            var_geo['z'][nz] = nz * var_geo['Lz'] / var_geo['Nz']  
        
        
        dict_geometry['is_physical'] = True
        dict_geometry['physical'] = var_geo
    
    else:
        dict_geometry['is_physical'] = False
        
    
    if dict_files['spectral'] == True:
        # READ
        ut.printSectionTitle('Reading the spectral geometry file')
        
        f = ut.openFile(direc + '/' + dict_files['spc_geo'])
        
        ut.message('Constructing geometry dictionary')
        spec_geo = {}

        # CONSTRUCT
        for i, line in enumerate(f):
            values = line.split()
            
            if values[1] == '%kx':
                spec_geo['kx'] = int(values[0])
    
            if values[1] == '%kz':
                spec_geo['kz'] = int(values[0])
                
            if values[1] == '%y':
                spec_geo['Ny'] = int(values[0])
                
            if values[1] == '%Nd':
                spec_geo['Nd'] = int(values[0])
    
        ut.message('Closing spectral geometry file')
        f.close()
        
        dict_geometry['is_spectral'] = True
        dict_geometry['spectral'] = spec_geo
    
    else:
        dict_geometry['is_spectral'] = False


    return dict_geometry    


def read_construct_flow_field(direc, dict_files, dict_geometry):
    """
    Construct a dictionary with the contents of the ASCII files split into
    relevant physical/spectral keys.
    
    
    INPUTS:
         direc:  the directory where the solution files are
     dict_files:  a dictionary containing the file names of the ascii files
 dict_geometry:  a dictionary of the geometrical variables
     
     
    OUTPUTS:
 dict_flowField:  a dictionary of the flow field stored in physical and spectral
                 states.
                 
                     The physical flow field is stored as a 4D array: 
                      - (i, nx, ny, nz)
                      
                     The spectral flow field is stored as a 4D array:
                      - (i, kx, ny, kz)
    
    """

    dict_flowField = {}
    dict_flowField['is_physical'] = False
    dict_flowField['is_spectral'] = False
    
    if dict_files['physical'] == True:
        dict_flowField['is_physical'] = True
        # READ
        ut.printSectionHeader()
        ut.printSectionTitle('Reading physical ASCII file')
        f = ut.openFile(direc + '/' + dict_files['phy_asc'])
        
        # CONSTRUCT
        ut.message('Creating the flow field vector')
        
        # Create an empty 4D array to store the velocity data
        U = np.zeros((dict_geometry['physical']['Nd'], 
                      dict_geometry['physical']['Nx'], 
                      dict_geometry['physical']['Ny'], 
                      dict_geometry['physical']['Nz']))
        
        for i, line in enumerate(f):
            values = line.split()
            nx = int(values[0])
            ny = int(values[1])
            nz = int(values[2])
            nd = int(values[3])
            vel = float(values[4])
            
            U[nd, nx, ny, nz] = vel
            
        ut.message('Closing the physical ASCII file')
        f.close()
        
        dict_flowField['physical'] = U
        
        
    if dict_files['spectral'] == True:
        dict_flowField['is_spectral'] = True
        # READ
        ut.printSectionTitle('Reading spectral ASCII file')
        f = ut.openFile(direc + '/' + dict_files['spc_asc'])
        
        # CONSTRUCT
        ut.message('Creating the flow field vector')
        
        # Create an empty 4D array to store the velocity data
        U_hat = np.zeros((dict_geometry['spectral']['Nd'], 
                          dict_geometry['spectral']['kx'], 
                          dict_geometry['spectral']['Ny'], 
                          dict_geometry['spectral']['kz']), 
                          dtype=np.complex128)
                            
                            
        for i, line in enumerate(f):
            values = line.split()
            nd = int(values[0])
            kx = int(values[1])
            kz = int(values[2])
            ny = int(values[3])
            coeff = complex(values[4])
            
            U_hat[nd, kx, ny, kz] = coeff
    
        ut.message('Closing the spectral ASCII file')
        f.close()
        
        dict_flowField['spectral'] = U_hat
        
    
    return dict_flowField


def get_other_vels(nd):
    """
    The other velocity components are used to plot the quiver plots. So that 
    the velocity vectors of the arrows matches the axes that are being looked
    at. 
    
    
    INPUTS:
  four_d_array:  the velocity component passed in four_d_array

    
    OUTPUTS:
            v0:  the velocity of axis_0 (x axis of plot)
            v1:  the velocity of axis_1 (y axis of plot)
                  
    """
    
    if nd == 0:         # U
        v0, v1 = 2, 1   # W, V
        
    elif nd == 1:       # V
        v0, v1 = 0, 2   # U, W
        
    elif nd == 2:       # W
        v0, v1 = 0, 1   # U, V
        
    return v0, v1

    
def get_data_slice(ff, four_d_array, var_geo):
    """
    We extract the relevant data and setup a dictionary that can be used
    later to plot a velocity slice out of the 3D flowfield.
    
    
    INPUTS:
            ff:  a 4D matrix with the velocity data organized as (i, nx, ny, nz)
  four_d_array:  the plane co-ordinates of data you want to plot
       var_geo:  the dictionary with the details of the physical flow field
       
  
    OUTPUTS:
    slice_data:  a dictionary of all the necassary information require to plot
                 the velocity slice specified in four_d_array
    
    """
    
    ut.printSectionHeader()
    ut.printSectionTitle('Setup velocity slice for plotting')

    a0title = a1title = title = ''
    v0 = v1 = a0 = a1 = 0
    
    U_u2 = np.power(ff[0,:,:,:], 2)
    U_v2 = np.power(ff[1,:,:,:], 2)
    U_w2 = np.power(ff[2,:,:,:], 2)

    U_mag = U_u2 + U_v2 + U_w2
    
    U_mag = np.sqrt(U_mag)
    
    
    if four_d_array[1] != 'all': #ZY plane
        contour_data = U_mag[four_d_array[1], :, :]
        
        v0, v1 = get_other_vels(four_d_array[0])
        v0 = ff[v0, four_d_array[1], :, :]
        v1 = ff[v1, four_d_array[1], :, :]
        a0, a0title = var_geo['z'], 'z'
        a1, a1title = var_geo['y'], 'y'
        title = 'ZY plane at X: ' + str(four_d_array[1])

    elif four_d_array[2] != 'all': # XZ plane
        # we transpose the matrices because otherwise we would be plotting
        # X data on the y axis of the plot. 
        contour_data = U_mag[:, four_d_array[2], :]
        
        v0, v1 = get_other_vels(four_d_array[0])
        v0 = ff[v0, :, four_d_array[2], :].T
        v1 = ff[v1, :, four_d_array[2], :].T
        a0, a0title = var_geo['x'], 'x'
        a1, a1title = var_geo['z'], 'z'
        title = 'XZ plane at Y: ' + str(four_d_array[2])
        
    elif four_d_array[3] != 'all': # XY plane
        contour_data = U_mag[:, :, four_d_array[3]]
        
        v0, v1 = get_other_vels(four_d_array[0])
        v0 = ff[v0, :, :, four_d_array[3]]
        v1 = ff[v1, :, :, four_d_array[3]]
        a0, a0title = var_geo['x'], 'x'
        a1, a1title = var_geo['y'], 'y'
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
    