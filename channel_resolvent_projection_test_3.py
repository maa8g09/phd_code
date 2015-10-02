import tools_pseudospectral as ps
#import utils_plots as up
import utils as ut
import math
import numpy as np
import numpy.fft as fourier
from scipy.linalg import inv
from scipy.linalg import solve
from scipy.linalg import pinv
from scipy.linalg import svd
from pylab import *
from matplotlib import pyplot as plt
from matplotlib import animation


np.set_printoptions(precision=4)

"""

The inputs for the "main_resolvent_analysis" function are given at the bottom of
the file.

"""


def main_resolvent_analysis(N, Re, kx, kz, c):
    """
    The full resolvent formulation is given here, from the generation of the
    modes to the projection of a channelflow solution onto the modes.

    INPUTS:
             N:  resolution in y axis
            Re:  Reynolds number
            kx:  vector of streamwise wavenumbers
            kz:  vector of spanwise wavenumbers
             c:  phase speed

    OUTPUTS:

    """


    t = 0





    # read in the spectral flowfield:
    
    # Mac
    direc = '/Users/arslan/Desktop/phd_coding/solutions/equilibria/nagata1' 
    
    # Linux
    direc = '/home/arslan/Documents/phd/code/channelflow-1.4.2/solutions/equilibria/nagata1'
    

    # read the spectral geometry file:
    f = ut.openFile(direc + '/eq1_spec.geom')
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

    f.close()
    
    
    kx = np.arange((spec_geo['kx'] / -2.0) + 1 , (spec_geo['kx'] / 2.0) + 1)
    kz = np.arange(0,spec_geo['kz'])


    # CONSTRUCT the flow field
    f = ut.openFile(direc + '/eq1_spec.asc')
    
    # a 4D array to store it all:
    u_hat = np.zeros((spec_geo['Nd'], spec_geo['kx'], spec_geo['Ny'], spec_geo['kz']), dtype=np.complex128)
                        
                        
    for i, line in enumerate(f):
        values = line.split()
        
        # Nd
        nd = int(values[0])
        
        # kx
        alpha = int(values[1])
        
        # kz
        beta = int(values[2])
        
        # y
        why = int(values[3])
        
        # coefficient
        cef = complex(values[4])
        
        u_hat[nd, alpha, why, beta] = cef

    u_hat_u = np.zeros((spec_geo['kx'], spec_geo['Ny'], spec_geo['kz']), dtype=np.complex128)
    u_hat_v = np.zeros((spec_geo['kx'], spec_geo['Ny'], spec_geo['kz']), dtype=np.complex128)
    u_hat_w = np.zeros((spec_geo['kx'], spec_geo['Ny'], spec_geo['kz']), dtype=np.complex128)

    u_hat_u = u_hat[0, :, :, :]
    u_hat_v = u_hat[1, :, :, :]
    u_hat_w = u_hat[2, :, :, :]


    f.close()
    # here we have to read in the physical ascii files.
    # what we can do is, fft the ascii file and see if we get the same coeffs as the channelflow code.
    geom = {}
    
    f = ut.openFile(direc + '/eq1.geom')
    for i, line in enumerate(f):
        values = line.split()
        
        if values[1] == '%Nx':
            geom['Nx'] = int(values[0])

        elif values[1] == '%Ny':
            geom['Ny'] = int(values[0])

        elif values[1] == '%Nz':
            geom['Nz'] = int(values[0])

        elif values[1] == '%Nd':
            geom['Nd'] = int(values[0])

        elif values[1] == '%Lx':
            geom['Lx'] = float(values[0])

        elif values[1] == '%Lz':
            geom['Lz'] = float(values[0])

        elif values[1] == '%lx=Lx/(2pi)':
            geom['lx'] = float(values[0])

        elif values[1] == '%lz=Lz/(2pi)':
            geom['lz'] = float(values[0])

        elif values[1] == '%alpha=2pi/Lx':
            geom['alpha'] = float(values[0])

        elif values[1] == '%gamma=2pi/Lz':
            geom['gamma'] = float(values[0])


    f.close()
    
    f = ut.openFile(direc + '/eq1.asc')
    var_ff_list = []
    
    for i, line in enumerate(f):
        var_ff_list.append(float(line))

    print('\n    Closing the ASCII file')
    f.close()
    
    
    # CONSTRUCT
    ut.printSectionHeader()
    ut.printSectionTitle('Creating the flow field vector')
    
    # Create an empty 4D array to store the velocity 
    var_ff = np.zeros((geom['Nd'], geom['Nx'], geom['Ny'], geom['Nz']))
    
    allU = []
    allV = []
    allW = []

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
        
    
    allU = np.asarray(allU)
    allV = np.asarray(allV)
    allW = np.asarray(allW)
    
    var_ff[0,:,:,:] = allU.reshape((geom['Nx'], geom['Ny'], geom['Nz']))
    var_ff[1,:,:,:] = allV.reshape((geom['Nx'], geom['Ny'], geom['Nz']))
    var_ff[2,:,:,:] = allW.reshape((geom['Nx'], geom['Ny'], geom['Nz']))
    var_ff_u = var_ff[0,:,:,:]
    allU = allV = allW = []
    
    resultant_vel = np.zeros((geom['Nx'], geom['Ny'], geom['Nz']))
    for nx in range(0, geom['Nx']):
        for ny in range(0, geom['Ny']):
            for nz in range(0, geom['Nz']):
                a2 = math.pow(var_ff[0,nx,ny,nz], 2.0)
                b2 = math.pow(var_ff[1,nx,ny,nz], 2.0)
                c2 = math.pow(var_ff[2,nx,ny,nz], 2.0)
                resultant_vel[nx,ny,nz] = math.sqrt((a2) + (b2) + (c2))
    
    
    x_grid_points = np.zeros((geom['Nx']))
    for nx in range(0, geom['Nx']):
        x_grid_points[nx] = nx * geom['Lx'] / geom['Nx']
        
        
    y_grid_points = np.zeros((geom['Ny']))
    for ny in range(0, geom['Ny']):
        y_grid_points[ny] = math.cos(ny*math.pi/(geom['Ny']-1))
        
        
    z_grid_points = np.zeros((geom['Nz']))
    for nz in range(0, geom['Nz']):
        z_grid_points[nz] = nz * geom['Lz'] / geom['Nz']
    
    

    # now that we have the physical flow field and the spectral flowfield. we can decompose the
    # physical ff to see if we get the same answers as the spectral ff that was read in.
    # so...
    transformed_signal = np.zeros((spec_geo['kx'], 
                                   spec_geo['Ny'],
                                   spec_geo['kz']),
                                   dtype=np.complex128)
    
    transformed_signal2 = np.zeros((geom['Nx'], 
                                    geom['Ny'],
                                    geom['Nz']),
                                    dtype=np.complex128)
    
    
    for i_kx in range(0, len(kx)):          # for every stream-wise wavenumber
        for i_kz in range(0, len(kz)):          # for every span-wise wavenumber        
            for i_c in range(0, len(c)):            # for every phase speed
#                for iy in range(0, geom['Ny']):
                    # U
                tmpu = var_ff_u[i_kx,:,i_kz]
                transformed_signal[i_kx, :, i_kz] = fft(tmpu, kx[i_kx], kz[i_kz], c[i_c], x_grid_points, z_grid_points, t, geom['Lx'], geom['Lz'])
                
                print('kx:', kx[i_kx], '    kz:', kz[i_kz], '    c:', c[i_c])
                    
    arslan = transformed_signal
    
    test01 = u_hat[:,:,11]
    test02 = transformed_signal[:,:,11]

    # now I am going to take the 2D fourier transform using python numpy library.
    for iy in range(0, geom['Ny']):
        xz_plane = resultant_vel[:,iy,:]
        xz_plane_fft = fourier.fft(xz_plane)
        transformed_signal2[:, iy, :] = xz_plane_fft


    arslan2 = transformed_signal2

    
    
    
    
    z = np.asmatrix(np.arange(-2.0, 2.1, 0.1)).T       # -2.0 to 2.0 in 0.1 steps
    x = np.arange(0, 2.0*math.pi, 0.1)
   # x = np.asmatrix(np.linspace(0, 2.0*math.pi, 20)).T    # 0 to 2pi in 20 steps (every 0.314)

    Lx = 4.0 * math.pi
    Lz = 2.0 * math.pi
    
    
    # run-time
    runtime = 10
    t = 0
    
    
    # m = number of modes generated.
    m = N - 2



    # Store the generated modes for each wavenumber phase speed combination
    mode = np.zeros((3.0 * m, len(kx), len(kz), len(c)), dtype=np.complex128)



    # Initialize a 3D matrix for keeping track of the first singular value
    # per mode
    singular_value = np.zeros((len(kx), len(kz), len(c)))



    # Initialize generated flow field variable
    generated_flow_field = 0



    # Scalars
    scalars = np.zeros((len(kx), spec_geo['Ny'], len(kz)))
#    scalars = np.zeros((3.0 * m, 3.0 * m))
    scalars[:, 0, :] = 1.0
    



    # Now we go through all combinations of wavenumber triplet,
    # i.e. (kx, kz, c) triplet

    for i_kx in range(0, len(kx)):          # for every stream-wise wavenumber
        for i_kz in range(0, len(kz)):          # for every span-wise wavenumber
            for i_c in range(0, len(c)):            # for every phase speed
                print('Wavenumber triplet:')
                print('kx:', kx[i_kx], '    kz:', kz[i_kz], '    c:', c[i_c])
                
                # account for kx = kz = 0 (turbulent mean)
                if kx[i_kx] == 0 and kz[i_kz] == 0:
                    continue

                
                omega = kx[i_kx] * c[i_c]
                
                
                resolvent, clencurt_quad, y = get_resolvent_operator(N, Re, 
                                                                     kx[i_kx] * (2.0 * math.pi) / Lx, 
                                                                     kz[i_kz] * (2.0 * math.pi) / Lz, 
                                                                     omega)

                # Perform SVD on the resolvent operator to get the forcing, 
                # singluar and response modes. This also gives the flow field 
                # in spectral space.
                #
                # U_spectral: resolvent modes in fourier space
                #      sigma: signular values
                # V_spectral: forcing modes in fourier space
                #
                U_spectral, sigma, V_spectral = svd(resolvent)


                # Solve linear equation to get the physical velocities in
                # streamwise direction.
                #
                #      Ax = b
                #
                #  A: clencurt_quad.T
                #  x: resolvent modes
                #  b: U_spectral
                #
                resolvent_modes = solve(clencurt_quad.T, U_spectral)
#                resolvent_modes = U_spectral
                
                
                # now we get the scalars used to re-construct the flow field
                scalars = get_scalars(u_hat, resolvent_modes, clencurt_quad)
                

                # Store the generated modes and the first singular values
                # for each wavenumber triplet.
                #
                mode[:, i_kx, i_kz, i_c] = np.squeeze(np.asarray(resolvent_modes[:, 0])) # take the first column of the resolvent modes
                singular_value[i_kx, i_kz, i_c] = sigma[0]


                # Convert the response modes to physical space to be used to 
                # generate a physical flowfield.
                #
                # ifft
                fou_field = resolvent_modes[:, 0] # take the first column of the resolvent modes matrix
                phys_flow_field   = ifft(resolvent_modes[:, 0],
                                         kx[i_kx], kz[i_kz], c[i_c], np.matrix([0]), z, t, Lx, Lz)

#                ## phys_flow_field is 3D
#                # fft
#                fourier_flow_field = fft(phys_flow_field[0, :, :],
#                                         kx[i_kx], kz[i_kz], c[i_c], np.matrix([0]), z, t, Lx, Lz)
#                
#                ## ^^^ The code works                
                
                # Generated flow field is the velocity vector U
                # which contains (u,v,w)
#                sc = np.asmatrix(scalars)
#                phy = np.asmatrix(phys_flow_field[0,:,:])
                a0 = sigma[0]
                a1 = scalars[0,:,:]
                a2 = a0 * a1
                a3 = phys_flow_field[0,:,:]
                a4 = a2 * a3
                
                generated_flow_field += (sigma[0] * scalars[0,:,:] * phys_flow_field[0,:,:])
                
                
#                u_hat = fft(generated_flow_field.real, 
#                            kx[i_kx],
#                            kz[i_kz],
#                            c[i_c],
#                            np.matrix([0]), z, t,
#                            Lx, Lz)
                
#                coeff = projection(U_spectral, clencurt_quad, u_hat)
                
#                coeff2 = get_scalars(u_hat, resolvent_modes, clencurt_quad)

                
                
    
#    
#    # fourier transform the generated flowfield.
#    # Now we are going to iterate over time.
#    for t in range(0, runtime):
#        for i_kx in range(0, len(kx)):          # for every stream-wise wavenumber
#            for i_kz in range(0, len(kz)):          # for every span-wise wavenumber
#                for i_c in range(0, len(c)):            # for every phase speed
#                    print 'Wavenumber triplet:'
#                    print 'kx:', kx[i_kx], '    kz:', kz[i_kz], '    c:', c[i_c], '    t:', t
#                    
#                    a = phys_flow_field   = ifft(mode[:, i_kx, i_kz, i_c],
#                                             kx[i_kx], kz[i_kz], c[i_c], np.matrix([0]), z, t, Lx, Lz)
#                    b = singular_value[i_kx, i_kz, i_c]
#                    
#                    c = coeff2[i_kx, i_kz, i_c]
#                    
##                    Ut += 
    
    
    return 0 



def get_resolvent_operator(N, Re, alpha, beta, omega):
    """
    INPUTS:
             N:  the number of grid points in the y-axis.
            Re:  Reynolds number
         alpha:  streamwise wavenumber
          beta:  spanwise wavenumber
         omega:  temporal frequency

    OUTPUTS:
             H:  the resolvent operator
             W:  ClenCurt matrix (Clenshaw-Curtis quadrature)
             y:  y-axis co-ordinates
    """

    # We need to get the boundary conditions before we can get the modes
    # for this we use ChebDiff.
    x, DM = ps.chebdiff(N, 2, False)

    # Second derivative matrix
    D1 = DM[0, 1:-1, 1:-1]
    # Boundary condition
    D2 = DM[1, 1:-1, 1:-1]

    # Fourth derivative matrix
    # and clamped boundary conditions
    y, D4 = ps.cheb4c(N, False)
    # returns y without endpoints

    # For the Orr-Sommerfeld equations we need to calculate the derivates
    # D = partial_dy
    # del_hat_2 = D**2.0 - K**2.0,
    #                             where K**2.0 = alpha**2.0 + beta**2.0
    # dUdy  = the  first derivative of Uo(y)
    # dU2dy = the second derivative of Uo(y)
    # f denotes the time derivative.

    # Number of modes
    m = N - 2

    I = np.eye(m)
    Z = np.zeros(shape=(m, m))
    K2 = alpha**2.0 + beta**2.0

    del_hat_2 = D2 - K2 * I
    del_hat_4 = D4 - 2.0 * D2 * K2 + K2 * K2 * I

    U = np.identity(m)  # Mean flow Uo(y), 1 at centreline
    np.fill_diagonal(U, 1 - y**2.0)
#    up.plotMatrix(U)
    dU_dy = np.identity(m)
    np.fill_diagonal(dU_dy, -2.0 * y)
#    up.plotMatrix(dU_dy)
    dU2_dy = -2.0

    # pg 60 Schmid Henningson eq3.29 and 3.30
    # -1j*M + L = 0
    OS_operator = ((1.0 / Re) * del_hat_4) + \
        (1j * alpha * dU2_dy * I) - (1j * alpha * U * del_hat_2)
    SQ_operator = ((1.0 / Re) * del_hat_2) - (1j * alpha * U)
    C_operator = -1.0j * beta * dU_dy

    # Moarreff (2013) - Model-based scaling of the streamwise energy density
    x0 = solve(del_hat_2, OS_operator)

    # State-Operator
    A = np.vstack((np.hstack((x0, Z)), np.hstack((C_operator, SQ_operator))))

    # C maps state vector to the velocity vector
    C_row0 = np.hstack(((1j / K2) * (alpha * D1), (-1j / K2) * (beta * I)))
    C_row1 = np.hstack((I, Z))
    C_row2 = np.hstack(((1j / K2) * (beta * D1), (1j / K2) * (alpha * I)))
    C = np.vstack((C_row0, C_row1, C_row2))
    C = np.asmatrix(C)

    tmp, dy = ps.clencurt(N)
    W = np.diag(np.sqrt(dy[1:-1]))
    W = np.vstack((np.hstack((W, Z, Z)),
                   np.hstack((Z, W, Z)),
                   np.hstack((Z, Z, W))))
    C = W * C

    # Adjoint of C maps the forcing vector to the state vector.
    C_adjoint_2 = pinv(C)

    # The resolvent of A
    RA = inv(-1j * omega * np.eye(m * 2) - A)

    # Transfer Function
    H = C * RA * C_adjoint_2

    return H, W, y




def fft(signal, kx, kz, c, x, z, t, Lx, Lz):
    """
    INPUTS:
     signal : Physical signal
         kx : wavenumber in x direction
         kz : wavenumber in z direction
          c : temporal frequency
          x : streamwise vector
          z : spanwise vector
          t : time instant to probe
         Lx : length of channel
         Lz : width of channel

    OUTPUT:
 fft_signal : spectral signal
    """
    
    kx *= (2.0 * math.pi) / Lx
    kz *= (2.0 * math.pi) / Lz
    
    omega = kx * c
    
    fft_signal = np.zeros((len(x), signal.shape[0], len(z)), dtype=np.complex128)
#    fft_signal = np.zeros((len(x), len(z)), dtype=np.complex128)# shape could simply be put as signal.shape

#    for i_z in range(0, len(signal)):
#        for i_x in range(0, x.shape[0]):
#            exp_component = np.exp(1j * ((kx * x[i_x]) + (kz * z[i_z]) - (omega * t)))
#            a = signal[:, i_z]
#            b = a * exp_component[0,0]
#            b = np.squeeze(b)
#            fft_signal[i_x, :, i_z] = b
            
            
    
# I should be iterating through kx, kz rather than x and z
# and also along wiht that,
# I should initialize the fourier transformed solution with the correct dimensions.
# also I think that there are u_hat, v_hat and w_hat, each should be transformed, then
# then vertically stacked so that you have a singular row of dimension 105 (3 * Ny) 
# that way when you multiply it by the resolvent modes, the dimensions match!
    
    
    
    for i_z in range(0, len(z)):
        for i_x in range(0, len(x)):
            exp_component = np.exp(-1j * ((kx * x[i_x]) + (kz * z[i_z]) - (omega * t)))
            a = signal[:]
            #b = a * exp_component[0,0]
            b = a * exp_component
            b = np.squeeze(b)
            fft_signal[i_x, :, i_z] = b
    
    fft_signal2 = fft_signal[0,:,0] # just for testing.
    
    return fft_signal2
    
    
def ifft(fft_signal, kx, kz, c, x, z, t, Lx, Lz):
    """
    INPUTS:
 fft_signal : fourier signal
         kx : wavenumber in x direction
         kz : wavenumber in z direction
          c : temporal frequency
          x : streamwise vector
          z : spanwise vector
          t : time instant to probe
         Lx : length of channel
         Lz : width of channel

    OUTPUT:
     signal : spectral signal
    """
    
    kx *= (2.0 * math.pi) / Lx
    kz *= (2.0 * math.pi) / Lz
    
    omega = kx * c
    
    signal = np.zeros((len(x), len(fft_signal), len(z)), dtype=np.complex128)
    
    
    for i_z in range(0, z.shape[0]):
        for i_x in range(0, x.shape[0]):
            exp_component = np.exp(1j * ((kx * x[i_x]) + (kz * z[i_z]) - (omega * t)))
            a = fft_signal[:, :]
            b = a * exp_component[0,0]
            b = np.squeeze(b)
            signal[i_x, :, i_z] = b
            
    
    return signal



def get_scalars(fourier_transformed_flow_field, res_modes, W):
   

    res_modes = np.asmatrix(res_modes)
    res_modes_star = res_modes.getH() # complex conjugate of the modes.
    
    W = np.asmatrix(W)

    scalars = np.zeros((fourier_transformed_flow_field.shape), dtype=np.complex128)
    
#    for i in range(0, scalars.shape[1]):
#        for z in range(0, scalars.shape[2]):
#            for x in range(0, scalars.shape[0]):
#                # column res_modes_star
#                # W
#                # u = vector of u(y), i.e. the first column of flow field post-fourier-transform
#                res_modes_star_column = np.asmatrix(res_modes_star[0, :])
#                b = np.asmatrix(fourier_transformed_flow_field[x, :, z]).T
#                
#                tmp = res_modes_star_column * W * b
#                scalars[x, i, z] = tmp[0,0]
    
    
    
    
#    
#    
#    for i_kx in range(0, len(kx)):          # for every stream-wise wavenumber
#        for i_kz in range(0, len(kz)):          # for every span-wise wavenumber
#            for i_c in range(0, len(c)):            # for every phase speed
    
    
    
    
    
    
    return scalars



# INPUTS to the main function:
N = 37
Re = 200
U_centreline = 1.0
c = np.array([1]) * U_centreline
#kx = np.array([1])
#kz = np.array([0])
kx = np.arange(-15,17)
kz = np.arange(0,17)

# Execute main function:
main_resolvent_analysis(N, Re, kx, kz, c)
