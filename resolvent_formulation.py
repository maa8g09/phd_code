"""
RESOLVENT FORMULATION

Resolvent formulation and projection code. This is the main routine that calls 
all of the functions and files that perform the whole routine of generating 
resolvent modes and then projecting them onto a channel flow solution.


Author details:
    Muhammad Arslan Ahmed
    maa8g09@soton.ac.uk
    
    Aerodynamics and Flight Mechanics Research Group
    Faculty of Engineering and the Environment
    University of Southampton
"""

import tools_pseudospectral as ps
import utils_plots as up
import math
import numpy as np

from scipy.linalg import inv
from scipy.linalg import solve
from scipy.linalg import pinv
from scipy.linalg import svd


def main_resolvent_analysis(N, Re, kx, kz, c, modesOnly, data, fourdarray):
    """
    The full resolvent formulation is given here, from the generation of the 
    modes to the projection of a channelflow solution onto the modes.

    
    INPUTS:
             N:  resolution in y axis
            Re:  Reynolds number
            kx:  vector of streamwise wavenumbers
            kz:  vector of spanwise wavenumbers
             c:  phase speed
     modesOnly:  boolean that controls whether or not to project modes onto a 
                 solution
          data:  a dictionary with all of the necassary flow field solutions in
                 it classed by physical and spectral:
                     The physical flow field is stored as a 4D array: 
                      - (i, nx, ny, nz)
                     The spectral flow field is stored as a 4D array:
                      - (i, kx, ny, kz)
    fourdarray:  a 1D array with 4 variables to use to plot projected solution

    
    OUTPUTS:

    """
    
    
    u_hat = 0
    
    if modesOnly:
        # Can use independent geometrical values
        x = np.arange(0.0, 2.0*math.pi, 0.1)
        z = np.arange(-2.0, 2.1, 0.1)
        Nx = len(kx)
        Nz = 2*len(kz) - 1
        Lx = 4.0*math.pi #2 periods
        Lz = 2.0*math.pi #1 period
        
        
    else:
        # Or can use channelflow geometry dimensions to generate resolvent modes
        if data['flowField']['is_spectral'] == True:
            u_hat = data['flowField']['spectral']
            u_hat = np.concatenate((u_hat[0,:,:,:],
                                    u_hat[1,:,:,:],
                                    u_hat[2,:,:,:]), 
                                    axis=1)
            
            Nx = data['geometry']['spectral']['kx']
            N = data['geometry']['spectral']['Ny'] + 2
            Nz = data['geometry']['spectral']['kz']*2 - 2
            Lx = data['geometry']['physical']['Lx']
            Lz = data['geometry']['physical']['Lz']
            x = data['geometry']['physical']['x']
            z = data['geometry']['physical']['z']
            kx = np.arange(-0.5*Nx + 1, 0.5*Nx + 1)
            kz = np.arange(data['geometry']['spectral']['kz'])
            
            
        elif data['flowField']['is_physical'] == True:
            x = data['geometry']['physical']['x']
            z = data['geometry']['physical']['z']
            Nx = data['geometry']['physical']['Nx']
            Nz = data['geometry']['physical']['Nz']
            N = data['geometry']['physical']['Ny'] + 2
            Lx = data['geometry']['physical']['Lx']
            Lz = data['geometry']['physical']['Lz']
            kx = np.arange(Nx)
            kz = np.arange(-0.5*Nz, 0.5*Nz + 1)
            
            # Fourier transform the channelflow solution
            u_hat = fft_flowField(data)
            
        


    t = 1  # run-time
    
    delta_c = np.abs(c[0] - c[1])
    m = N - 2 # number of modes
    
    
    # Initialize an object array to store the generated modes for each
    # wavenumber phase speed combination
    mode = np.zeros((3.0*m, len(kx), len(kz), len(c)), dtype=np.complex128)
    
    # Initialize a 3D matrix for keeping track of the first 
    # singular value per mode
    singular_value = np.zeros((len(kx), len(kz), len(c)))
    
    U = 0
    
    
    # Scalars
    if not modesOnly:
        if data['flowField']['is_spectral']:
            scalars = np.zeros((len(kx), 3.0*m, len(kz)), dtype=np.complex128)
            
    elif modesOnly:
        scalars = np.ones((len(x), 3.0*m, len(z)), dtype=np.complex128)
    
    
    for i_kx in range(0, len(kx)):
        for i_kz in range(0, len(kz)):
            
            # account for kx = kz = 0 (turbulent mean)
            if kx[i_kx] == 0 or kz[i_kz] == 0:
                continue
            
            C, C_adjoint_2, A, clencurt_quad, y = get_state_vectors(N, Re,
                                                                         Lx, Lz,
                                                                         kx[i_kx],
                                                                         kz[i_kz])
                                                                 
            for i_c in range(0, len(c)):
                print 'kx:',kx[i_kx],'    kz:',kz[i_kz],'    c:',c[i_c]

                
                # Phase speed
                omega = kx[i_kx] * c[i_c]


                # The resolvent of A
                resolvent_operator = inv(-1j*omega*np.eye(m*2) - A)


                # Transfer Function
                transfer_function = C * resolvent_operator * C_adjoint_2
                
                
                # Perform SVD on the resolvent operator to get the forcing, 
                # singluar and response modes. This also gives a flow field 
                # in spectral space.
                #
                # U_spectral: resolvent modes in fourier space
                #      sigma: signular values
                # V_spectral: forcing modes in fourier space
                #
                U_spectral, sigma, V_spectral = svd(transfer_function)
                

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
                
                
                if not modesOnly:
                    # now we get the scalars used to re-construct the flow field
                    tmp_scalar = get_scalars(u_hat[i_kx, :, i_kz], resolvent_modes, clencurt_quad)
                    scalars[i_kx, :, i_kz] = tmp_scalar
                    sc_tmp = np.delete(scalars, (0), axis=2)
                    mirroredSc = sc_tmp[:,:,::-1]
                    doublesSc = np.concatenate((mirroredSc, sc_tmp), axis=2)
                    scaShape = doublesSc.shape


                # Store the generated modes and the first singular values
                # for each wavenumber triplet.
                mode[:, i_kx, i_kz, i_c] = np.squeeze(np.asarray(resolvent_modes[:, 0])) # take the first column of the resolvent modes
                singular_value[i_kx, i_kz, i_c] = sigma[0]


                # Convert the response modes to physical space to be used to 
                # generate a physical flowfield.
                #
                # ifft
                fou_field = resolvent_modes[:, 0] # take the first column of the resolvent modes matrix
                phys_flow_field_0   = ifft(resolvent_modes[:, 0],
                                         kx[i_kx], kz[i_kz], c[i_c], np.array([0]), z, t, Lx, Lz, omega)
                                         
                phys_flow_field   = ifft(resolvent_modes[:, 0],
                                         kx[i_kx], kz[i_kz], c[i_c],x, z, t, Lx, Lz, omega)
                
                
                # Generated flow field is the velocity vector U
                # which contains (u,v,w)
#                U += (sigma[0] * scalars[0,:,:] * phys_flow_field[0,:,:])
              
              
                # Generated flow field
                if modesOnly:
                    U += (scalars * sigma[0] * delta_c * phys_flow_field)
                else:
                    U += (doublesSc * sigma[0] * delta_c * phys_flow_field)

                
    generated_flowField = {}
    generated_flowField['resolvent_flowField'] = U.real
    generated_flowField['U'] = U[:,   0:m   ,:].real
    generated_flowField['V'] = U[:,   m:2*m ,:].real
    generated_flowField['W'] = U[:, 2*m:3*m ,:].real
    generated_flowField['X'] = x
    generated_flowField['Y'] = y
    generated_flowField['Z'] = z
    generated_flowField['Nx'] = Nx
    generated_flowField['Ny'] = m
    generated_flowField['Nz'] = Nz
        
    data['resolvent_flowField'] = U
    
    return generated_flowField


def get_state_vectors(N, Re, Lx, Lz, alpha, beta):
    """
    We are calculating the state vectors in this function. The methodology
    followed here is given in the following reference in the "Formulation" 
    section (2. Low-rank approximation to channel flow), Moarreff, 2013, 
    Model-based scaling of the streamwise energy density.
    
    
    INPUTS:
             N:  the number of grid points in the y-axis.
            Re:  Reynolds number
            Lx:  length of solution box
            Lz:  width of solution box
         alpha:  streamwise wavenumber
          beta:  spanwise wavenumber
            
     
    OUTPUTS:
             C:  this operator maps the state vector onto the velocity vector
   C_adjoint_2:  adjoint of C, maps the forcing vector to the state vector
             A:  state operator
             W:  ClenCurt matrix (Clenshaw-Curtis quadrature)
             y:  grid-points in the y-axis
    """
    
    
    alpha = alpha * (2.0 * math.pi) / Lx
    beta = beta * (2.0 * math.pi) / Lz
    
    
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
    #          D = partial_dy
    #  del_hat_2 = D**2.0 - K**2.0, where K**2.0 = alpha**2.0 + beta**2.0
    #       dUdy = first derivative of Uo(y)
    #      dU2dy = second derivative of Uo(y)
    #          f = time derivative.
    
    
    # Number of modes
    m = N - 2 
    
    I = np.eye(m)
    Z = np.zeros(shape=(m, m))
    K2 = alpha**2.0 + beta**2.0

    del_hat_2 = D2 - K2*I
    del_hat_4 = D4 - 2.0*D2*K2 + K2*K2*I
    

    U = np.identity(m) # Mean flow Uo(y), 1 at centreline
    np.fill_diagonal(U, 1 - y**2.0)
#    up.plotMatrix(U)
    dU_dy  = np.identity(m)
    np.fill_diagonal(dU_dy, -2.0*y)
#    up.plotMatrix(dU_dy)
    dU2_dy = -2.0


    # pg 60 Schmid Henningson eq3.29 and 3.30
    # -1j*M + L = 0
    OS_operator = ((1.0 / Re) * del_hat_4) + (1j * alpha * dU2_dy * I) - (1j * alpha * U * del_hat_2)
    SQ_operator = ((1.0 / Re) * del_hat_2) - (1j * alpha * U)
    C_operator  = -1.0j * beta * dU_dy


    # Moarreff (2013) - Model-based scaling of the streamwise energy density
    x0 = solve(del_hat_2, OS_operator)
        
    # State-Operator
    A = np.vstack((np.hstack((x0,Z)), np.hstack((C_operator, SQ_operator))))

 
    # C maps state vector to the velocity vector
    C = np.vstack((np.hstack((( 1j / K2) * (alpha * D1), (-1j / K2) *  (beta * I))), 
                   np.hstack((                        I,                       Z)), 
                   np.hstack(( ( 1j / K2) * (beta * D1), ( 1j / K2) * (alpha * I)))))
    
    C = np.asmatrix(C)
    
    tmp, dy = ps.clencurt(N);
    W = np.diag(np.sqrt(dy[1:-1]))
    W = np.vstack((np.hstack((W,Z,Z)),
                   np.hstack((Z,W,Z)),
                   np.hstack((Z,Z,W))))
    C = W*C;
    
    
    # Adjoint of C maps the forcing vector to the state vector. 
    C_adjoint_2 = pinv(C)
    
    
    return C, C_adjoint_2, A, W, y
    
        
    
    
def ifft(fft_signal, kx, kz, c, x, z, t, Lx, Lz, omega):
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
    
    
    signal = np.zeros((len(x), len(fft_signal), len(z)), dtype=np.complex128)
    
    
    for i_z in range(0, z.shape[0]):
        for i_x in range(0, x.shape[0]):
            exp_component = np.exp(1j * ((kx * x[i_x]) + (kz * z[i_z]) - (omega * t)))
            a = fft_signal[:]
            b = a * exp_component
            b = np.squeeze(b)
            signal[i_x, :, i_z] = b
            
    
    return signal


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
    
    delta_x = np.abs(x[0] - x[1])
    delta_z = np.abs(z[0] - z[1])
    
#    fft_signal = np.zeros((len(x), len(z)), dtype=np.complex128)# shape could simply be put as signal.shape

#    for i_z in range(0, len(signal)):
#        for i_x in range(0, x.shape[0]):
#            exp_component = np.exp(1j * ((kx * x[i_x]) + (kz * z[i_z]) - (omega * t)))
#            a = signal[:, i_z]
#            b = a * exp_component[0,0]
#            b = np.squeeze(b)
#            fft_signal[i_x, :, i_z] = b

    for i_z in range(0, len(z)):
        for i_x in range(0, len(x)):
            exp_component = np.exp(-1j * ((kx * x[i_x]) + (kz * z[i_z]) - (omega * t)))
            a = signal[:, :]
            b = a * exp_component * delta_x * delta_z
            b = np.squeeze(b)
            fft_signal[i_x, :, i_z] = b
    
    
    reciprocal = (1.0 / (Lx * Lz))
    fft_signal *= reciprocal
    
    fft_signal2 = fft_signal[0,:,0] # just for testing.
    
    return fft_signal2

def fft_flowField(data):
    
    return data
    

def get_scalars(u_hat, resolvent_modes, clencurt_quad):

    resolvent_modes = np.asmatrix(resolvent_modes)
    resolvent_modes_star = resolvent_modes.getH() # complex conjugate of the modes.
    
    clencurt_quad = np.asmatrix(clencurt_quad)

    scalars = np.zeros((u_hat.shape), dtype=np.complex128) # (3*Ny)
    
    u_hat = np.asmatrix(u_hat)

    
    
    for i in range(0, u_hat.shape[1]): # from 0 to 3*Ny - 1
        resolvent_modes_col = np.asmatrix(resolvent_modes_star[:,i]) # column of resolvent modes
        # Dimensions:
        #    1     =  1xN *     NxN     *       Nx1
        tmp = u_hat * clencurt_quad * resolvent_modes_col
        scalars[i] = tmp[0,0]
    
    a = scalars
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


    
    return scalars