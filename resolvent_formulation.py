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
        gen_ff:  the flowField that is generated from the resolvent 
                 formulation.
                      
    """

    # TIDY THIS FUNCTION
        
    # Geometric considerations:
    channel_half_height = 1.0
    u_hat = 0
    
    if modesOnly:
        # Can use independent geometrical values
        Nx = 80
        Nz = 70
        
        x = np.linspace(0.0, 2.0*math.pi, Nx)
        z = np.linspace(0.0, 2.0*math.pi, Nz)
        
        Ny = N - 2
        Nd = 3
        
        Lx = 2.0*math.pi
        Lz = 2.0*math.pi
        
        y = np.linspace(-channel_half_height, channel_half_height, Ny)

        dP_dx = 0.005 # I use this value so that the centreline velocity is 1
        laminar_baseflow = np.zeros((Nx, Ny, Nz))
        for iy in range(0, Ny):
            laminar_baseflow[:, iy, :] = (dP_dx * Re * 0.5) * (1.0 - y[iy]**2)
        
        yz_component = np.zeros((Nx, Ny, Nz))
        
        laminar_baseflow = np.concatenate((laminar_baseflow, yz_component, yz_component), axis=1)
        
#        f = np.zeros((Nx, Nd*Ny, Nz))

        
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
            
            
            
    # according to Gibson's channel flow code, kinematic_visc = 1/Re
    Ubulk = (1.0) / (channel_half_height * 2.0)
    
    # Ubulk is actually calculated as follows:
    # take the integral from -1 to 1 of laminar flow:
    # i.e.
    # integrate 1-y^2 to get y - y^3/3
    # putting in the limits [-1,1]
    # we get
    # Ubulk = 4/3

    t = 1  # run-time
    m = N - 2 # number of modes
    
    # Initialize an object array to store the generated modes for each
    # wavenumber phase speed combination
    mode = np.zeros((3.0*m, len(kx), len(kz)), dtype=np.complex128)
    
    # Initialize a 3D matrix for keeping track of the first 
    # singular value per mode
    singular_value = np.zeros((len(kx), len(kz)))
    
    
    # Scalars
    if not modesOnly:
        if data['flowField']['is_spectral']:
            scalars = np.zeros((len(x), 3.0*m, len(z)), dtype=np.complex128)
            
    elif modesOnly:
        scalars = np.ones((len(kx), 3.0*m, len(kz)), dtype=np.complex128)
    
        



    U = 0
    laminarOnly = True
    if laminarOnly:
        lam = 'laminar'
    else:
        lam = ''
    
    
    # Amplitudes
    # A = sigma * khai
    amplitude = 1.0e-5
    
    
    for ikx in range(0, len(kx)):
        for ikz in range(0, len(kz)):
            # account for kx = kz = 0 (turbulent mean)
            if kx[ikx] == 0 or kz[ikz] == 0:
                continue
            
            print 'kx:',kx[ikx],'    kz:',kz[ikz], '    A:',amplitude
            
            
            # Calculate the temporal frequency
            omega = kx[ikx] * c
            
            
            # Get the state vectors so that you can compute the resolvent operator
            # and hence the transfer function.
            # 
            # The calculations given below (and variable names) are outlined in:
            #     Moarref, Model-based scaling of the streamwise energy 
            #     density in high-Reynolds number turbulent channels, 2013.
            C, C_adj, A, w, y_cheb = get_state_vectors(N, Re, Lx, Lz, kx[ikx], kz[ikz])
            
            
            # The resolvent of state operator A
            R_A = inv(-1j*omega*np.eye(m*2) - A)
            
            
            # Transfer function
            H = C*R_A*C_adj
            
            
            # Perform SVD on the resolvent operator (R_A) to get the forcing, 
            # singluar and response modes. This also gives a flow field 
            # in spectral space.
            # 
            #   U_spectral: response modes in fourier space (resolvent modes)
            #        sigma: signular values
            #   V_spectral: forcing modes in fourier space
            U_spectral, sigma, V_spectral = svd(H)
            
    
            # Solve linear equation to get the physical velocities in
            # streamwise direction.
            #     Ax = b
            #      A: w.T
            #      x: resolvent modes
            #      b: U_spectral
            resolvent_modes = solve(w.T, U_spectral) # Make sure you disassociate the modes from the clencurt quadrature.
            
            
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            # Re-check the following block .. .. ..
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            
#            if not modesOnly:
#                # now we get the scalars used to re-construct the flow field
#                tmp_scalar = get_scalars(u_hat[ikx, :, ikz], resolvent_modes, w)
#                scalars[ikx, :, ikz] = tmp_scalar
#                sc_tmp = np.delete(scalars, (0), axis=2)
#                mirroredSc = sc_tmp[:,:,::-1]
#                doublesSc = np.concatenate((mirroredSc, sc_tmp), axis=2)
#                scaShape = doublesSc.shape
            
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
            
            
            # Store the generated modes and the first singular values
            # for each wavenumber triplet.
            # take the first column of the resolvent modes
            mode[:, ikx, ikz] = np.squeeze(np.asarray(resolvent_modes[:, 0])) 
            singular_value[ikx, ikz] = sigma[0]
            
            
            # Calculate u_tilde
            if modesOnly:
                u_tilde = amplitude * resolvent_modes
            else:
                u_tilder = sigma.T * scalars * resolvent_modes
            
            
            # Convert the resolvent modes to physical space
            physical_ff = np.zeros((len(x), 3*m, len(z)), dtype=np.complex128)


            # Rank of resolvent modes, this can be lowered to use less of the 
            # modes to calculate the physical flow field.
            rank = 13
            rank = min(rank, 3*m)
            for iy in range(0, rank):
                physical_ff += ifft(u_tilde[:,iy], kx[ikx], kz[ikz], omega, x, z, t, Lx, Lz)
                
            
            # Generated flow field is the velocity vector U
            # which contains (u,v,w)
            # U += (sigma[0] * scalars[0,:,:] * physical_ff[0,:,:])
            
          
            # Generated flow field
            if modesOnly:
                if laminarOnly:
                    U = laminar_baseflow
                else:
                    U += physical_ff
                
            else:
                U += physical_ff
        
        
        fmt = '{0:<20} {1:<20}'
        string_kx = str(kx[ikx])
        string_kz = str(kz[ikz])
        string_c = format(c, '.4f')
        string_A = str(amplitude)
        



    U = U.real
#    U = f
    U_u = U[:,   0:m  , :]
    U_v = U[:,   m:2*m, :]
    U_w = U[:, 2*m:3*m, :]
    
    
    # Output the flow field as an ASCII file for channelflow to read in.
    
    # Here I need ot construct my geometry file and pack the dictionary such that
    # when it's opened in utils, the writing is as easy as channelflow writing 
    # to file.
    # I store the 4D matrix as follows:
    # u[i, nx, ny, nz]
    # so:
    U = np.zeros((Nd, Nx, Ny, Nz))
    
    for i in range(0, Nd):
        for nx in range(0, Nx):
            for ny in range(0, Ny):
                for nz in range(0, Nz):
                    if i == 0: # u direction
                        U[i, nx, ny, nz] = U_u[nx, ny, nz]
                    elif i == 1: # v direction
                        U[i, nx, ny, nz] = U_v[nx, ny, nz]
                    
                    elif i == 2: # w direction
                        U[i, nx, ny, nz] = U_w[nx, ny, nz]
                        
                    
                    
    
    gen_ff = {}
    gen_ff['resolvent_flowField'] = U
    gen_ff['U'] = U_u
    gen_ff['V'] = U_v
    gen_ff['W'] = U_w
    gen_ff['X'] = x
    gen_ff['Y'] = y
    gen_ff['Z'] = z

    gen_ff['Nx'] = Nx
    gen_ff['Ny'] = m
    gen_ff['Nz'] = Nz
    gen_ff['Nd'] = Nd

    gen_ff['Lx'] = Lx
    gen_ff['Lz'] = Lz
    
    gen_ff['kx'] = string_kx
    gen_ff['kz'] = string_kz
    gen_ff['c'] = string_c
    gen_ff['A'] = string_A
    gen_ff['lam'] = lam
    
    
    return gen_ff



def chebyshev_points(N, a, b):
    """
    INPUTS:
    N:  number of grid points
    a:  upper bound of domain
    b:  lower bound of domain
    
    OUTPUTS:
    y:  Chebyshev points
    """

    y_chebheb = np.arange(N, dtype=np.float)
    pi_N = math.pi / (N - 1.0)
    radius = (b - a) / 2.0
    center = (b + a) / 2.0
    for j in range(0,N):
        tmp0=j*pi_N
        tmp1=math.cos(j*pi_N)
        y_chebheb[j] = center + radius*tmp1
    
    return y_chebheb



def get_state_vectors(N, Re, Lx, Lz, kx, kz):
    """
    We are calculating the state vectors in this function. The methodology
    followed here is given in the following reference in the "Formulation" 
    section, 
        2. Low-rank approximation to channel flow,
        
        (Moarref, Model-based scaling of the streamwise energy 
        density in high-Reynolds number turbulent channels, 2013)
    
    
    INPUTS:
             N:  number of grid points in the y-axis.
            Re:  Reynolds number
            Lx:  length of solution
            Lz:  width of solution
            kx:  streamwise wavenumber
            kz:  spanwise wavenumber
            
     
    OUTPUTS:
             C:  this operator maps the state vector onto the velocity vector
         C_adj:  adjoint of C, maps the forcing vector to the state vector
             A:  state operator
             W:  ClenCurt matrix (Clenshaw-Curtis quadrature)
             y:  grid-points in the y-axis
    """
    
    
    kx = kx * (2.0 * math.pi) / Lx
    kz = kz * (2.0 * math.pi) / Lz
    
    
    # Calculate the differentiation matrix, DM, for the resolution in the 
    # y-axis, N. 
    # y_cheb are the interpolated y co-ordinates, i.e. Chebyshev interior points.
    y_cheb, DM = ps.chebdiff(N, 2)
    
    
    # First derivative matrix
    D1 = DM[0, 1:-1, 1:-1]
    
    # Second derivative matrix
    D2 = DM[1, 1:-1, 1:-1]


    # Fourth derivative matrix and clamped boundary conditions
    tmp, D4 = ps.cheb4c(N, False)
    # tmp is the same as y_cheb without endpoints, i.e. no [-1,1]
    
    
    # For the Orr-Sommerfeld equations we need to calculate the derivates
    #            D:  partial_dy
    #    del_hat_2:  D**2.0 - K**2.0         (where K**2.0 = kx**2.0 + kz**2.0)
    #         dUdy:  first derivative of Uo(y)
    #        dU2dy:  second derivative of Uo(y)
    #            f:  time derivative
    
    
    # Number of modes
    m = N - 2 # calculate without the endpoints, otherwise m = N
    I = np.identity(m)
    Z = np.zeros(shape=(m, m))
    K2 = kx**2.0 + kz**2.0
    del_hat_2 = D2 - K2*I
    del_hat_4 = D4 - 2.0*D2*K2 + K2*K2*I
    
    
    # Laminar Base flow 
    U = np.identity(m) 
    np.fill_diagonal(U, 1 - y_cheb**2.0) # 1 at centreline
    
    dU_dy  = np.identity(m)
    np.fill_diagonal(dU_dy, -2.0*y_cheb)
    
    dU2_dy = -2.0


    # pg 60 Schmid Henningson eq3.29 and 3.30
    # -1j*M + L = 0
    SQ_operator = ((1.0/Re) * del_hat_2) - (1j * kx * U)
    C_operator  = -1.0j * kz * dU_dy
    
    OS_operator = ((1.0/Re) * del_hat_4) + (1j * kx * dU2_dy * I) - (1j * kx * U * del_hat_2)
    x0 = solve(del_hat_2, OS_operator)
    
    # Equation 2.7
    # (Moarref, Model-based scaling of the streamwise energy density in 
    # high-Reynolds number turbulent channels, 2013)
    #
    # State operator
    # A = | x0   Z  |
    #     |  C   SQ |
    A = np.vstack((np.hstack((x0,Z)), np.hstack((C_operator, SQ_operator))))

 
    # C maps state vector to the velocity vector
    C = np.vstack((np.hstack(((1j/K2) * (kx*D1), (-1j/K2) * (kz*I))), 
                   np.hstack((                I,                Z)), 
                   np.hstack(((1j/K2) * (kz*D1), ( 1j/K2) * (kx*I)))))
    
    C = np.asmatrix(C)
    
    
    tmp, dy = ps.clencurt(N);
    W = np.diag(np.sqrt(dy[1:-1]))
    W = np.vstack((np.hstack((W,Z,Z)),
                   np.hstack((Z,W,Z)),
                   np.hstack((Z,Z,W))))
    C *= W
    
    
    # Adjoint of C maps the forcing vector to the state vector. 
    C_adj = pinv(C)
    
    
    return C, C_adj, A, W, y_cheb
        
    
    
def ifft(fft_signal, kx, kz, omega, x, z, t, Lx, Lz):
    """
    INPUTS:
 fft_signal : fourier signal
         kx : wavenumber in x direction
         kz : wavenumber in z direction
      omega : temporal frequency
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
    

def get_scalars(u_hat, resolvent_modes, w):

    resolvent_modes = np.asmatrix(resolvent_modes)
    resolvent_modes_star = resolvent_modes.getH() # complex conjugate of the modes.
    
    w = np.asmatrix(w)

    scalars = np.zeros((u_hat.shape), dtype=np.complex128) # (3*Ny)
    
    u_hat = np.asmatrix(u_hat)

    
    
    for i in range(0, u_hat.shape[1]): # from 0 to 3*Ny - 1
        resolvent_modes_col = np.asmatrix(resolvent_modes_star[:,i]) # column of resolvent modes
        # Dimensions:
        #    1     =  1xN *     NxN     *       Nx1
        tmp = u_hat * w * resolvent_modes_col
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