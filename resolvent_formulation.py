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
import utils as ut
import utils_plots as up
import numpy as np


from scipy.interpolate import interp1d
from colorama import Fore, Back, Style
from math import pi
from numpy.linalg import inv
from numpy.linalg import solve
from numpy.linalg import pinv
from numpy.linalg import svd


import matplotlib.pyplot as plt




def main_resolvent_analysis(N, Re, kx, kz, c, amplitudes, modesOnly, data, fourdarray):
    """
    The full resolvent formulation is given here, from the generation of the 
    modes to the projection of a channelflow solution onto the modes.
    
    
    INPUTS:
             N:  resolution in y axis
            Re:  Reynolds number
            kx:  vector of streamwise wavenumbers
            kz:  vector of spanwise wavenumbers
             c:  phase speed
             A:  amplitude
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
    
    
    # Geometric considerations:
    channel_half_height = 1.0
    u_hat = 0
    
    if modesOnly:
        # Can use independent geometrical values
        Nx = 40
        Nz = 30
        
        alpha = kx[0]
        beta = kz[0]
        
        
        lamda_x = 2.0*pi / alpha
        lamda_z = 2.0*pi / beta
        
        
        Lx = lamda_x
        Lz = lamda_z
        
        
        x = np.linspace(0.0, Lx, Nx)
        z = np.linspace(-Lz/2.0, Lz/2.0, Nz)
        
        
        Mx = alpha * np.arange(Nx + 1)
        Mz = beta  * np.arange(Nz + 1)
        
        
        Ny = N - 2
        Nd = 3
        
        
        y = np.linspace(-channel_half_height, channel_half_height, Ny)
        
        
        
    else:
        # Or can use channelflow geometry dimensions to generate resolvent modes
        if data['flowField']['is_spectral'] == True:
            u_hat = data['flowField']['spectral']
            u_hat = np.concatenate((u_hat[0,:,:,:],
                                    u_hat[1,:,:,:],
                                    u_hat[2,:,:,:]), 
                                    axis=1)
            
            Nx = data['geometry']['spectral']['kx']
            Ny = data['geometry']['spectral']['Ny']
            N = Ny + 2
            Nz = data['geometry']['spectral']['kz']*2 - 2
            Nd = data['geometry']['spectral']['Nd']
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
            Ny = data['geometry']['spectral']['Ny']
            N = Ny + 2
            Nd = data['geometry']['physical']['Nd']
            Lx = data['geometry']['physical']['Lx']
            Lz = data['geometry']['physical']['Lz']
            kx = np.arange(Nx)
            kz = np.arange(-0.5*Nz, 0.5*Nz + 1)
            
            # Fourier transform the channelflow solution
            u_hat = fft_flowField(data)
            


    # Ubulk is calculated as follows:
    # take the integral from -1 to 1 of laminar flow:
    # i.e.
    # integrate 1-y^2 to get y - y^3/3
    # putting in the limits [-1,1]
    # we get
    # Ubulk = 4/3


    t = 0  # run-time
    m = N - 2 # number of modes
    
    # Initialize an object array to store the generated modes for each
    # wavenumber phase speed combination
#    mode = np.zeros((kx.shape[0], 3.0*m, len(kx), len(kz)), dtype=np.complex128)
    
    # Initialize a 3D matrix for keeping track of the first 
    # singular value per mode
#    singular_value = np.zeros((kx.shape[0], len(kx), len(kz)))
    
    

    if not modesOnly:
        if data['flowField']['is_physical'] == True:
            alpha = data['geometry']['physical']['alpha']
            beta = data['geometry']['physical']['gamma']
    else:
        alpha = 1.0
        beta  = 1.0
    
    
    physical_ff = np.zeros((len(x), 3*m, len(z)), dtype=np.complex128)
    generated_ff = np.zeros((len(x), 3*m, len(z)), dtype=np.complex128)

    # Amplitudes
    # A = S * khi

#    for mode_num in range(0, kx.shape[0]):
#        # the number of modes we are generating to create a packet
#        # So we will loop through each set of modes...
#        kx_mode_num = kx[mode_num]
#        kz_mode_num = kz[mode_num]
#        
#        for ikx in range(0, len(kx_mode_num)):
#            for ikz in range(0, len(kz_mode_num)):
#                # account for kx = kz = 0 (turbulent mean)
#                if kx_mode_num[ikx] == 0 or kz_mode_num[ikz] == 0:
#                    continue
#                
#                print('kx:',kx_mode_num[ikx],'    kz:',kz_mode_num[ikz], '    A:',amplitude[mode_num])
#                
#                
#                # Calculate the temporal frequency
##                omega = 2.0*pi*kx_mode_num[ikx]/Lx * c
#                
#                
#                # Get the state vectors so that you can compute the resolvent operator
#                # and hence the transfer function.
#                # 
#                # The calculations given below (and variable names) are outlined in:
#                #     Moarref, Model-based scaling of the streamwise energy 
#                #     density in high-Reynolds number turbulent channels, 2013.
#                iter_kx=kx_mode_num[ikx]*alpha
#                iter_kz=kz_mode_num[ikz]*beta
#                
#                omega = kx_mode_num[ikx] * c
#                omega = iter_kx * c
#                
#    
#                C, C_adj, A, w, y_cheb = get_state_vectors(N, Re, Lx, Lz, 
#                                                           iter_kx,
#                                                           iter_kz)
#                
#                I = np.eye(A.shape[0])
#                # The resolvent of state operator A
#                L = -1j*omega*I - A
#                RA = inv(L)
#                
#                
#                # Transfer function
#                H = C*RA*C_adj
#                
#                # Transfer function testing output...
##                Hreal = H.real
##                Himag = H.imag
#                
##                fileName='/Hreal.txt'
##                file = open('/home/arslan/Documents' + fileName, "w")
##                for i in range(0, Hreal.shape[0]):
##                    tmp = Hreal[i,:]#take each row
##                    for q in range(0, tmp.shape[1]):#write each element
##                        tmp0 = str(tmp[0,q])
##                        file.write(tmp0 + " ")
##                    file.write("\n")
##                file.close()
##                
##                fileName='/Himag.txt'
##                file = open('/home/arslan/Documents' + fileName, "w")
##                for i in range(0, Hreal.shape[1]):
##                    tmp = Himag[i,:]#take each row
##                    for q in range(0, tmp.shape[1]):#write each element
##                        tmp0 = str(tmp[0,q])
##                        file.write(tmp0 + " ")
##                    file.write("\n")
##                file.close()
#                
#                
#                
#                # Perform SVD on the resolvent operator (R_A) to get the forcing, 
#                # singluar and response modes. This also gives a flow field 
#                # in spectral space.
#                # 
#                #   U_spectral: response modes in fourier space (resolvent modes)
#                #            S: signular values
#                #   V_spectral: forcing modes in fourier space
#                #
#                # In matlab U S V  = svd(H)
#                # In python U S Vh = svd(H), where V = Vh.T
#                #
#                U_spectral, S, V_spectral = svd(H)
#                
#                if np.linalg.norm(np.dot( np.dot(U_spectral, np.diag(S)), V_spectral) - H) >= 1e-10:
#                    nrm = str(np.linalg.norm(np.dot( np.dot(U_spectral, np.diag(S)), V_spectral) - H))
#                    err = 'Something went wrong with the SVD, norm is '+nrm
#                    ut.error(err)
#                    
#                
#                # Solve linear equation to recover the velocities
#                #     Ax = b
#                #      A: w.T (Clencurt vector has been converted to a diagonal matrix)
#                #      x: resolvent modes
#                #      b: U_spectral
#                resolvent_modes = solve(w.T, U_spectral)
#                
#                
#                # Store the generated modes and the first singular values
#                # for each wavenumber triplet.
#                # take the first column of the resolvent modes
##                mode[mode_num, :, ikx, ikz] = np.squeeze(np.asarray(resolvent_modes[:, 0])) 
##                singular_value[mode_num, ikx, ikz] = S[0]
#                
#                
#                # Calculate u_tilde
#                if modesOnly:
#                    tmp= w * resolvent_modes[:, 0]
#                    u_tilde = amplitude[mode_num] * resolvent_modes[:, 0] # np.diag(tmp)
#    #                u_tilde = amplitude * resolvent_modes[:, 0]
#                    
#                else:
#                    # To get the weighting we can use the scalars that I can 
#                    # get from Gibson's solutions.
#                    # Otherwise I use a constant coefficient as the weighting
#                    # factor.
#                    scalars = get_scalars(u_hat[ikx, :, ikz], resolvent_modes[:, 0], w)
#                    amplitude = np.asmatrix(S[0]) * np.asmatrix(scalars).T
#                    u_tilde = amplitude[0,0] * resolvent_modes
#                
#                
#                u_tilde += np.conjugate(u_tilde)
#                
#                # add the complex conjugate of u_tilde to u_tilde before you inverse fourier it.
#                # alternatively we could add the complex conjugate after the physical flowfield has been 
#                # generated.
#                
#                u_tilde = np.asmatrix(u_tilde)
#                # Convert the resolvent modes to physical space
#                physical_ff = np.zeros((len(x), 3*m, len(z)), dtype=np.complex128)
#    
#    
#                # Number of resolvent modes to use
#                num_modes = 1
#                num_modes = min(num_modes, 3*m)
#                
#                for iy in range(0, num_modes):
#                    physical_ff += ifft(u_tilde[:,iy], iter_kx, iter_kz, c, x, z, t, Lx, Lz, m)
##                    physical_ff += ifft(u_tilde[:,iy], kx_mode_num[ikx], kz_mode_num[ikz], c, x, z, t, Lx, Lz, m)
#                    
#                
#                # Generated flow field is the velocity vector U
#                # which contains (u,v,w)
#                # U += (S[0] * scalars[0,:,:] * physical_ff[0,:,:])
#                #
#                # Generated flow field is just the physical_ff
#                generated_ff += physical_ff.real
#    
#            
#            string_kx = str(kx_mode_num[ikx])
#            string_kz = str(kz_mode_num[ikz])
#            string_c = format(c, '.4f')
#            string_A = str(amplitude[mode_num])
#        
#
#
#
#
#    U = np.zeros((Nd, Nx, Ny, Nz))
#    U_u = generated_ff[:,   0:m  , :]
#    U_v = generated_ff[:,   m:2*m, :]
#    U_w = generated_ff[:, 2*m:3*m, :]












    lambda_x = 2.0*pi / kx[0]
    lambda_z = 2.0*pi / kz[0]
    
    Lx = lambda_x
    Lz = lambda_z
    
    
    

    # Number of grid points
    Nx = 32 # even
    Nz = 33 # even

    # Stationary nodes along each axis
    # X axis
    Mx = Nx
    Mx = np.arange((-Mx/2.0), (Mx/2.0)+1)
#    Mx = np.arange(-5.0, 6.0) # 5 harmonics
    Mx = np.arange(-1.0, 2.0) # 1 harmonic
    
    # Z axis
    Mz = 0.5*(Nz) + 1.0
    Mz = np.arange(0, Mz)
#    Mz = np.arange(6.0) # 5 harmonics
    Mz = np.arange(2.0) # 1 harmonic



    # X & Z axes
    x = np.linspace(0.0, Lx, Nx)
    z = np.linspace(-Lz/2.0, Lz/2.0, Nz)
    
    
    
    m = N - 2 # number of modes in y axis

    
    generated_ff = np.zeros((Nx, 3*m, Nz), dtype=np.complex128) # multiply y axis by 3 to take [u v w] into account
    sing_vals = np.zeros((len(kx), len(Mx),len(Mz)))
    
    
    for index in range(0, len(kx)):
        # fundamental wavenumbers from wavenumber triplets
        fund_alpha = kx[index]
        fund_beta  = kz[index]
        
        # the stationary wave modes being used to calculate spectral flow field
        # these are also known as the harmonics
        streamwise_stationary_modes = fund_alpha * Mx
        spanwise_stationary_modes = fund_beta * Mz
        
        
        
        if index == 0:
            string_kx = str(fund_alpha)
            string_kz = str(fund_beta)
            string_c = format(c, '.4f')
            string_A = str(amplitudes[index])
        
        
        text01='alpha:'+ str(fund_alpha)+ '  beta:'+ str(fund_beta)+ '    amplitude:'+ str(amplitudes[index])
        print(Fore.RED + text01 + Style.RESET_ALL)
        
        print('kx = mx * alpha        kz = mz * beta')
        
        # loop through the stationary modes
        for ia in range(0, len(streamwise_stationary_modes)):
            for ib in range(0, len(spanwise_stationary_modes)):

                    # now we have the modes.
                    # streamwise_stationary_modes[ia] is the streamwise stationary mode/harmonic
                    # spanwise_stationary_modes[ib] is the spanwise stationary mode/harmonic
                    alpha = streamwise_stationary_modes[ia]
                    beta = spanwise_stationary_modes[ib]
                    
                    # wavenumbers
                    if alpha == 0 or beta == 0:
                        continue
                    
                    
                    # kx = mx * fund_alpha
                    # mx = Mx[ia]
                    # same for Mz
                    text02='(mx)kx: ('+str(Mx[ia])+') '+ str(alpha)+'    (mz)kz: ('+str(Mz[ib])+') '+ str(beta)
                    print(Fore.BLUE + text02 + Style.RESET_ALL)
                    
                    omega = alpha * c
                    
                    C, C_adj, A, w, y_cheb, D1, ResolventA2 = get_state_vectors(N, Re, Lx, Lz, alpha, beta, c)
                    I = np.eye(A.shape[0])
                    L = 1.0j*omega*I + A
                    Linv = inv(L)
                    ResolventA = -1.0* Linv # resolvent
                    H = C*ResolventA*C_adj # transfer function
                    U_spectral, S, V_spectral = svd(H)
                    
                    
                    #===========================================================
                    sing_vals[index, Mx[ia], Mz[ib]] = S[0]
                    # ideally one would have a percentage checker to only 
                    # hold onto the singular values which are within 5/10/15%
                    # of the highest singular value
                    #===========================================================
                    
                    if np.linalg.norm(np.dot( np.dot(U_spectral, np.diag(S)), V_spectral) - H) >= 1e-10:
                        nrm = str(np.linalg.norm(np.dot( np.dot(U_spectral, np.diag(S)), V_spectral) - H))
                        err = 'Something went wrong with the SVD, norm is '+nrm
                        ut.error(err)
                    
                    
                    resolvent_modes = solve(w.T, U_spectral)
                    
                    divergence = 1.0j*alpha*resolvent_modes[0:m, 0] + np.dot(D1, resolvent_modes[m:2*m, 0]) + 1.0j*beta*resolvent_modes[2*m:3*m, 0]
                    div_norm = np.linalg.norm(divergence)
                    
                    if div_norm >= 1e-10:
                        err = 'Something went wrong with the divergence criteria, norm is '+div_norm
                        ut.error(err)
                    
                    
#                    resolvent_modes2 = solve(w.T, U_spectral2)
#                    divergence2 = 1.0j*alpha*resolvent_modes2[0:m, 0] + np.dot(D1, resolvent_modes2[m:2*m, 0]) + 1.0j*beta*resolvent_modes2[2*m:3*m, 0]
#                    div_norm = np.linalg.norm(divergence2)
                    
                    
                    
                    
                    u_tilde = amplitudes[index] * resolvent_modes[:, 0] 
#                    u_tilde = amplitudes[index] * resolvent_modes[:, 0] * S[0]
#                    u_tilde = resolvent_modes[:, 0] * S[0]
#                    u_tilde += np.conjugate(u_tilde)
                    u_tilde = np.asmatrix(u_tilde)
                    
#                    spectral_ff[index, :, ia, ib] = np.squeeze(u_tilde[:,0])
                    
                    
                    physical_ff = np.zeros((len(x), 3*m, len(z)), dtype=np.complex128)
                    physical_ff += ifft(u_tilde[:,0], alpha, beta, c, x, z, t, Lx, Lz, m) # inverse fourier transform
            
                    generated_ff += physical_ff


                    U_u = physical_ff.real[:,   0:m  , :]
                    U_v = physical_ff.real[:,   m:2*m, :]
                    U_w = physical_ff.real[:, 2*m:3*m, :]
                    U = np.zeros((Nd, Nx, Ny, Nz))

        print('')

    
    # Output the flow field as an ASCII file for channelflow to read in.
    
    # Here I need ot construct my geometry file and pack the dictionary such that
    # when it's opened in utils, the writing is as easy as channelflow writing 
    # to file.
    # I store the 4D matrix as follows:
    # u[i, nx, ny, nz]
    # so:
    U = np.zeros((Nd, Nx, Ny, Nz))
    U_u = generated_ff.real[:,   0:m  , :]
    U_v = generated_ff.real[:,   m:2*m, :]
    U_w = generated_ff.real[:, 2*m:3*m, :]
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

    
    
    
    
    L2Norm = np.linalg.norm(U)
    print(np.allclose(L2Norm, np.sqrt(np.sum(np.square(U[:,:,:,:])))))

#    magn = 10.0
#    U *= magn / L2Norm
    

    
    
    
    # Interpolation to go from y_cheb toy_uniform
    Ny = m
    y_uniform = np.linspace(1.0, -1.0, Ny*1.0)
    y_cheb = np.asarray(y_cheb)
    y_cheb = np.squeeze(y_cheb)
    
    
    U_u_uniform = np.zeros((Nx, m, Nz))
    U_v_uniform = np.zeros((Nx, m, Nz))
    U_w_uniform = np.zeros((Nx, m, Nz))
    
    for nx in range(0, Nx):
        for nz in range(0, Nz):
            uprofile = U_u[nx, :, nz] # 1-d vector
            # fill value is the no-slip boundary condition
            fu = interp1d(y_cheb, uprofile, bounds_error=False, fill_value=0.0, kind='cubic') 
            fu = fu(y_uniform)
            U_u_uniform[nx, :, nz] = fu
            
            
            vprofile=U_v[nx, :, nz] # 1-d vector
            fv = interp1d(y_cheb, vprofile, bounds_error=False, fill_value=0.0, kind='cubic') 
            fv = fv(y_uniform)
            U_v_uniform[nx, :, nz] = fv
            
            wprofile=U_w[nx, :, nz] # 1-d vector
            fw = interp1d(y_cheb, wprofile, bounds_error=False, fill_value=0.0, kind='cubic') 
            fw = fw(y_uniform)
            U_w_uniform[nx, :, nz] = fw
            
            

    
#    
#    plt.plot(y_cheb, uprofile, 'r-', y_uniform, fu, 'g--')
#    plt.legend(['data', 'cubic'], loc='best')
#    plt.grid(True)
#    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
        
    gen_ff = {}
    gen_ff['resolvent_flowField'] = U
    
    
#    gen_ff['U'] = U_u_uniform
#    gen_ff['V'] = U_v_uniform
#    gen_ff['W'] = U_w_uniform
    
    gen_ff['U'] = U_u
    gen_ff['V'] = U_v
    gen_ff['W'] = U_w
    
    
    gen_ff['X'] = x
#    gen_ff['Y'] = y_uniform
    gen_ff['Y'] = y_cheb
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
    
    
    
    
    
    
    return gen_ff



def get_state_vectors(N, Re, Lx, Lz, alpha, beta, c):
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
            alpha:  streamwise wavenumber already in 2pialpha/Lx state
            beta:  spanwise wavenumber already in 2pibeta/Lz state
            
     
    OUTPUTS:
             C:  this operator maps the state vector onto the velocity vector
         C_adj:  adjoint of C, maps the forcing vector to the state vector
             A:  state operator
             W:  ClenCurt matrix (Clenshaw-Curtis quadrature)
             y:  grid-points in the y-axis
    """
    
    ## use the fundamental wavenumbers to multiply by the iterative tuples you use
    # i.e. if .geom file says:alpha = 1.14
    # then each alpha is a multiple of alpha
    
    
    
    # Calculate the differentiation matrix, DM, for the resolution in the 
    # y-axis, N. 
    # y_cheb are the interpolated y co-ordinates, i.e. Chebyshev interior points.
    tmp, DM = ps.chebdiff(N, 2)
    
    
    # First derivative matrix
    D1 = DM[0, 1:-1, 1:-1]
    
    # Second derivative matrix
    D2 = DM[1, 1:-1, 1:-1]


    # Fourth derivative matrix and clamped boundary conditions
    y_cheb, D4 = ps.cheb4c(N, False)
    # tmp is the same as y_cheb without endpoints, i.e. no [-1,1]
    
    
    # For the Orr-Sommerfeld equations we need to calculate the derivates
    #            D:  partial_dy
    #          Lap:  D**2.0 - K**2.0         (where K**2.0 = alpha**2.0 + beta**2.0)
    #         dUdy:  first derivative of Uo(y)
    #        dU2dy:  second derivative of Uo(y)
    #            f:  time derivative
    
    
    # Number of modes
    m = N - 2.0 # calculate without the endpoints, otherwise m = N
    I = np.identity(m)
    Z = np.zeros(shape=(m, m))
    K2 = (alpha**2.0) + (beta**2.0)
    Lap = D2 - K2*I #Laplacian
    del_hat_4 = D4 - 2.0*D2*K2 + K2*K2*I
    
    laminarBaseFlow = True
    
    if laminarBaseFlow:
        # Laminar Base flow 
        U = np.identity(m) 
        np.fill_diagonal(U, 1.0 - y_cheb**2.0) # 1 at centreline
        
        dU_dy  = np.identity(m)
        np.fill_diagonal(dU_dy, -2.0*y_cheb)
        
        d2U_dy2 = -2.0

    else:
        # Couette Base flow
        U = np.identity(m)
        np.fill_diagonal(U, y_cheb)
        
        dU_dy  = np.identity(m)
        np.fill_diagonal(dU_dy, 1.0)
        
        d2U_dy2 = 0.0
    
    

    # pg 60 Schmid Henningson eqns 3.29 and 3.30
    SQ_operator = ((Lap/Re) - (1.0j * alpha * U))
    C_operator  = -1.0j*beta*dU_dy
    
    a0=(del_hat_4 / Re)
    a1=( 1.0j * alpha * d2U_dy2 * I)
    a2=(-1.0j * alpha * np.asmatrix(U) * np.asmatrix(Lap))
    
    OS_operator = a0 + a1 + a2
    x0 = solve(Lap, OS_operator)
    
    # Equation 2.7
    # (Moarref, Model-based scaling of the streamwise energy density in 
    # high-Reynolds number turbulent channels, 2013)
    #
    # State operator
    # A = | x0   Z  |
    #     |  C   SQ |
    A = np.vstack((np.hstack((x0,Z)), np.hstack((C_operator, SQ_operator))))

 
    # C maps state vector to the velocity vector
    C = np.vstack((np.hstack(((1.0j/K2) * (alpha*D1), (-1.0j/K2) * (beta*I))), 
                   np.hstack((                  I,                   Z)), 
                   np.hstack(((1.0j/K2) * (beta*D1), ( 1.0j/K2) * (alpha*I)))))
    
    C = np.asmatrix(C)
    
    
    tmp, w = ps.clencurt(N)
    w = np.diag(np.sqrt(w[1:-1]))
    w = np.vstack((np.hstack((w,Z,Z)),
                   np.hstack((Z,w,Z)),
                   np.hstack((Z,Z,w))))
    C = w*C
    
    
    # Adjoint of C maps the forcing vector to the state vector
    C_adj = pinv(C)
    C_adj = C.getH()
    
    
    I = np.eye(A.shape[0])
    omega = alpha * c
    L = 1.0j*omega*I + A
    Linv = inv(L)
    ResolventA = -1.0*inv(L) # resolvent
    
    
    
    # we can alternatively construct the resolvent:
    UminusC2 = np.diag(U) - c
    UminusC = np.identity(m) 
    np.fill_diagonal(UminusC, UminusC2)
    flooby = 1.0j*alpha*UminusC*Lap - 1.0j*alpha*d2U_dy2 - (del_hat_4 / Re)
    
    topleft = solve(Lap, flooby)
    topright = Z
    bottomleft = 1.0j*beta*dU_dy
    bottomright = 1.0j*alpha*UminusC - (Lap / Re)
    
    
    ResolventA2 = np.vstack((np.hstack((topleft,topright)), np.hstack((bottomleft, bottomright))))
    
    
    return C, C_adj, A, w, y_cheb, D1, ResolventA2
        
    
    
def ifft(u_tilde, kx, kz, c, x, z, t, Lx, Lz, m):
    """
    INPUTS:
    u_tilde : velocity as a function of kx y kz omega
         kx : wavenumber in x direction
         kz : wavenumber in z direction
          c : phase speed
          x : streamwise vector
          z : spanwise vector
          t : time instant to probe
         Lx : length of channel
         Lz : width of channel

    OUTPUT:
          u : velocity as a function of x y z t
    """
#    print('FOURIER')
#    for i in range(0, len(fft_signal)):
#        string = str(fft_signal[i].real) + '+' + str(fft_signal[i].imag) + 'j'
#        print(string)
#        tmp = fft_signal[i]
#        
#    print('')
    
    
    u = np.zeros((len(x), 3*m, len(z)), dtype=np.complex128)
    u_tilde = np.asarray(u_tilde)
    u_tilde = u_tilde[:, 0]
    
    
    for iz in range(0, z.shape[0]):
        for ix in range(0, x.shape[0]):
            u[ix, :, iz] = u_tilde * np.exp(1j * ((kx * x[ix]) + (kz * z[iz]) - (kx*c * t)))

    ur = u.real
    ui = u.imag
    
    return u


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
    
    kx *= (2.0 * pi) / Lx
    kz *= (2.0 * pi) / Lz
    
    omega = kx * c
    
    fft_signal = np.zeros((len(x), signal.shape[0], len(z)), dtype=np.complex128)
    
    delta_x = np.abs(x[0] - x[1])
    delta_z = np.abs(z[0] - z[1])
    
#    fft_signal = np.zeros((len(x), len(z)), dtype=np.complex128)# shape could simply be put as signal.shape

#    for iz in range(0, len(signal)):
#        for ix in range(0, x.shape[0]):
#            exp_component = np.exp(1j * ((kx * x[ix]) + (kz * z[iz]) - (omega * t)))
#            a = signal[:, iz]
#            b = a * exp_component[0,0]
#            b = np.squeeze(b)
#            fft_signal[ix, :, iz] = b

    for iz in range(0, len(z)):
        for ix in range(0, len(x)):
            exp_component = np.exp(-1j * ((kx * x[ix]) + (kz * z[iz]) - (omega * t)))
            a = signal[:, :]
            b = a * exp_component * delta_x * delta_z
            b = np.squeeze(b)
            fft_signal[ix, :, iz] = b
    
    
    reciprocal = (1.0 / (Lx * Lz))
    fft_signal *= reciprocal
    
    fft_signal2 = fft_signal[0,:,0] # just for testing.
    
    return fft_signal2




def get_scalars(u_hat, resolvent_modes, w):

    resolvent_modes = np.asmatrix(resolvent_modes)
    
    # Get the complex conjugate of the modes.
    resolvent_modes_star = resolvent_modes.H
    resolvent_modes_star = np.asmatrix(resolvent_modes_star)
    
    # Initialize the scalars vector (shape = Nd*Ny, long vector for u, v, w)
    scalars = np.zeros((u_hat.shape), dtype=np.complex128)
    
    
    # Convert from array to matrix for multiplication later
    w = np.asmatrix(w)    
    u_hat = np.asmatrix(u_hat)
    
    # loop through each resolvent mode
    for i in range(0, resolvent_modes.shape[0]): # from 0 to Nd*Ny - 1
        # take each column of the resolvent modes
        resolvent_modes_col = np.asmatrix(resolvent_modes_star[:,i]) 
        # Dimensions:
        #    1     =  1xN *     NxN     *       Nx1
        tmp = u_hat * w * resolvent_modes_col
        tmp2= u_hat * resolvent_modes_col
        scalars[i] = tmp[0,0]
    
    return scalars
    
#            Mz = beta  * np.arange((np.ceil(-Nz/2.0)+1.0), (np.ceil(Nz/2.0)+1.0), 1.0)
#            Mz = beta  * np.arange((np.ceil(-Nz/2.0)+1.0), (np.ceil(Nz/2.0)+1.0), 1.0)