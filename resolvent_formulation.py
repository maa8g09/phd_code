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
import tests
import utils as ut
import utils_plots as up
import numpy as np


from colorama import Fore, Back, Style
from math import pi
from numpy.linalg import inv
from numpy.linalg import solve
from numpy.linalg import pinv
from numpy.linalg import svd
from scipy.interpolate import interp1d


import matplotlib.pyplot as plt




def resolvent_analysis(geom, Re, kx, kz, c, amplitudes, data, fourdarray):
    """
    The full resolvent formulation is given here, from the generation of the 
    modes to the projection of a channelflow solution onto the modes.
    
    
    INPUTS:
          geom:  grid discretization in X Y Z and # dimensions
            Re:  Reynolds number
            kx:  vector of streamwise wavenumbers
            kz:  vector of spanwise wavenumbers
             c:  phase speed
             A:  amplitude
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
    Mx = geom['Mx']
    Mz = geom['Mz']
    
    generated_ff = np.zeros((len(geom['x']), 3*geom['m'], len(geom['z'])), dtype=np.complex128)     # Complex array
    sing_vals = np.zeros((len(kx), len(Mx),len(Mz)))
    
    
    # Loop through wavenumber triplets
    for index in range(0, len(kx)):
        # Fundamental wavenumbers from wavenumber triplets
        fund_alpha = kx[index]
        fund_beta  = kz[index]
        
        # The stationary wave modes being used to calculate spectral flow field
        # These are also known as the harmonics
        streamwise_stationary_modes = fund_alpha * Mx
        spanwise_stationary_modes   = fund_beta * Mz
        
        if index == 0:
            string_kx = str(fund_alpha)
            string_kz = str(fund_beta)
            string_c = format(c, '.4f')
            string_A = str(amplitudes[index])
        
        
        text01='alpha:'+ str(fund_alpha)+ '  beta:'+ str(fund_beta)+ '    amplitude:'+ str(amplitudes[index])
        print(Fore.RED + text01 + Style.RESET_ALL)
        print('kx = mx * alpha        kz = mz * beta')
        
        # Loop through the stationary modes
        for ia in range(0, len(streamwise_stationary_modes)):
            for ib in range(0, len(spanwise_stationary_modes)):

                    # Wavenumbers
                    alpha = streamwise_stationary_modes[ia]
                    beta = spanwise_stationary_modes[ib]
                    
                    if alpha == 0 or beta == 0:
                        continue
                    
                    # mx = Mx[ia]
                    # same for Mz
                    text02='(mx)kx: ('+str(Mx[ia])+') '+ str(alpha)+'    (mz)kz: ('+str(Mz[ib])+') '+ str(beta)
                    print(Fore.BLUE + text02 + Style.RESET_ALL)
                    
                    omega = alpha * c
                    
                    laminar_mean_flow = True # True: channel, False: couette
                    C, C_adj, A, clencurt_quadrature, y_cheb, D1 = get_state_vectors(geom['Ny'], Re, geom['Lx'], geom['Lz'], alpha, beta, c, laminar_mean_flow)
                    
                    
                    I = np.eye(A.shape[0])
                    L = 1.0j*omega*I + A
                    Linv = inv(L)
                    ResolventA = -1.0* Linv # resolvent of A
                    H = C*ResolventA*C_adj # transfer function
                    
                    vel_modes, singular_values, forcing_modes = svd(H)
                    tests.SVDNorm(vel_modes, singular_values, forcing_modes, H)
                    
                    
                    #===========================================================
                    sing_vals[index, Mx[ia], Mz[ib]] = singular_values[0] # Storing the leading value
                    # Ideally one would have a percentage checker to only 
                    # hold onto the singular values which are within 5/10/15%
                    # of the highest singular value
                    #===========================================================
                    
                    # Non-weighted resolvent modes
                    resolvent_modes = solve(clencurt_quadrature.T, vel_modes)
                    tests.divergence(resolvent_modes, alpha, beta, geom['m'], D1)
                    
                    
                    # u_tilde = chi * Psi
                    u_tilde = amplitudes[index] * resolvent_modes[:, 0] # Rank 1
                    u_tilde = np.asmatrix(u_tilde)
#                    u_tildeH = u_tilde.H
#                    u_tildeH = u_tildeH.T
#                    u_tilde = u_tilde + u_tildeH
                    
                    
                    # Inverse fourier transform
                    physical_ff = np.zeros((len(geom['x']), 3*geom['m'], len(geom['z'])), dtype=np.complex128)
                    physical_ff += ifft(u_tilde[:,0], alpha, beta, c, geom['x'], geom['z'], geom['t'], geom['Lx'], geom['Lz'], geom['m'])
                    generated_ff += physical_ff

        print('\n\n')


    y_uniform=False
    
    outputDic = ut.makeOutputDictionary(generated_ff, geom, y_cheb, y_uniform, string_kx, string_kz, string_c, string_A)

    return outputDic



def resolvent_approximation(data, c, Re, rank):

    string_c = format(c, '.4f')

    gen_ff = {}
    
    
    
    if data['flowField']['is_spectral'] == True:
        u_hat = data['flowField']['spectral']
        u_hat = np.concatenate((u_hat[0,:,:,:],
                                u_hat[1,:,:,:],
                                u_hat[2,:,:,:]), 
                                axis=1)
        
        Mx    = data['geometry']['spectral']['kx']
        Nx    = Mx
        Mz    = data['geometry']['spectral']['kz']
        Nz    = 2*(Mz - 1)
        
        Ny    = data['geometry']['spectral']['Ny']
        N     = Ny + 2
        
        Nd    = data['geometry']['spectral']['Nd']
        
        Lx    = data['geometry']['spectral']['Lx']
        Lz    = data['geometry']['spectral']['Lz']
        
        x     = np.linspace(0.0, Lx, Nx)
        z     = np.linspace(-Lz/2.0, Lz/2.0, Nz)
        
        fund_alpha = data['geometry']['spectral']['alpha']
        fund_beta  = data['geometry']['spectral']['gamma']
        gen_ff['fund_alpha'] = fund_alpha
        gen_ff['fund_beta']  = fund_beta
        
        
        
    elif data['flowField']['is_physical'] == True:
        Nx    = data['geometry']['physical']['Nx']
        Nz    = data['geometry']['physical']['Nz']
        Ny    = data['geometry']['physical']['Ny']
        N     = Ny + 2
        Nd    = data['geometry']['physical']['Nd']
        Lx    = data['geometry']['physical']['Lx']
        Lz    = data['geometry']['physical']['Lz']
        alpha = data['geometry']['physical']['alpha']
        beta  = data['geometry']['physical']['gamma']
        

    
    
    # Stationary nodes along each axis
    # X axis
    # Mx = np.arange(-1.0, 2.0) # 1 harmonic
    Mx = np.arange((-Nx/2.0)+1, (Nx/2.0)+1)
    
    # Z axis
    # Mz = np.arange(2.0) # 1 harmonic
    # Mz = np.arange((-Nz/2.0)+1, (Nz/2.0)+1) # only if youre looking at the whole spectral field.
    Mz = np.arange(0, (Nz/2.0 + 1))
    
    kx = Mx * fund_alpha      # list of wavenumbers to use (modes multiplied by fundamental alpha)
    kz = Mz * fund_beta       # list of wavenumbers to use (modes multiplied by fundamental beta)
    m  = N - 2                # Ny
    
    
    
    u_hat_approx  = np.zeros((len(kx),    3*Ny, len(kz)), dtype=np.complex128)
    sing_vals = np.zeros((      1, len(kx), len(kz)))



    for ikx in range(0, len(kx)):
        for ikz in range(0, len(kz)):
            
            alpha = kx[ikx]
            beta  = kz[ikz]
            
            text02='(mx)kx: ('+str(Mx[ikx])+') '+ str(alpha)+'    (mz)kz: ('+str(Mz[ikz])+') '+ str(beta)
            print(Fore.BLUE + text02 + Style.RESET_ALL)
            
            
            u_hat_approx[ikx, :, ikz] = u_hat[ikx, :, ikz]
            
            if alpha == 0 or beta == 0:
                continue
            
            
            omega=c*alpha
            
            C, C_adj, A, clencurt_quadrature, y_cheb, D1 = get_state_vectors(N, Re, Lx, Lz, alpha, beta, c, False)
            
            # Now make the resolvent.
            I = np.eye(A.shape[0])
            L = 1.0j*omega*I + A
            Linv = inv(L)
            ResolventA = -1.0* Linv                 # resolvent
            H = C*ResolventA*C_adj                  # transfer function
            psi, sigma, phi_star = svd(H)           # SVD
            
            tests.SVDNorm(psi, sigma, phi_star, H)
            
            #===========================================================
            sing_vals[0, Mx[ikx], Mz[ikz]] = sigma[0]
            # ideally one would have a percentage checker to only 
            # hold onto the singular values which are within 5/10/15%
            # of the highest singular value
            #===========================================================
            
            resolvent_modes = solve(clencurt_quadrature.T, psi) # unweighted resolvent modes
            # Divergence test
            tests.divergence(resolvent_modes, alpha, beta, m, D1)
            # Orthogonality test
            tests.orthogonality(psi)
        
            
            
            
            rank = min(rank, 3*m)
            
            # chi  = sigma * eta
            chi = get_scalars(u_hat[ikx, :, ikz], psi, clencurt_quadrature, m, rank)
            chi = np.asarray(chi)
            chi = np.asmatrix(chi)

            
            
            # eta 
            sigma = sigma[:rank]
            sigma = np.asmatrix(np.diag(sigma))
            
            u_tilde_approx = psi[:,:rank] * chi         
            
            
            result = np.asmatrix(u_hat[ikx, :, ikz]).T - u_tilde_approx
            result = np.linalg.norm(result)
            text03='The norm is: '+str(result)
            
            if result <= 1e-10:
                print(Back.GREEN + text03 + Style.RESET_ALL)
            elif result >= 1e-10 and result <= 1e-5:
                print(Back.YELLOW + text03 + Style.RESET_ALL)
            else:
                print(Back.RED + text03 + Style.RESET_ALL)
                
            
            u_hat_approx[ikx, :, ikz] = np.squeeze(np.asarray(u_tilde_approx))

            
            
            
            
    
    # Difference between My Fourier flow field
    # compared to the Gibson flow field.
    diff = np.abs(u_hat - u_hat_approx)
    diff = np.linalg.norm(diff)
    text04='The total flow field norm is: '+str(diff)
    print('\n',Fore.WHITE + Back.MAGENTA + text04 + Style.RESET_ALL,'\n')


    U = np.zeros((Nd, Nx, Ny, Nz))
    generated_ff = np.zeros((Nx, 3*m, Nz))
    
    U_u = generated_ff[:,   0:m  , :]
    U_v = generated_ff[:,   m:2*m, :]
    U_w = generated_ff[:, 2*m:3*m, :]
#    for i in range(0, Nd):
#        for nx in range(0, Nx):
#            for ny in range(0, Ny):
#                for nz in range(0, Nz):
#                    if i == 0: # u direction
#                        U[i, nx, ny, nz] = U_u[nx, ny, nz]
#                    elif i == 1: # v direction
#                        U[i, nx, ny, nz] = U_v[nx, ny, nz]
#                    elif i == 2: # w direction
#                        U[i, nx, ny, nz] = U_w[nx, ny, nz]

    
    
    
    
    L2Norm = np.linalg.norm(U)

    
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
            
            

#    plt.plot(y_cheb, uprofile, 'r-', y_uniform, fu, 'g--')
#    plt.legend(['data', 'cubic'], loc='best')
#    plt.grid(True)
#    plt.show()
    
    
        
    gen_ff['resolvent_flowField'] = U
    
    gen_ff['U'] = U_u_uniform
    gen_ff['V'] = U_v_uniform
    gen_ff['W'] = U_w_uniform
    
#    gen_ff['U'] = U_u
#    gen_ff['V'] = U_v
#    gen_ff['W'] = U_w
    
    
    gen_ff['X'] = x
    gen_ff['Z'] = z
    
    gen_ff['Y'] = y_uniform
#    gen_ff['Y'] = y_cheb

    gen_ff['Nx'] = Nx
    gen_ff['Ny'] = m
    gen_ff['Nz'] = Nz
    gen_ff['Nd'] = Nd

    gen_ff['Lx'] = Lx
    gen_ff['Lz'] = Lz
    
    gen_ff['kx'] = alpha
    gen_ff['kz'] = beta
    gen_ff['c'] = string_c
    gen_ff['A'] = 0
    
    
    U_hat = np.zeros((Nd, Nx, Ny, Nz), dtype=np.complex128)
    U_hat_u = u_hat_approx[:,   0:m  , :]
    U_hat_v = u_hat_approx[:,   m:2*m, :]
    U_hat_w = u_hat_approx[:, 2*m:3*m, :]
    for i in range(0, Nd):
        for mx in range(0, len(Mx)):
            for ny in range(0, Ny):
                for mz in range(0, len(Mz)):
                    if i == 0: # u direction
                        U_hat[i, mx, ny, mz] = U_hat_u[mx, ny, mz]
                    elif i == 1: # v direction
                        U_hat[i, mx, ny, mz] = U_hat_v[mx, ny, mz]
                    elif i == 2: # w direction
                        U_hat[i, mx, ny, mz] = U_hat_w[mx, ny, mz]
    
    
    gen_ff['spectral_ff'] = U_hat
    gen_ff['Mx'] = len(Mx)
    gen_ff['Mz'] = len(Mz)
    gen_ff['Rank'] = rank
    
    return gen_ff




def get_state_vectors(N, Re, Lx, Lz, alpha, beta, c, laminarBaseFlow):
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
    
    
    tmp, clencurt_quadrature = ps.clencurt(N)
    clencurt_quadrature = np.diag(np.sqrt(clencurt_quadrature[1:-1]))
    clencurt_quadrature = np.vstack((np.hstack((clencurt_quadrature,Z,Z)),
                   np.hstack((Z,clencurt_quadrature,Z)),
                   np.hstack((Z,Z,clencurt_quadrature))))
    C = clencurt_quadrature*C
    
    
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
    
    
    return C, C_adj, A, clencurt_quadrature, y_cheb, D1
        
    
    
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









def get_scalars(u_hat, resolvent_modes, clencurt_quadrature, m, rank):

    #========================================================================
    # Projecting with the required amount of column vectors==================
    #========================================================================
    resolvent_modes = np.asmatrix(resolvent_modes)
    psi = resolvent_modes[: , :rank] # column vectors
    
    # Get the complex conjugate of the modes.
    psi_star = psi.H
    
    # Initialize the scalars vector (shape = Nd*Ny, long vector for u, v, w)
    chi = np.zeros((rank, 1), dtype=np.complex128)
    
    # Convert from array to matrix for multiplication later
    clencurt_quadrature = np.asmatrix(clencurt_quadrature)    
    u_hat = np.asmatrix(u_hat).T
    
    chi = psi_star * u_hat
    
    
    #========================================================================
    # Projecting with the full number of column vectors======================
    #========================================================================
#    psi_full = resolvent_modes
#    
#    psi_full_star = psi_full.H
#    
#    chi_full = np.zeros((psi_full.shape[0], 1), dtype=np.complex128)
#    
#    chi_full = psi_full_star * u_hat
    
    
    
    #========================================================================
    # Conclusion: It doesn't make a difference===============================
    #========================================================================
    
    
    return chi
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
def fft_test(Nx, Ny, Nz, U_u_gibson, u_hat_u):
    
        # Lets see if we can get u_hat form u_phy of Gibson
    # Only looking at the streamwise compoenent of the velocity:

    # Going through each y grid point taking each xz-plane:
    test_u_hat = np.zeros((Nx, Ny, Nz), dtype=np.complex128)
    for y in range(0, Ny):
        xzplane = U_u_gibson[:, y, :]
        xzplane_hat = np.fft.fft2(xzplane)
        test_u_hat[:, y, :] = xzplane_hat
        

    test_res      = np.abs(u_hat_u - test_u_hat)
    test_res_norm = np.linalg.norm(test_res)

    
    
    # Alternatively we could Fourier transform it in one go:
    test_u_hat2 = np.zeros((Nx, Ny, Nz), dtype=np.complex128)
    test_u_hat2 = np.fft.fft2(U_u_gibson, axes=(0,2))
    
    test_res2      = np.abs(u_hat_u - test_u_hat2)
    test_res_norm2 = np.linalg.norm(test_res2)
    
    
    return 0
    
    
    
    