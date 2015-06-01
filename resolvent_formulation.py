import tools_pseudospectral as ps
import utils_plots as up
import math
import numpy as np

from scipy.linalg import inv
from scipy.linalg import solve
from scipy.linalg import pinv
from scipy.linalg import svd


from pylab import *
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import cm as cm
from mpl_toolkits import mplot3d as axes3D




np.set_printoptions(precision=4)


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
        Lx = 4.0*math.pi #2 periods
        Lz = 2.0*math.pi #1 period
        
        
    else:
        # Or can use channelflow geometry dimensions to generate resolvent modes
        if data['is_physical'] == True:
            x = data['geometry']['physical']['x']
            z = data['geometry']['physical']['z']
            N = data['geometry']['physical']['Ny'] + 2
            Lx = data['geometry']['physical']['Lx']
            Lz = data['geometry']['physical']['Lz']
            kx = np.array(data['geometry']['physical']['Nx'])
            kz = np.array([data['geometry']['physical']['gamma']])
            # Fourier transform the channelflow solution
            u_hat = decompose_flow_field(flowfield, geom)
        elif data['is_spectral'] == True:
            u_hat = data['flowField']['spectral']


    t = 1  # run-time
    
    delta_c = np.abs(c[0] - c[1])
    m = N-2 # number of modes
    
    
    # Initialize an object array to store the generated modes for each
    # wavenumber phase speed combination
    mode = np.zeros((3.0*m, len(kx), len(kz), len(c)), dtype=np.complex128)
    
    # Initialize a 3D matrix for keeping track of the first 
    # singular value per mode
    singular_value = np.zeros((len(kx), len(kz), len(c)))
    
    U = 0
    
    
    # Scalars
    scalars = np.ones((len(x), 3.0*m, len(z)))
    #scalars = np.zeros((len(kx), 3.0*m, len(kz)))
    
    
    for i_kx in range(0, len(kx)):
        for i_kz in range(0, len(kz)):
            
            # account for kx = kz = 0 (turbulent mean)
            if kx[i_kx] == 0 or kz[i_kz] == 0:
                continue
            
            
            C, A, C_adjoint_2, clencurt_quad, y = get_resolvent_operator(N, Re,
                                                                         Lx, Lz,
                                                                         kx[i_kx],
                                                                         kz[i_kz])
                                                                 
            for i_c in range(0, len(c)):
                print 'kx:',kx[i_kx],'    kz:',kz[i_kz],'    c:',c[i_c]
                
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
                    scalars = get_scalars(u_hat, resolvent_modes, clencurt_quad)
                

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
                U += (scalars * sigma[0] * delta_c * phys_flow_field)
                
    generated_flowField = {}
    generated_flowField['resolvent_flowField'] = U
    generated_flowField['U'] = U[0,:,:]
    
    
    
    data['resolvent_flowField'] = U
    
    return U


def get_resolvent_operator(N, Re, Lx, Lz, alpha, beta):
    """
    We are calculating the resolvent operator in this function. The methodology
    followed here is given in the following reference in the "Formulation" 
    section, "Moarreff, 2013, Model-based scaling of ..."
    
    
    INPUTS:
             N:  the number of grid points in the y-axis.
            Re:  Reynolds number
            
     
    OUTPUTS:
             H:  the resolvent operator
             W:  ClenCurt matrix (Clenshaw-Curtis quadrature)
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
    C_row0 = np.hstack((( 1j / K2) * (alpha * D1), (-1j / K2) *  (beta * I)))
    C_row1 = np.hstack((                        I,                       Z))
    C_row2 = np.hstack(( ( 1j / K2) * (beta * D1), ( 1j / K2) * (alpha * I)))
    C = np.vstack((C_row0, C_row1, C_row2))
    C = np.asmatrix(C)
    
    tmp, dy = ps.clencurt(N);
    W = np.diag(np.sqrt(dy[1:-1]))
    W = np.vstack((np.hstack((W,Z,Z)),
                   np.hstack((Z,W,Z)),
                   np.hstack((Z,Z,W))))
    C = W*C;
    
    
    # Adjoint of C maps the forcing vector to the state vector. 
    C_adjoint_2 = pinv(C)
    
    
    return C, A, C_adjoint_2, W, y
    
        
    
    
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
            a = fft_signal[:, :]
            b = a * exp_component
            b = np.squeeze(b)
            signal[i_x, :, i_z] = b
            
    
    return signal

    


def get_scalars(fourier_transformed_flow_field, geom, frdary, U_spectral, sigma, m, W):
   

    U_spectral = np.asmatrix(U_spectral)
    U_spectral_star = U_spectral.getH()
    
    W = np.asmatrix(W)
    
    scalars = np.zeros((geom['Nx'], 3*m, geom['Nz']))
    
    
    
    
    
    return scalars
    
    
    
def decompose_flow_field(flowfield, geom):
    
    #    Think of Fourier Transforms as looking at the 
    #    individual notes that make up a chord.
    #    So if a chord is reala function,
    #    by fourier transforming a chord you are left with the 
    #    notes that made it, not just the final result (the chord)
    fourier_transformed_flow_field = np.zeros((3, geom['Nx'], geom['Ny'], geom['Nz']), dtype=np.complex128)

    
    
    
    # Fourier transform in x and z
    # So I take a 2D fourier transform of each XZ plane.
    for ny in range(0, geom['Ny']):
        fourier_transformed_flow_field[0,:,ny,:] = np.fft.fft2(flowfield['flow_field'][0,:,ny,:])
        fourier_transformed_flow_field[1,:,ny,:] = np.fft.fft2(flowfield['flow_field'][1,:,ny,:])
        fourier_transformed_flow_field[2,:,ny,:] = np.fft.fft2(flowfield['flow_field'][2,:,ny,:])

    
    return fourier_transformed_flow_field
    
