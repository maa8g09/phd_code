import tools_pseudospectral as ps
import utils_plots as up
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
from matplotlib import cm as cm
from mpl_toolkits import mplot3d as axes3D




np.set_printoptions(precision=4)


def main_resolvent_analysis(flowfield, geom, frdary, N, Re, kx, kz, c, independent):
    """
    The full resolvent formulation is given here, from the generation of the 
    modes to the projection of a channelflow solution onto the modes.
    
    INPUTS:
       flowfield:  downloaded solution from channelflow.org
          geom:  geometry variables of the channel flow solution read in
        frdary:  a 4D array used to plot slices from channel flow solutions
             N:  resolution in y axis
            Re:  Reynolds number
            kx:  vector of streamwise wavenumbers
            kz:  vector of spanwise wavenumbers
             c:  phase speed
    
    OUTPUTS:
         fourd:  projected solution / generated modes.
    """
    independent = False
    
    
    
    
    
    if independent:
        # Can use independent geometrical values
        x = np.arange(0.0, 2.0*math.pi, 0.1)
        z = np.arange(-2.0, 2.1, 0.1)
        Lx = 4.0*math.pi
        Lz = 2.0*math.pi
        
        
    else:
        # Or can use channelflow geometry dimensions to generate resolvent modes
        x = flowfield['X']
        z = flowfield['Z']
        N = geom['Ny'] + 2
        Lx = geom['Lx']
        Lz = geom['Lz']
        kx = np.array([geom['alpha']])
        kz = np.array([geom['gamma']])



    flowfield = np.random(random((geom['Nx'], geom['Ny'], geom['Nz'])))
    

    t = 1  # run-time
    
    
    
    delta_c = np.abs(c[0] - c[1])
    m = N-2 # m = number of modes generated.
    
    
    # Initialize an object array to store the generated modes for each
    # wavenumber phase speed combination
    mode = np.zeros((3.0*m, len(kx), len(kz), len(c)), dtype=np.complex128)
    
    # Initialize a 3D matrix for keeping track of the first 
    # singular value per mode
    singular_value = np.zeros((len(kx), len(kz), len(c)))
    
    
    
    U = 0
    
    # Fourier transform the channelflow solution
    fourier_transformed_flow_field = decompose_flow_field(flowfield, geom)
    scalars = np.zeros((len(kx), 3.0*m, len(kz)))
    
    
    
    for i_kx in range(0, len(kx)):
        for i_kz in range(0, len(kz)):
            for i_c in range(0, len(c)):
                omega = kx[i_kx] * c[i_c]

                # account for kx = kz = 0 (turbulent mean)
                if kx[i_kx] == 0 or kz[i_kz] == 0:
                    continue
                
                resolvent, clencurt_quad, y = get_resolvent_operator(N, 
                                                                     Re, 
                                                                     geom, 
                                                                     kx[i_kx], 
                                                                     kz[i_kz], 
                                                                     omega)
                
                # Perform SVD on the resolvent to get the forcing, singluar and
                # response modes. This also gives the flow field in spectral 
                # space.
                U_spectral, sigma, V_spectral = svd(resolvent)
                

                if independent:
                    scalars = np.ones((len(x), 3.0*m, len(z)))
                else:
                    scalars = get_scalars(fourier_transformed_flow_field, geom, frdary, U_spectral, clencurt_quad)
                    
                    
                # Solve linear equation to get the physical velocities in
                # streamwise and wall-normal directions. 
                u_physical = solve(clencurt_quad.T, U_spectral)
                
                # Store the generated modes and the first singular value 
                # (to be used below)
                mode[:, i_kx, i_kz, i_c] = np.asarray(np.squeeze(u_physical[:,0]))
                singular_value[i_kx, i_kz, i_c] = sigma[0]
                
                # Conver the flow field to physical space.
                phys_flow_field = fourier_to_physical(u_physical[:, 0],
                                                     kx[i_kx],
                                                     kz[i_kz],
                                                     c[i_c],
                                                     x, z, t,
                                                     Lx,
                                                     Lz)
                
                # generated flow field is the velocity vector U 
                # which contains (u,v,w)
                U = U + (scalars * sigma[0] * delta_c * phys_flow_field) # multiplied by a weighting coefficient, these are the scalras I get from the projcetion shit I gtta do.
                print 'kx:',kx[i_kx],'    kz:',kz[i_kz],'    c:',c[i_c]
                














#    
#    fig = plt.figure()
#    
#    ims = []
#
#    Ut = np.zeros((3.0*m, len(z)), dtype=np.complex128)
#    
#    for time in range(0, t):
#        time = time# * 0.1
#        Ut = Ut * 0.0
#        print '\n\n\n\n\n\n\n====================================Time is:',time
#        for i in range(0, len(kx)):
#            for j in range(0, len(kz)):
#                for q in range(0, len(c)):
#                    
#                    if kx[i] == 0 or kz[j] == 0:
#                        continue
#
#                    transformed_ff_t = fourier_to_physical(mode[:, i, j, q], kx[i], kz[j], c[q], [frdary[1]], z, time, Lx, Lz)
#                    a = singular_value[i,j,q]
#                    b = transformed_ff_t[0,:,:]
#                    e = scalars[frdary[1],:,:]
#                    Ut = Ut + a * delta_c * e * b
#                    print 'kx:',kx[i_kx],'    kz:',kz[i_kz],'    c:',c[i_c]
#    
#    Ut = Ut.real
#    u = Ut[0:m,:]
#    v = Ut[m:2*m,:]
#    w = Ut[2*m:3*m,:]
#    
#    fourd = {'X': x, 'Y': y, 'Z': z, 'U': u, 'V': v, 'W': w}
#        
    fourd = 0
    
    return fourd


def get_resolvent_operator(N, Re, geo_var, alpha, beta, omega):
    """
    INPUTS:
             N:  the number of grid points in the y-axis.
            Re:  Reynolds number
     
    OUTPUTS:
             H:  the resolvent operator
             W:  ClenCurt matrix (Clenshaw-Curtis quadrature)
    """
    # 
    alpha = alpha * (2.0 * math.pi) / geo_var['Lx']
    beta = beta * (2.0 * math.pi) / geo_var['Lz']
    
    
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
    
    # The resolvent of A
    RA = inv(-1j*omega*np.eye(m*2) - A)

    # Transfer Function
    H = C * RA * C_adjoint_2
    
    return H, W, y
    
        
    
    
def fourier_to_physical(u_hat, kx, kz, c, x, z, t, Lx, Lz):
    """
    INPUTS:
      u_hat  : Fourier mode vector (u, v, w), fourier coeffs
         kx  : wavenumber in x direction
         kz  : wavenumber in z direction
      omega  : temporal frequency
          x  : streamwise vector
          z  : spanwise vector
          t  : time instant to probe
         Lx  : length of channel
         Lz  : width of channel
    
    OUTPUT:
          U  : physical flowfield
          
    """
    
    omega = kx * c
    
    # initialize flow field with zeros
    U = np.zeros((len(x), u_hat.shape[0], len(z)))
    
    for k in range(0, len(z)):
        for i in range(0, len(x)):
            a = (kx * 2.0*math.pi/Lx * x[i])
            b = (kz * 2.0*math.pi/Lz * z[k])
            c = (omega*t)
            expo = np.exp(1j*(a + b - c))
            U[i, :, k] = np.asarray(np.squeeze(u_hat)) * expo
    
    return U
    


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
    
