import math
import numpy as np
from scipy.linalg import toeplitz
from scipy.fftpack import ifft


def chebdiff(N, M, debug):
    #  The function computes the differentiation
    #  matrices D1, D2, ..., DM on Chebyshev nodes.
    #
    #
    #  Input:
    #   -    N:  Size of differentation matrix
    #   -    M:  Number of derivatives required (integer).
    #   - Note:  0 < M <= N-1
    #
    #
    #  Output:
    #   -   DM:  DM(1:N, 1:N, e11) contains ell-th derivative
    #            matrix, e11 = 1..M
    #
    #
    #  The code implements two strategies for enhanced accuracy suggested by
    #  W. Don and S. Solomonoff in SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268
    #  (1994).
    #
    #  The two strategies are:
    #    (a) the use of trigonometric identities to avoid the computation of
    #        differences x(k)-x(j)
    #
    #    (b) the use of the "flipping trick" which is necessary since
    #        sin t can be computed to high relative precision when t is small
    #        whereas sin (pi-t) cannot.
    #
    #  Note added May 2003: It may, in fact, be slightly better not to
    #  implement the strategies (a) and (b).
    #
    #  Please consult the following paper for details:
    #    "Spectral Differencing with a Twist", by R. Baltensperger and
    #    M.R. Trummer, to appear in SIAM J. Sci. Comp.
    #
    #  J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified by
    #  JACW, May 2003.


    # Identity matrix with logical matrix
    I = np.identity(N)
    #I_logical = np.matrix(np.eye(N, dtype=bool))
    if debug:
        print '\n\nI:\n', I
        #print '\n\nI_logical:\n', I_logical


    # n1 n2 are indices used for the flipping trick
    n1 = math.floor(N / 2.0)
    n2 = math.ceil(N / 2.0)
    if debug:
        print '\n\nn1:\n', n1
        print '\n\nn2:\n', n2


    # Get the theta vector
    k = np.arange(N)[:, np.newaxis]
    theta = k * math.pi / (N - 1.0)
    if debug:
        print '\n\nk:\n', k
        print '\n\ntheta:\n', theta


    # Chebyshev points
    x = np.sin(math.pi * np.arange(N - 1, -1.0 * N, -2)[:, np.newaxis] / (2.0 * (N - 1)))
    if debug:
        print '\n\nx:\n', x

    T = np.tile(theta / 2.0, (1.0, N))
    if debug:
        print '\n\nT:\n', T


    # tplust = T.T + T
    # tminust= T.T - T
    # print '\n\nT.T + T:\n', tplust
    # print '\n\nT.T - T:\n', tminust


    # Trigonometric identity
    DX = 2.0 * np.sin(T.T + T) * np.sin(T.T - T)
    if debug:
        print '\n\nDX:\n', DX

    # Flipping trick
    DX = np.concatenate((DX[:n1], -1.0 * np.flipud(np.fliplr(DX[:n2]))), axis = 0)
    if debug:
        print '\n\nDX:\n', DX


    # Put 1s in the main diagonal of DX
    np.fill_diagonal(DX, 1.0)
    if debug:
        print '\n\nDX (after altering diagonals):\n', DX


    # C is the matrix with entries c(k)/c(j)
    C = toeplitz((-1.0)**k)
    C[0:1] = C[0:1] * 2.0       # multiply the first row by 2
    C[N-1:N] = C[N-1:N] * 2.0   # multiply the  last row by 2
    C[:,0] = C[:,0] / 2.0       #   divide the first col by 2
    C[:,N-1] = C[:,N-1] / 2.0   #   divide the  last col by 2
    if debug:
        print '\n\nC:\n', C

    # Z contains entries 1/(x(k)-x(j))
    Z = 1.0 / DX
    if debug:
        print '\n\nZ:\n', Z
    # with zeros on the diagonal.
    np.fill_diagonal(Z, 0.0)
    if debug:
        print '\n\nZ (after zeroing):\n', Z

    D = np.matrix(np.identity(N), copy=False)
    if debug:
        print '\n\nD:\n', D


    DM = np.empty((M, N, N))
    if debug:
        print '********************************************************************'
        print '\n\nEntering the for-loop:'
    for ell in range(0,M):
        if debug:
            print '\n__________________________________________________________________\nell:\n', ell
        A = np.tile(np.diag(D), N).reshape(N,N).T
        if debug:
            print '\n\nA:\n', A
        B = np.multiply(C, A) - D
        if debug:
            print '\n\nB:\n', B
        G = np.multiply(Z, B)
        if debug:
            print '\n\nG:\n', G

        D = (ell + 1) * G
        #D1 = (ell + 1) * np.multiply(Z, np.multiply(C, np.tile(np.diag(D), N).reshape(N,N).T) - D)
        if debug:
            print '\n\n D:\n', D
            #print '\n\n D1:\n', D1

        #sum up the rows of D
        row_sums = -1.0 * D.sum(axis=1)
        # then assign the diagonals the value of the sum of the row.
        np.fill_diagonal(D, row_sums)
        if debug:
            print '\n\nD (after diagonal stuff):\n', D
        DM[ell,:,:] = D # DM[ell::] = D
        if debug:
            print '\n\nDM:\n', DM
    
    
    return x, DM






def cheb4c(N, debug):
    #  The function [x, D4] =  cheb4c(N) computes the fourth 
    #  derivative matrix on Chebyshev interior points, incorporating 
    #  the clamped boundary conditions u(1)=u'(1)=u(-1)=u'(-1)=0.
    #
    #  Input:
    #   -  N:     N-2 = Order of differentiation matrix.  
    #                   (The interpolant has degree N+1.)
    #
    #  Output:
    #   -  x:     Interior Chebyshev points (vector of length N-2)
    #   - D4:     Fourth derivative matrix  (size (N-2)x(N-2))
    #
    #
    #  The code implements two strategies for enhanced accuracy suggested by 
    #  W. Don and S. Solomonoff in SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 
    #  (1994).
    #  
    #  The two strategies are: 
    #    (a) the use of trigonometric 
    #        identities to avoid the computation of differences 
    #        x(k)-x(j) and 
    #
    #    (b) the use of the "flipping trick"
    #        which is necessary since sin t can be computed to high
    #        relative precision when t is small whereas sin (pi-t) cannot.
    #  
    #  J.A.C. Weideman, S.C. Reddy 1998.
        
        
    # Identity matrix with logical matrix
    I = np.matrix(np.identity(N - 2.0), copy=False)
    I_logical = np.matrix(np.eye(N - 2.0, dtype=bool))
    if debug:
        print '\n\nI:\n', I
        print '\n\nI_logical:\n', I_logical

    
    # n1 n2 are indices used for the flipping trick
    n1 = math.floor((N / 2.0) - 1)
    n2 = math.ceil((N / 2.0) - 1)
    if debug:
        print '\n\nn1:\n', n1
        print '\n\nn2:\n', n2
    
    
    # Get the theta vector
    k = np.arange(1, N - 1.0)[:, np.newaxis]
    theta = k * math.pi / (N - 1.0)
    if debug:
        print '\n\nk:\n', k
        print '\n\ntheta:\n', theta
    
    
    # Chebyshev points
    x = np.sin(math.pi * np.arange(N - 3.0, 1.0 - N, -2.0)[:, np.newaxis] / (2.0 * (N - 1.0)))
    if debug:
        print '\n\nx:\n', x
    
    
    # s = sin(theta)
    s = np.concatenate((np.sin(theta[:n1]), np.flipud(np.sin(theta[:n2]))), axis = 0)
    if debug:
        print '\n\ns:\n', s
    
    
    # Compute weight function and its derivatives
    alpha = s**4.0
    beta1 = -4.0 * (s**2.0) * (x / alpha)
    beta2 = 4.0 * (3.0 * x**2.0 - 1.0) / alpha
    beta3 = 24.0 * x / alpha
    beta4 = 24.0 / alpha
    B = np.array([beta1.T, beta2.T, beta3.T, beta4.T])
    if debug:
        print '\n\nalpha:\n', alpha
        print '\n\nbeta1:\n', beta1
        print '\n\nbeta2:\n', beta2
        print '\n\nbeta3:\n', beta3
        print '\n\nbeta4:\n', beta4
        print '\n\nB:\n', B

    T = np.tile(theta / 2.0, (1.0, N - 2.0))
    if debug:
        print '\n\nT:\n', T
    

    # Trigonometric identity
    DX = 2.0 * np.sin(T.T + T) * np.sin(T.T - T)
    if debug:
        print '\n\nDX:\n', DX
    
    # Flipping trick
    DX = np.concatenate((DX[:n1], -1.0 * np.flipud(np.fliplr(DX[:n2]))), axis = 0)
    if debug:
        print '\n\nDX:\n', DX

    
    # Put 1s in the main diagonal of DX
    np.fill_diagonal(DX, 1.0)
    if debug:
        print '\n\nDX (after altering diagonals):\n', DX
    
    
    # Compute the matrix with entries c(k)/c(j)
    ss = (s**2.0) * ((-1.0)**k)
    if debug:
        print '\n\nss:\n', ss
    S = np.tile(ss, N - 2.0)
    if debug:
        print '\n\nS:\n', S
    C = S / S.T
    if debug:
        print '\n\nC:\n', C
    
    
    # Z contains entries 1/(x(k)-x(j))
    Z = 1.0 / DX
    if debug:
        print '\n\nZ:\n', Z
    # with zeros on the diagonal.
    np.fill_diagonal(Z, 0.0)
    if debug:
        print '\n\nZ:\n', Z
    
    
    # X is same as Z', but with diagonal entries removed
    X = Z.T
    if debug:
        print '\n\nX:\n', X

    
    # INSTEAD OF THIS, think of the arrays as vertical lists!!!, 
    # then when you get rid of the zeroes, 
    # you can concatenate the lists together side by side.
    Y = np.zeros((N-3, N-2.0))
    for i in range(0, int(N)-2):
        Y[:,i] = filter(lambda a: a != 0.0, X[:,i])
    
    X = Y
    if debug:
        print '\n\nX:\n', X


    # Initialize Y and D vectors.
    # Y contains matrix of cumulative sums,
    # D scaled differentiation matrices.
    Y = np.ones((N-3, N-2.0))
    if debug:
        print '\n\nY:\n', Y
    
    
    D = I
    A = np.ones((N - 2.0, N - 2.0))
    DM = np.empty((4, N - 2.0, N - 2.0))
    
    if debug:
        print '*************************************************************************\n\nEntering the for-loop:'
    for ell in range(0, 4):
        if debug:
            print '____________________________________________________________________________'
            print 'ell:\n', ell+1
            
        test_01 = B[ell, :]
        test_02 = (ell + 1) * Y[0:N-3, :] * X
        test_03 = np.concatenate((B[ell, :], (ell + 1) * Y[0:N-3, :] * X), axis=0)
        if debug:
            print '\ntest_03:\n', test_03
        
        # Recursion for diagonals
        Y = np.cumsum(np.concatenate((B[ell, :], (ell + 1) * Y[0:N-3, :] * X), axis=0), axis=0)
#        Y = A
#        print '\n\nA:\n', A
        if debug:
            print '\n\nY:\n', Y
        
        # Off-diagonal
        a = np.tile(np.diag(D), N - 2.0).reshape(N - 2.0, N - 2.0).T
        b = np.multiply(C, a) - D
        g = np.multiply(Z, b)
        D = (ell + 1) * g
        
        if debug:
            print '\n\nD:\n', D
        # Correct the diagonal, put the last row of Y into the empty diagonal elements in D
#        print '\n\nA last row:\n', A[N - 3, :]
        np.fill_diagonal(D, Y[-1, :])
        if debug:
            print '\n\nD (after filling diagonals:\n', D
        DM[ell,:,:] = D # DM[ell::] = D

        
    if debug:
        print '\n\nDM:\n', DM


    D4 = DM[3,:,:]
    if debug:
        print '\n____\nD4:\n', D4
      
    return x, D4
    
    
    


'''
CLENCURT

    nodes x (Chebyshev points) and weights w for Clenshaw-Curtis
    quadrature from Spectral Methods in MATLAB
'''
def clencurt(n1):
    # Trefethen's N is different from W&R's
    if n1 is 1:
        x = 0
        w = 2
    else:
        n = n1 - 1
        C = np.zeros((n1, 2))
        k = 2 * (1 + np.arange(np.floor(n / 2)))
        C[::2, 0] = 2 / np.hstack((1, 1 - k * k))
        C[1, 1] = -n
        V = np.vstack((C, np.flipud(C[1:n, :])))
        F = np.real(ifft(V, n=None, axis=0))
        x = -1.0 * F[0:n1, 1]  # I added the -1 at the beginning of this operation
        w = np.hstack((F[0, 0], 2 * F[1:n, 0], F[n, 0]))

    return x, w


def chebint(fk, x, debug):
    #  The function p = chebint(fk, x) computes the polynomial interpolant
    #  of the data (xk, fk), where xk are the Chebyshev nodes.  
    #  Two or more data points are assumed.
    #
    #  Input:
    #   - fk:  Vector of y-coordinates of data, at Chebyshev points 
    #          x(k) = cos((k-1)*pi/(N-1)), k = 1...N.
    #   - x:   Vector of x-values where polynomial interpolant is to be evaluated.
    #
    #  Output:
    #   - p:    Vector of interpolated values.
    #
    #
    #  The code implements the barycentric formula; see page 252 in
    #  P. Henrici, Essentials of Numerical Analysis, Wiley, 1982.
    #  (Note that if some fk > 1/eps, with eps the machine epsilon,
    #  the value of eps in the code may have to be reduced.)
    #
    #  J.A.C. Weideman, S.C. Reddy 1998

    print '\n\n  --== Start chebint ==--'
     
    
    # Write some code to make sure that thte data passed in is column vecotors
    fk = np.expand_dims(fk, axis=1)
    x = np.expand_dims(x, axis=1)
    if debug:
        print '\n\nfk:\n', fk
        print '\n\nx:\n', x
        

    N = len(fk)
    M = len(x)
    if debug:
        print '\n\nN:\n', N
        print '\n\nM:\n', M
        

    # Compute the chebyshev points
    xk = np.sin(math.pi * np.arange(N - 1, -1.0 * N, -2.0)[:, np.newaxis] / (2.0 * (N - 1)))
    if debug:
        print '\n\nxk:\n', xk

    
    powers = np.arange(N)[:, np.newaxis]
    w = np.ones((N,1)) * ((-1.0)**(powers))
    if debug:
        print '\n\nw:\n', w
    
    w[0] = w[0] / 2.0
    w[-1] = w[-1] / 2.0
    if debug:
        print '\n\nw:\n', w.shape
    
    
    # Compute quantities x-x(k)
#    var1 = np.tile(x, N)
#    var2 = np.tile(xk, M).T
#    if debug:
#        print'\n\nvar1:\n', var1, '\n\n', var1.shape
#        print'\n\nvar2:\n', var2, '\n\n', var2.shape
    D = np.tile(x, N) - np.tile(xk, M).T
    if debug:
        print '\n\nD:\n', D, '\n\n', D.shape


    # Compute the reciprocals
    D = 1./(D + np.spacing(1)*(D == 0.0))
    if debug:
        print '\n\nD:\n', D, '\n\n', D.shape


    # Evaluate interpolant as matrix-vector products
#    var3 = w * fk
#    if debug:
#        print'\n\nvar3:\n', var3, '\n\n', var3.shape  
#    var4 = np.dot(D,w)
#    if debug:
#        print'\n\nvar4:\n', var4, '\n\n', var4.shape
    p = np.dot(D, w * fk) / (np.dot(D,w))
    if debug:
        print '\n\np:\n', p, '\n\n', p.shape


    print '\n\n  --== End ==--\n\n'
    
    return p
