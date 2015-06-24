import utils as ut
import numpy as np
import math


def diff_check(directory, ff_1, ff_2):
    
    
    u1 = np.zeros((3, 32, 35, 32))
    u2 = u1
    
    
    # ORIGINAL
    f = ut.openFile(directory + '/' + ff_1)
    for i, line in enumerate(f):
            values = line.split()
            nx = int(values[0])
            ny = int(values[1])
            nz = int(values[2])
            nd = int(values[3])
            vel = float(values[4])
            
            u1[nd, nx, ny, nz] = vel
    f.close()


    # AFTER SOLVER
    f = ut.openFile(directory + '/' + ff_2)
    for i, line in enumerate(f):
            values = line.split()
            nx = int(values[0])
            ny = int(values[1])
            nz = int(values[2])
            nd = int(values[3])
            vel = float(values[4])
            
            u2[nd, nx, ny, nz] = vel
    f.close()
    
    threshold = 1e-8
    
    Nd = u1.shape[0]
    Nx = u1.shape[1]
    Ny = u1.shape[2]
    Nz = u1.shape[3]
    
    differences = {}
    
    for nd in range(0, Nd):
        for nx in range(0, Nx):
            for ny in range(0, Ny):
                for nz in range(0, Nz):
                    diff = np.abs( u1[nd, nx, ny, nz] - u2[nd, nx, ny, nz] )
                    
                    if diff >= threshold:
                        label=str(nd)+'.'+str(ny)+'.'+str(ny)+'.'+str(nz)
                        differences[label] = diff
                        
    if not differences:
        print 'Done checking all differences, there were...'
        print 'No Differences :)'
        
    else:
        print 'Done checking all differences, there were...'
        print 'Some differences...'
        print 'Check the dictionary in debug mode'
        
    
    return 0
          

directory='/home/arslan/Documents/work/channelflow-related/database_solns/HKW/equilibria/eq4/nonlinear_solver/best_soln_ascii/differences'
ff_1='eq4.asc'
ff_2='ubest.asc'
diff_check(directory, ff_1, ff_2)