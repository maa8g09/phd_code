import numpy as np
import math
import utils
import read_construct as rc
import utils_math as um

from pylab import *
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter
#from matplotlib import animation
from matplotlib import cm as cm
from mpl_toolkits import mplot3d as axes3D

def plot2D_modes(flowField, four_d_array, three_d):

    quiver = True

#    generated_flowField['U'] = U1_u
#    generated_flowField['V'] = U1_v
#    generated_flowField['W'] = U1_w
#    generated_flowField['X'] = x
#    generated_flowField['Y'] = y
#    generated_flowField['Z'] = z
#    generated_flowField['Nx'] = Nx
#    generated_flowField['Ny'] = m
#    generated_flowField['Nz'] = Nz
#    generated_flowField['Nd'] = Nd
#    generated_flowField['Lx'] = Lx
#    generated_flowField['Lz'] = Lz



    # depending on the fourdarray, we plot the following:
    vel = 0
    if four_d_array[0] == 0:
        vel = flowField['U']
        cbar_t = 'u'
    elif four_d_array[0] == 1:
        vel = flowField['V']
#        vel = np.fliplr(vel)
        cbar_t = 'v'
    elif four_d_array[0] == 2:
        vel = flowField['W']
        cbar_t = 'w'
    
    
    
    
    if four_d_array[1] != 'all': #ZY plane
        axis_x = flowField['Z']
        axis_x_t= 'Z'
        axis_y = flowField['Y']
        axis_y_t= 'Y'
        axis_z = flowField['X']
        axis_z_t= 'X'
        
        data_to_plot = vel[four_d_array[1],:,:]

        xcoord=flowField['X'][four_d_array[1]]
        title = 'Generated flow field at x = ' + str(xcoord)
        
        quiv_2D_u = flowField['W'][four_d_array[1], :, :]
        quiv_2D_v = flowField['V'][four_d_array[1], :, :]
        
        
    elif four_d_array[2] != 'all': # XZ plane
        axis_x = flowField['X']
        axis_x_t= 'X'
        axis_y = flowField['Z']
        axis_y_t= 'Z'
        axis_z = flowField['Y']
        axis_z_t= 'Y'
        
        data_to_plot = vel[:,four_d_array[2],:].T
        
        ycoord=flowField['Y'][four_d_array[2]]
        title = 'Generated flow field at y = ' + str(ycoord)
        
        quiv_2D_u = flowField['U'][:, four_d_array[2], :]
        quiv_2D_v = flowField['W'][:, four_d_array[2], :]

        
        
    elif four_d_array[3] != 'all': # XY plane
        axis_x = flowField['X']
        axis_x_t= 'X'
        axis_y = flowField['Y']
        axis_y_t= 'Y'
        axis_z = flowField['Z']
        axis_z_t= 'Z'
        
        data_to_plot = np.fliplr(vel[:,:,four_d_array[3]].T)
        
        
        zcoord=flowField['Z'][four_d_array[3]]
        
        title = 'Generated flow field at z = ' + str(zcoord)
        
        quiv_2D_u = flowField['U'][:, :, four_d_array[3]].T
        quiv_2D_v = flowField['V'][:, :, four_d_array[3]].T

    
    
    
    three_d = False
    
    if three_d:
        # 3D plot, subplot with 2D as well
        
        ################
        # First subplot
        ################
        fig = plt.figure(figsize=plt.figaspect(.5))
        fig.suptitle('2 subplots tutorial')
        # "234" means "2x3 grid, 4th subplot".
        plot1 = fig.add_subplot(1, 2, 1)
        axis_meshg_x, axis_meshg_y = np.meshgrid(axis_x, axis_y)
        origin = 'lower'
        CS = plot1.contourf(axis_meshg_x, 
                            axis_meshg_y, 
                            data_to_plot, 
                            32, 
                            # levels 
                            # [-1, -0.1, 0, 0.1], #alpha=0.5, 
                            cmap=cm.seismic, 
                            origin=origin)
            
        xticks = np.linspace(axis_x[0], axis_x[-1], 6)
        yticks = np.linspace(axis_y[0], axis_y[-1], 8)
        
        plot1.quiver(axis_x, axis_y,quiv_2D_u, quiv_2D_v)
        plot1.set_xlabel(axis_x_t)
        plot1.set_ylabel(axis_y_t)
        plot1.set_xticks(xticks)#, fontsize = 15)
        plot1.set_yticks(yticks)#, fontsize = 15)
        
#        cbar = colorbar(CS)
#        cbar.plot1.set_ylabel(cbar_t)
        
        
        #################
        # Second subplot
        #################
        plot2 = fig.add_subplot(1, 2, 2, projection='3d')
        axis_3D_x, axis_3D_y, axis_3D_z = np.meshgrid(axis_x, axis_y, axis_z)
        plot2.quiver3D(axis_3D_x,
                       axis_3D_y,
                       axis_3D_z,
                       quiv_3D_u,
                       quiv_3D_v,
                       quiv_3D_w,
                       )
        
#        plot2.contour3D(quiv_3D_u,
#                        quiv_3D_v,
#                        quiv_3D_w)
        
#        plot2.plot_surface(quiv_3D_u,
#                            quiv_3D_v,
#                            quiv_3D_w)
        
        plot2.set_xlabel(axis_x_t)
        plot2.set_ylabel(axis_y_t)
        
        xticks = np.linspace(axis_x[0], axis_x[-1], 6)
        yticks = np.linspace(axis_y[0], axis_y[-1], 8)
        
        
        plot1.set_xlabel(axis_x_t)
        plot1.set_ylabel(axis_y_t)
        plot1.set_xticks(xticks)#, fontsize = 15)
        plot1.set_yticks(yticks)#, fontsize = 15)
        
        
    else: 
        xticks = np.linspace(axis_x[0], axis_x[-1], 6)
        yticks = np.linspace(axis_y[0], axis_y[-1], 8)
    
        axis_x, axis_y = np.meshgrid(axis_x, axis_y)
        origin = 'lower'
        plt.figure()
        CS = plt.contourf(axis_x, 
                          axis_y, 
                          data_to_plot, 
                          102, # levels
                          # [-1, -0.1, 0, 0.1], #alpha=0.5,
                          cmap=cm.jet,
                          origin=origin)
        cbar = plt.colorbar(CS)
        cbar.ax.set_ylabel(cbar_t)
    
        if quiver:
            # We will be plotting the quiver plots,
            Q = plt.quiver(axis_x, axis_y,quiv_2D_u, quiv_2D_v), 
#                           angles='xy', scale_units='xy', scale=0.5)
                           
    
        plt.xlabel(axis_x_t)
        plt.ylabel(axis_y_t)
        
        plt.title(title)
        

        
            
        
        plt.xticks(xticks)#, fontsize = 15)
        plt.yticks(yticks)#, fontsize = 15)
    


    plt.show()
    
    
    return
    
    
def plot2D(data_slice, outputDir, fileName, iteration_time):
#            
#    data_slice = {'contourData': contour_data,
#                  'axis_0': a0,
#                  'axis_1': a1,
#                  'axis_0_title': a0title,
#                  'axis_1_title': a1title,
#                  'vel_0': v0,
#                  'vel_1': v1,
#                  'plotTitle': title}
#    
    utils.printSectionHeader()
    utils.printSectionTitle('Plotting a 2D plane')






    x = data_slice['axis_0']
    y = data_slice['axis_1']

    xticks = np.linspace(x[0], x[-1], 6)
    yticks = np.linspace(y[0], y[-1], 8)
    
    
    x, y = np.meshgrid(x, y)
    fig = plt.figure()
#    plt.subplots(1,1)

    vel_min = np.amin(data_slice['contourData'])
    vel_max = np.amax(data_slice['contourData'])

    v = np.linspace(vel_min, vel_max, 10, endpoint=True)
    
    CS = plt.contourf(x, 
                      y, 
                      data_slice['contourData'], 
                      v, 
#                      301, # levels
                      cmap=cm.jet
                      )

    plt.quiver(x,
               y,
               data_slice['vel_0'],
               data_slice['vel_1'],
               color='y'
               )
               
    # Plot text
    plt.xlabel(data_slice['axis_0_title'])
    plt.ylabel(data_slice['axis_1_title'])
    title = ' (fluctuations) at t = ' + iteration_time
    plt.title(data_slice['plotTitle'] + title)
    # the first two value are the x-axis limits 
    # and then the next two values are the y-axis limits
#    plt.axis([0, data_slice['max_0']-1, 0, data_slice['max_1']-1])
    
    
        
    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel(data_slice['cbar_t'])
        
    # Plot grid 
#    plt.grid(True)
    plt.xticks(xticks)#, fontsize = 15)
    plt.yticks(yticks)#, fontsize = 15)
#    plt.show()
    plt.savefig(outputDir + '/' + fileName + '.png')
    plt.close(fig)
    
    return



def plot3D(flow_field, four_d_array, geo, verbosity, quiv, contur):

    utils.printSectionHeader()
    utils.printSectionTitle('Plotting in 3D')
    slice0 = rc.get_data_slice(flow_field, four_d_array , verbosity)
    slice1 = rc.get_data_slice(flow_field, [20, 'all', 'all', 2], verbosity)
    
    
    plane_0 = get_plane_info(slice0, four_d_array)
    
    
    set0, set1, U = convert_List_to_Array(slice0, plane_0, 'U')
    set0, set1, V = convert_List_to_Array(slice0, plane_0, 'V')
    set0, set1, W = convert_List_to_Array(slice0, plane_0, 'W')
    
    set0 = set(plane_0["ax0_coords"])
    set1 = set(plane_0["ax1_coords"])
    set2 = set(plane_0["ax2_coords"])
    print(set1)
    plane_0["ax0_coords"] = list(set0)
    plane_0["ax1_coords"] = list(set1)
    plane_0["ax2_coords"] = list(set2)
    
    X, Y, Z = np.meshgrid(plane_0["ax0_coords"], plane_0["ax1_coords"], plane_0["ax2_coords"])
    

    W = np.zeros(W.shape)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    Q = ax.quiver(X, Y, Z, U, V, W, length=0.1)
    xlabel(plane_0["ax0_label"])
    ylabel(plane_0["ax1_label"])
#    zlabel(plane_0["ax2_label"])
    title(plane_0["fig_title"])
    
    
    return


def plotMatrix(matrix):
    f, (ax1, ax2) = plt.subplots(1, 2)
    ax1.imshow(np.real(matrix), cmap=cm.seismic)
    ax2.imshow(np.imag(matrix), cmap=cm.seismic)
    return
    
def simplePlot(X, Y):
    plt.scatter(X, Y)
    
    
    return
