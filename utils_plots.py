import numpy as np
import utils
import read_construct as rc
import utils_math as um

from pylab import *
from matplotlib import pyplot as plt
from matplotlib import animation
from matplotlib import cm as cm
from mpl_toolkits import mplot3d as axes3D

def plot2D_modes(res_field, fourdarray):


    Z = len(res_field['Z'])
    Z = np.arange(Z)


    Y = len(res_field['Y'])
    Y = np.arange(Y)

    
    Z, Y = np.meshgrid(Z, Y)
    plane_image = plt.subplots(1,1)

#    plane_image = plt.contourf(
#                               Z, 
#                               Y,
#                               four_d['U']
#                               )
                               
    plane_image = plt.quiver(Z,
                             Y,
                             res_field['W'][0,:,:],
                             res_field['V'][0,:,:])
                             
    plane_image = plt.imshow(res_field['U'][0,:,:])
    plane_image.set_cmap('coolwarm')
    plt.colorbar()
    plt.xlabel('Z')
    plt.ylabel('Y')
    plt.title('Resolvent Modes')
    
    return
    
    
def plot2D(data_slice):
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

    x = np.arange(len(data_slice['axis_0']))
    y = np.arange(len(data_slice['axis_1']))

    
    x, y = np.meshgrid(x, y)

    plt.subplots(1,1)
    plt.quiver(x,
               y,
               data_slice['vel_0'],
               data_slice['vel_1']
               )
#    plt.xlabel(data_slice['axis_0_title'])
#    plt.ylabel(data_slice['axis_1_title'])
#    plt.title(data_slice['plotTitle'])
#               colour,
#               clim=[min(colour), max(colour)],
#               pivot='mid')

#    v_cont_normalized = um.normalize(v_cont)
    
    # find the maximum value in the array and normalize everything according to it
#    max_in_each_row_list = data_slice['vel_data'].max(axis=1)
#    max_vel = max_in_each_row_list.max()
#    min_in_each_row_list = data_slice['vel_data'].min(axis=1)
#    min_vel = min_in_each_row_list.min()
#    
#    delta = max_vel - min_vel
#    # now we do teh normalizing...
#    for i in range(0, data_slice['vel_data'].shape[0]):
#        for j in range(0, data_slice['vel_data'].shape[1]):
#            # range [-1, 1]
##            data_slice['vel_data'][i,j] = ((data_slice['vel_data'][i,j] - min_vel) / (delta) - 0.5) * 2.0
#            # range [0, 1]
#            data_slice['vel_data'][i,j] = (data_slice['vel_data'][i,j] - min_vel) / (delta)
#            

    plane_image = plt.imshow(data_slice['contourData'])
    plane_image.set_cmap('coolwarm')
    plt.colorbar()
    
    # Plot text
    plt.xlabel(data_slice['axis_0_title'])
    plt.ylabel(data_slice['axis_1_title'])
    plt.title(data_slice['plotTitle'])
    # the first two value are the x-axis limits 
    # and then the next two values are the y-axis limits
#    plt.axis([0, data_slice['max_0']-1, 0, data_slice['max_1']-1])
    
    # Plot grid 
    plt.grid(True)
    
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
    print set1
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
    ax1.imshow(np.real(matrix), cmap=cm.coolwarm)
    ax2.imshow(np.imag(matrix), cmap=cm.coolwarm)
    return
    
def simplePlot(X, Y):
    plt.scatter(X, Y)
    
    
    return
