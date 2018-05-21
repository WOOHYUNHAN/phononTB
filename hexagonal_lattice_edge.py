import os.path
import time
from math import *
from math import sqrt
from numpy.linalg import *
import numpy as np
import matplotlib.pyplot as plt
import cmath
from matplotlib.collections import LineCollection
import matplotlib.cm as cm
from draw_phononTB_HAN import *
#from numeric import *




lat_hex = [[1.2333638807431999,   -2.1362489056675500], [1.2333638807431999,   2.1362489056675500]]
orb_hex = [[1.0/3, 2.0/3], [2.0/3, 1.0/3]]
#lat_hex = [[1.0,0.0],[0.5,np.sqrt(3.0)/2.0]]
#orb_hex = [[1./3.,1./3.],[2./3.,2./3.]]
mass_hex = [12.0, 12.0]
fc_nn  =  [[-1.0,0.0],[0.0,-1.0]]
fc_nn1 = [[-22.19684176,   -7.699687431],[-7.699687431,   -13.306]]
fc_nn2 = [[-22.19684176 ,    7.699687431],[ 7.699687431,   -13.306]]
fc_nn3 = [[-8.860588573805234,    -0.000000000000000],[-0.000000000000000 ,  -26.642257254446765]]
fc_nnn = [[0.1,0.0],[0.0,0.3]]
V_info_hex = [-1.0, -0.1]

###############################################################################
FC_hex_edge = ForceConstant(2, 2)

FC_hex_edge.set_geometry(lat_hex, orb_hex, mass_hex)
FC_hex_edge.set_hopping(0,1,[0.0,0.0],V_info_hex)
FC_hex_edge.set_hopping(0,1,[-1.0,0.0],V_info_hex)
FC_hex_edge.set_hopping(0,1,[0.0,1.0],V_info_hex)
#FC_hex_edge.set_hopping(0,1,[0.0,0.0],V_info_hex)
#FC_hex_edge.set_hopping(1,0,[1.0,0.0],V_info_hex)
#FC_hex_edge.set_hopping(1,0,[0.0,1.0],V_info_hex)



FC_hex_edge.make_edge(31,0)
FC_hex_edge.set_acoustic_sum_rule()


#FC_hex_edge.set_fc_direct(0, 0, [0, 0], np.array([[ 0.775     ,  -0.38971143],   [ -0.38971143,  0.325    ]]))
#FC_hex_edge.set_fc_direct(61, 61, [0, 0], np.array([[ 0.775     ,  -0.38971143],   [ -0.38971143,  0.325    ]]))

FC_hex_edge.print_info()

###############################################################################
alpha0 = [[0.0,0.0, 0.00]]
alpha2 = [[0.0,0.0,0.01],[0.0,0.0,0.01]]
q_path_hex_edge = [[0.0, -0.5], [0.0, 0.0], [0.0, .5]]

q_spacing = 100
q_spacing_edge = 20


q_grid = ['slice',[51, 51, 1], 0.0]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
q_path_K = [[1.0/3+0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3-0.01]]
q_path_Kp = [[2.0/3+0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3-0.01]]
q_path_X = [[1.0/2+0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2-0.0001]]
q_grid_berry_K = ['berryphase', q_path_K, 50]
q_grid_berry_X = ['berryphase', q_path_X, 50]

###############################################################################
DM_hex_edge = DynamicalMatrix('hexagonal_edge_test', FC_hex_edge, alpha0)
#DM_hex_edge = DynamicalMatrix('hexagonal_edge_test', FC_hex_edge.dimension, FC_hex_edge.num_atom, FC_hex_edge.latt_vec,FC_hex_edge.atom_pos,FC_hex_edge.atom_mas,FC_hex_edge.fc_info,alpha0)

DM_hex_edge.get_phonon_band(q_path_hex_edge,q_spacing_edge)
#DM_hex_edge.draw_phonon_band()
DM_hex_edge.draw_phonon_band_edge_projection()