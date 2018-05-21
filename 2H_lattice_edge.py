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



lat = [[1.59129249942129 , -2.75619945870095,  0.00000000000000], [1.59129249942129 ,  2.75619945870095,  0.00000000000000], [0.00000000000000 ,  0.00000000000000, 10.00000000000000]]
orb = [[0.66666600000000 ,  0.33333300000000,  0.71006677700000 ], [0.00000000000000 ,  0.00000000000000,  0.86644056100000], [0.00000000000000 ,  0.00000000000000,  0.55369938800000]]
mass = [36.0, 12.0, 12.0]
fc_nn  =  [[-1.0,0.0],[0.0,-1.0]]
fc_nn1 = [[-22.19684176,   -7.699687431],[-7.699687431,   -13.306]]
fc_nn2 = [[-22.19684176 ,    7.699687431],[ 7.699687431,   -13.306]]
fc_nn3 = [[-8.860588573805234,    -0.000000000000000],[-0.000000000000000 ,  -26.642257254446765]]
fc_nnn = [[0.1,0.0],[0.0,0.3]]
V_info = [-1.0, -0.1]
num_repeat = 31
###############################################################################
FC_edge = ForceConstant(3, 3)

FC_edge.set_geometry(lat, orb, mass)
FC_edge.set_hopping(0,1,[0.0,0.0,0.0],V_info)
FC_edge.set_hopping(0,1,[1.0,0.0,0.0],V_info)
FC_edge.set_hopping(0,1,[1.0,1.0,0.0],V_info)

FC_edge.set_hopping(0,2,[0.0,0.0,0.0],V_info)
FC_edge.set_hopping(0,2,[1.0,0.0,0.0],V_info)
FC_edge.set_hopping(0,2,[1.0,1.0,0.0],V_info)
#FC_dice.set_hopping(0,0,[0.0,0.0],V_info_hex2)
#FC_dice.set_hopping(1,1,[0.0,0.0],V_info_hex2)

#FC_dice.set_hopping(1,1,[-1.0,-1.0],V_info2_hex)
#FC_dice.set_hopping(1,1,[0.0,1.0],V_info2_hex)
#FC_dice.set_hopping(1,1,[1.0,0.0],V_info2_hex)
#FC_dice.set_hopping(0,0,[-1.0,0.0],V_info2_           hex)
#FC_dice.set_hopping(0,0,[1.0,1.0],V_info2_hex)
#FC_dice.set_hopping(0,0,[0.0,-1.0],V_info2_hex)
FC_edge.make_edge(num_repeat,0)
FC_edge.set_acoustic_sum_rule()


#FC_dice.set_acoustic_sum_rule()

FC_edge.print_info()

###############################################################################
alpha0 = [[0.0,0.0, 0.0]]
alpha3 = [[0.0,0.0,0.04],[0.0,0.0,0.04],[0.0,0.0,0.04]]
alpha3_edge = []
for i in range(len(alpha3)):
	for j in range(num_repeat):
		alpha3_edge.append(alpha3[i])

q_path_edge = [[0.0, -0.5, 0.0], [0.0, 0.0, 0.0], [0.0, .5, 0.0]]

q_spacing = 100
q_spacing_edge = 40


q_grid = ['slice',[51, 51, 1], 0.0]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
q_path_K = [[1.0/3+0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3-0.01]]
q_path_Kp = [[2.0/3+0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3-0.01]]
q_path_X = [[1.0/2+0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2-0.0001]]
q_grid_berry_K = ['berryphase', q_path_K, 50]
q_grid_berry_X = ['berryphase', q_path_X, 50]

###############################################################################
DM_edge = DynamicalMatrix('2H_edge_test', FC_edge, alpha0)
#DM_hex_edge = DynamicalMatrix('hexagonal_edge_test', FC_hex_edge.dimension, FC_hex_edge.num_atom, FC_hex_edge.latt_vec,FC_hex_edge.atom_pos,FC_hex_edge.atom_mas,FC_hex_edge.fc_info,alpha0)

DM_edge.get_phonon_band(q_path_edge,q_spacing_edge)
#DM_hex_edge.draw_phonon_band()
DM_edge.draw_phonon_band_edge_projection()