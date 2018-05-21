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




lat_hex_bilayer = [[1.2333638807431999,   -2.1362489056675500, 0.0], [1.2333638807431999,   2.1362489056675500, 0.0], [0.0,0.0,15.0]]
orb_hex_bilayer = [[0.5000000000000000,  0.5000000000000000,  0.3919019252748578], [0.8333340000000007 , 0.1666659999999993,  0.3918741643673117], [0.1666669999999968,  0.8333330000000032 , 0.6081258356326885], [0.5000000000000000,  0.5000000000000000 , 0.6080963917251414]]
mass_hex_bilayer = [12.0, 12.0, 12.0, 12.0]
fc_nn  =  [[-1.0,0.0],[0.0,-1.0]]
fc_nn1 = [[-22.19684176,   -7.699687431],[-7.699687431,   -13.306]]
fc_nn2 = [[-22.19684176 ,    7.699687431],[ 7.699687431,   -13.306]]
fc_nn3 = [[-8.860588573805234,    -0.000000000000000],[-0.000000000000000 ,  -26.642257254446765]]
fc_nnn = [[0.1,0.0],[0.0,0.3]]
V_info_hex = [-1.0, -0.2]
V_info_hex_inter = np.array([-1.0, -1.0]) / 50.0

###############################################################################
FC_hex_bilayer = ForceConstant(3, 4)

FC_hex_bilayer.set_geometry(lat_hex_bilayer, orb_hex_bilayer, mass_hex_bilayer)

FC_hex_bilayer.set_hopping(0,1,[0.0,0.0,0.0],V_info_hex)
FC_hex_bilayer.set_hopping(0,1,[0.0,1.0,0.0],V_info_hex)
FC_hex_bilayer.set_hopping(0,1,[-1.0,0.0,0.0],V_info_hex)

FC_hex_bilayer.set_hopping(2,3,[0.0,0.0,0.0],V_info_hex)
FC_hex_bilayer.set_hopping(2,3,[0.0,1.0,0.0],V_info_hex)
FC_hex_bilayer.set_hopping(2,3,[-1.0,0.0,0.0],V_info_hex)

FC_hex_bilayer.set_hopping(0,2,[0.0,0.0,0.0],V_info_hex_inter)
FC_hex_bilayer.set_hopping(0,2,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_bilayer.set_hopping(0,2,[1.0,0.0,0.0],V_info_hex_inter)

FC_hex_bilayer.set_hopping(0,3,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_bilayer.set_hopping(1,2,[1.0,0.0,0.0],V_info_hex_inter)
FC_hex_bilayer.set_hopping(1,2,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_bilayer.set_hopping(1,2,[1.0,-1.0,0.0],V_info_hex_inter)

FC_hex_bilayer.set_hopping(1,3,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_bilayer.set_hopping(1,3,[1.0,0.0,0.0],V_info_hex_inter)
FC_hex_bilayer.set_hopping(1,3,[0.0,0.0,0.0],V_info_hex_inter)



FC_hex_bilayer.set_acoustic_sum_rule()



###############################################################################
alpha0 = [[0.0,0.0,0.0]]
alpha4 = [[0.0,0.0,0.01],[0.0,0.0,0.01],[0.0,0.0,-0.01],[0.0,0.0,-0.01]]
q_path_hex_bilayer = [[0, 0,0], [0.5, 0.0,0], [1.0/3, 1.0/3,0], [0.0, 0,0]]

q_spacing = 100


q_grid = ['slice',[51, 51, 1], 0.0]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
q_path_K = [[1.0/3+0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3-0.01]]
q_path_Kp = [[2.0/3+0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3-0.01]]
q_path_X = [[1.0/2+0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2-0.0001]]
q_grid_berry_K = ['berryphase', q_path_K, 50]
q_grid_berry_X = ['berryphase', q_path_X, 50]

###############################################################################
DM_hex_bilayer = DynamicalMatrix('hexagonal_bilayer', FC_hex_bilayer, alpha4)

DM_hex_bilayer.get_phonon_band(q_path_hex_bilayer,q_spacing)
DM_hex_bilayer.draw_phonon_band()
#DM_hex_bilayer.make_phband_PROCAR_format(q_grid)


###############################################################################
#band_range = [int(i) for i in range(4,6)] ; print '# of bands = ' + str(len(band_range)) + ' Detail bands = ' + str(band_range)
#CTI_hex = ComputeTopologicalInvariants('hexagonal_test',band_range, q_grid)
#CTI_hex.get_Willsons_loop()
#CTI_hex.calculate_Berry_phase()

