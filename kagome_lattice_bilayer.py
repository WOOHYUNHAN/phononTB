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




lat = [[2.0,  0.0, 0.0], [-1.0,  np.sqrt(3), 0.0], [0.0, 0.0, 10.0]]
orb = [[0.5,0.0,0.4], [0.0,0.5,0.4], [0.5,0.5,0.4], [0.5,0.0,0.6], [0.0,0.5,0.6], [0.5,0.5,0.6]]
mass = [12.0, 12.0, 12.0 , 12.0, 12.0, 12.0]
fc_nn  =  [[-1.0,0.0],[0.0,-1.0]]
fc_nn1 = [[-22.19684176,   -7.699687431],[-7.699687431,   -13.306]]
fc_nn2 = [[-22.19684176 ,    7.699687431],[ 7.699687431,   -13.306]]
fc_nn3 = [[-8.860588573805234,    -0.000000000000000],[-0.000000000000000 ,  -26.642257254446765]]
fc_nnn = [[0.1,0.0],[0.0,0.3]]
V_info = [-1.0, -0.3]
V_info_inter = np.array([-1.0, -1.0]) / 50.0

###############################################################################
FC = ForceConstant(3, 6)

FC.set_geometry(lat, orb, mass)

FC.set_hopping(0,1,[1.0, 0.0, 0.0],V_info)
FC.set_hopping(0,1,[0.0, -1.0, 0.0],V_info)

FC.set_hopping(0,2,[0.0, 0.0, 0.0],V_info) 
FC.set_hopping(0,2,[0.0, -1.0, 0.0],V_info) 
  
FC.set_hopping(1,2,[0.0, 0.0, 0.0],V_info) 
FC.set_hopping(1,2,[-1.0, 0.0, 0.0],V_info)

FC.set_hopping(3,4,[1.0, 0.0, 0.0],V_info)
FC.set_hopping(3,4,[0.0, -1.0, 0.0],V_info)

FC.set_hopping(3,5,[0.0, 0.0, 0.0],V_info) 
FC.set_hopping(3,5,[0.0, -1.0, 0.0],V_info) 
  
FC.set_hopping(4,5,[0.0, 0.0, 0.0],V_info) 
FC.set_hopping(4,5,[-1.0, 0.0, 0.0],V_info) 

FC.set_hopping(0,3,[0.0, 0.0, 0.0],V_info_inter) 
FC.set_hopping(1,4,[0.0, 0.0, 0.0],V_info_inter) 
FC.set_hopping(2,5,[0.0, 0.0, 0.0],V_info_inter) 

#
FC.set_acoustic_sum_rule()

#FC_hex.set_acoustic_sum_rule()

FC.print_info()

###############################################################################
alpha0 = [[0.0,0.0,0.04]]
alpha6 = [[0.0,0.0,0.05],[0.0,0.0,0.05],[0.0,0.0,0.05],[0.0,0.0,-0.05],[0.0,0.0,-0.05],[0.0,0.0,-0.05]]
q_path_hex_bilayer = [[0, 0,0], [0.5, 0.0,0], [1.0/3, 1.0/3,0], [0.0, 0,0]]

q_spacing = 100


q_grid = ['slice',[51, 51, 1], 0.0]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
q_path_K = [[1.0/3+0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3-0.01]]
q_path_Kp = [[2.0/3+0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3-0.01]]
q_path_X = [[1.0/2+0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2-0.0001]]
q_grid_berry_K = ['berryphase', q_path_K, 50]
q_grid_berry_X = ['berryphase', q_path_X, 50]

###############################################################################
DM = DynamicalMatrix('kagome_bilayer', FC, alpha6)

DM.get_phonon_band(q_path_hex_bilayer,q_spacing)
DM.draw_phonon_band()
#DM_hex_bilayer.make_phband_PROCAR_format(q_grid)


###############################################################################
#band_range = [int(i) for i in range(4,6)] ; print '# of bands = ' + str(len(band_range)) + ' Detail bands = ' + str(band_range)
#CTI_hex = ComputeTopologicalInvariants('hexagonal_test',band_range, q_grid)
#CTI_hex.get_Willsons_loop()
#CTI_hex.calculate_Berry_phase()

