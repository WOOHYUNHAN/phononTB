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
import sys
sys.path.insert(0, 'D:\github_localrep\phononTB')
from draw_phononTB_HAN import *
#from numeric import *




lat_hex = [[1.2333638807431999,   -2.1362489056675500], [1.2333638807431999,   2.1362489056675500]]
orb_hex = [[1.0/3, 2.0/3], [2.0/3, 1.0/3]]
mass_hex = [12.0, 12.0]
fc_nn  =  [[-1.0,0.0],[0.0,-1.0]]
fc_nn1 = [[-22.19684176,   -7.699687431],[-7.699687431,   -13.306]]
fc_nn2 = [[-22.19684176 ,    7.699687431],[ 7.699687431,   -13.306]]
fc_nn3 = [[-8.860588573805234,    -0.000000000000000],[-0.000000000000000 ,  -26.642257254446765]]
fc_nnn = [[0.1,0.0],[0.0,0.3]]
V_info_hex = [-1.0, -0.1]
V_info_hex2 = [0.0,1.6499999999999977]

###############################################################################
FC_hex = ForceConstant(2, 2)

FC_hex.set_geometry(lat_hex, orb_hex, mass_hex)
FC_hex.set_hopping(0,1,[0.0,0.0],V_info_hex)
FC_hex.set_hopping(0,1,[-1.0,0.0],V_info_hex)
FC_hex.set_hopping(0,1,[0.0,1.0],V_info_hex)
#FC_hex.set_hopping(0,0,[0.0,0.0],V_info_hex2)
#FC_hex.set_hopping(1,1,[0.0,0.0],V_info_hex2)

#FC_hex.set_hopping(1,1,[-1.0,-1.0],V_info2_hex)
#FC_hex.set_hopping(1,1,[0.0,1.0],V_info2_hex)
#FC_hex.set_hopping(1,1,[1.0,0.0],V_info2_hex)
#FC_hex.set_hopping(0,0,[-1.0,0.0],V_info2_           hex)
#FC_hex.set_hopping(0,0,[1.0,1.0],V_info2_hex)
#FC_hex.set_hopping(0,0,[0.0,-1.0],V_info2_hex)

FC_hex.set_acoustic_sum_rule()
#FC_hex.set_acoustic_sum_rule()

FC_hex.print_info()


###############################################################################
alpha0 = [[0.0,0.0,0.04]]
alpha2 = [[0.0,0.0,0.01],[0.0,0.0,0.01]]
q_path_hex = [[0, 0], [0.5, 0.0], [1.0/3, 1.0/3], [0.0, 0]]

q_spacing = 100


q_grid = ['slice',[51, 51, 1], 0.0]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
q_path_K = [[1.0/3+0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3-0.01]]
q_path_Kp = [[2.0/3+0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3-0.01]]
q_path_X = [[1.0/2+0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2-0.0001]]
q_grid_berry_K = ['berryphase', q_path_K, 50]
q_grid_berry_X = ['berryphase', q_path_X, 50]

###############################################################################
DM_hex = DynamicalMatrix('hexagonal_test', FC_hex ,alpha0)


#DM_hex.get_phonon_band(q_path_hex,q_spacing)
#DM_hex.draw_phonon_band()
#DM_hex.make_phband_PROCAR_format(q_grid)

###############################################################################
band_range = [int(i) for i in range(6,8)] ; print '# of bands = ' + str(len(band_range)) + ' Detail bands = ' + str(band_range)
CTI_hex = ComputeTopologicalInvariants('hexagonal_test',band_range, q_grid)
CTI_hex.get_Willsons_loop()
#CTI_hex.calculate_Berry_phase()

