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




lat_hex_hexalayer = [[1.2333638807431999,   -2.1362489056675500, 0.0], [1.2333638807431999,   2.1362489056675500, 0.0], [0.0,0.0, 40.0]]
orb_hex_hexalayer = [[0.5,	0.5,	0.146963222], [0.833334	,0.166666,	0.146952812], [0.166667,	0.833333	,0.228047188], [0.5	,0.5	,0.228036147],[0.5	,0.5	,0.309130524],[0.833334,	0.166666,	0.309119482],[0.166667	,0.833333	,0.390213859],[0.5,	0.5	,0.390202817],[0.5,	0.5,	0.471297194], [0.833334,	0.166666,	0.471286153],[0.166667	,0.833333	,0.552380529],[0.5	,0.5	,0.552369488]]
mass_hex_hexalayer = [12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0]
fc_nn  =  [[-1.0,0.0],[0.0,-1.0]]
fc_nn1 = [[-22.19684176,   -7.699687431],[-7.699687431,   -13.306]]
fc_nn2 = [[-22.19684176 ,    7.699687431],[ 7.699687431,   -13.306]]
fc_nn3 = [[-8.860588573805234,    -0.000000000000000],[-0.000000000000000 ,  -26.642257254446765]]
fc_nnn = [[0.1,0.0],[0.0,0.3]]
V_info_hex = [-1.0, -0.2]
V_info_hex_inter = np.array([-1.0, -0.2]) / 40.0

###############################################################################
FC_hex_hexalayer = ForceConstant(3, 12)

FC_hex_hexalayer.set_geometry(lat_hex_hexalayer, orb_hex_hexalayer, mass_hex_hexalayer)

FC_hex_hexalayer.set_hopping(0,1,[0.0,0.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(0,1,[0.0,1.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(0,1,[-1.0,0.0,0.0],V_info_hex)

FC_hex_hexalayer.set_hopping(2,3,[0.0,0.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(2,3,[0.0,1.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(2,3,[-1.0,0.0,0.0],V_info_hex)

FC_hex_hexalayer.set_hopping(4,5,[0.0,0.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(4,5,[0.0,1.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(4,5,[-1.0,0.0,0.0],V_info_hex)

FC_hex_hexalayer.set_hopping(6,7,[0.0,0.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(6,7,[0.0,1.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(6,7,[-1.0,0.0,0.0],V_info_hex)

FC_hex_hexalayer.set_hopping(8,9,[0.0,0.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(8,9,[0.0,1.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(8,9,[-1.0,0.0,0.0],V_info_hex)

FC_hex_hexalayer.set_hopping(10,11,[0.0,0.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(10,11,[0.0,1.0,0.0],V_info_hex)
FC_hex_hexalayer.set_hopping(10,11,[-1.0,0.0,0.0],V_info_hex)



FC_hex_hexalayer.set_hopping(0,2,[0.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(0,2,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(0,2,[1.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(1,2,[1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(1,2,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(1,2,[1.0,-1.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(1,3,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(1,3,[1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(1,3,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(0,3,[0.0,0.0,0.0],V_info_hex_inter)


FC_hex_hexalayer.set_hopping(4,2,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(4,2,[1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(4,2,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(4,3,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(5,2,[1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(5,2,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(5,2,[1.0,-1.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(5,3,[0.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(5,3,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(5,3,[1.0,0.0,0.0],V_info_hex_inter)


FC_hex_hexalayer.set_hopping(6,4,[0.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(6,4,[0.0,1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(6,4,[-1.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(6,5,[0.0,1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(6,5,[-1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(6,5,[-1.0,1.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(7,5,[0.0,1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(7,5,[-1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(7,5,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(7,4,[0.0,0.0,0.0],V_info_hex_inter)


FC_hex_hexalayer.set_hopping(8,6,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(8,6,[1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(8,6,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(8,7,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(9,6,[1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(9,6,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(9,6,[1.0,-1.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(9,7,[0.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(9,7,[0.0,-1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(9,7,[1.0,0.0,0.0],V_info_hex_inter)


FC_hex_hexalayer.set_hopping(10,8,[0.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(10,8,[0.0,1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(10,8,[-1.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(10,9,[0.0,1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(10,9,[-1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(10,9,[-1.0,1.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(11,9,[0.0,1.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(11,9,[-1.0,0.0,0.0],V_info_hex_inter)
FC_hex_hexalayer.set_hopping(11,9,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_hopping(11,8,[0.0,0.0,0.0],V_info_hex_inter)

FC_hex_hexalayer.set_acoustic_sum_rule()



###############################################################################
alpha0 = [[0.0,0.0,0.0085]]
alpha4 = [[0.0,0.0,0.01],[0.0,0.0,0.01],[0.0,0.0,-0.01],[0.0,0.0,-0.01]]
q_path_hex_hexalayer = [[1.0/3, 1.0/3,0], [0.0, 0,0], [-1.0/3, -1.0/3, 0]]

q_spacing = 100


q_grid = ['slice',[51, 51, 1], 0.0]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
q_path_K = [[1.0/3+0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3-0.01]]
q_path_Kp = [[2.0/3+0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3-0.01]]
q_path_X = [[1.0/2+0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2-0.0001]]
q_grid_berry_K = ['berryphase', q_path_K, 50]
q_grid_berry_X = ['berryphase', q_path_X, 50]

###############################################################################
DM_hex_hexalayer = DynamicalMatrix('hexagonal_hexalayer', FC_hex_hexalayer, alpha0)

DM_hex_hexalayer.get_phonon_band(q_path_hex_hexalayer,q_spacing)
DM_hex_hexalayer.draw_phonon_band()
#DM_hex_hexalayer.make_phband_PROCAR_format(q_grid)


###############################################################################
band_range = [int(i) for i in range(30,50)] ; print '# of bands = ' + str(len(band_range)) + ' Detail bands = ' + str(band_range)
#CTI_hex = ComputeTopologicalInvariants('hexagonal_hexalayer',band_range, q_grid)
#CTI_hex.get_Willsons_loop()
#CTI_hex.calculate_Berry_phase()

