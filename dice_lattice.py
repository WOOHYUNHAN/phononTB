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




lat_dice = [[1.2333638807431999,   -2.1362489056675500,0.0], [1.2333638807431999,   2.1362489056675500,0.0], [0.0, 0.0, 12.0]]
orb_dice = [[0.0/3, 0.0/3, 0.5], [1.0/3, 2.0/3, 0.48], [2.0/3, 1.0/3, 0.52]]
mass_dice = [10.0, 12.0, 12.0]
fc_nn  =  [[-1.0,0.0],[0.0,-1.0]]
fc_nn1 = [[-22.19684176,   -7.699687431],[-7.699687431,   -13.306]]
fc_nn2 = [[-22.19684176 ,    7.699687431],[ 7.699687431,   -13.306]]
fc_nn3 = [[-8.860588573805234,    -0.000000000000000],[-0.000000000000000 ,  -26.642257254446765]]
fc_nnn = [[0.1,0.0],[0.0,0.3]]
V_info_dice = [-1.0, -0.1]

###############################################################################
FC_dice = ForceConstant(3, 3)

FC_dice.set_geometry(lat_dice, orb_dice, mass_dice)
FC_dice.set_hopping(0,1,[0.0,0.0,0.0],V_info_dice)
FC_dice.set_hopping(0,1,[-1.0,-1.0,0.0],V_info_dice)
FC_dice.set_hopping(0,1,[0.0,-1.0,0.0],V_info_dice)

FC_dice.set_hopping(0,2,[0.0,0.0,0.0],V_info_dice)
FC_dice.set_hopping(0,2,[-1.0,-1.0,0.0],V_info_dice)
FC_dice.set_hopping(0,2,[-1.0,0.0,0.0],V_info_dice)

#FC_dice.set_hopping(0,0,[0.0,0.0],V_info_hex2)
#FC_dice.set_hopping(1,1,[0.0,0.0],V_info_hex2)

#FC_dice.set_hopping(1,1,[-1.0,-1.0],V_info2_hex)
#FC_dice.set_hopping(1,1,[0.0,1.0],V_info2_hex)
#FC_dice.set_hopping(1,1,[1.0,0.0],V_info2_hex)
#FC_dice.set_hopping(0,0,[-1.0,0.0],V_info2_           hex)
#FC_dice.set_hopping(0,0,[1.0,1.0],V_info2_hex)
#FC_dice.set_hopping(0,0,[0.0,-1.0],V_info2_hex)

FC_dice.set_acoustic_sum_rule()
#FC_dice.set_acoustic_sum_rule()

FC_dice.print_info()


###############################################################################
alpha0 = [[0.0,0.0,0.0]]
alpha3 = [[0.0,0.0,0.],[0.0,0.0,0.0],[0.0,0.0,0.0]]
q_path_dice = [[0, 0, 0.0], [0.5, 0.0, 0.0], [1.0/3, 1.0/3, 0.0], [0.0, 0, 0.0]]

q_spacing = 400


q_grid = ['slice',[51, 51, 1], 0.0]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
q_path_K = [[1.0/3+0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3+0.01],[1.0/3-0.01,1.0/3-0.01],[1.0/3+0.01,1.0/3-0.01]]
q_path_Kp = [[2.0/3+0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3+0.01],[2.0/3-0.01,-1.0/3-0.01],[2.0/3+0.01,-1.0/3-0.01]]
q_path_X = [[1.0/2+0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2+0.0001],[1.0/2-0.0001,0.0/2-0.0001],[1.0/2+0.0001,0.0/2-0.0001]]
q_grid_berry_K = ['berryphase', q_path_K, 50]
q_grid_berry_X = ['berryphase', q_path_X, 50]
print FC_dice.recip_vec
q_grid = ['berrycurv', [-2.7, 2.7], [-1.6, 1.6], [2, 0.0], 201, 201]

###############################################################################
DM_dice = DynamicalMatrix('dice_test', FC_dice ,alpha0)


#DM_dice.get_phonon_band(q_path_dice,q_spacing)
DM_dice.draw_phonon_band()
#DM_dice.make_phband_PROCAR_format(q_grid)
#DM_dice.get_3Dplot_data(12)
#DM_dice.make_anime_file_for_vsim([0.0, 0.0, 0.0])


###############################################################################
band_range = [int(i) for i in range(9,16)] ; print '# of bands = ' + str(len(band_range)) + ' Detail bands = ' + str(band_range)
#CTI_dice = ComputeTopologicalInvariants('dice_test',band_range, q_grid)
#CTI_dice.get_Willsons_loop()
#CTI_hex.calculate_Berry_phase()

