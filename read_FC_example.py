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
sys.path.insert(0, '/Users/woohyunhan/githup_reposi/phononTB')
from draw_phononTB_HAN import *
#from numeric import *


FC_from_phonopy = Read_FC_from_other_calculators('PHONOPY')
FC_symmetry = False
#supercell =  [2,2,2]
supercell = [4,4,4]
#mass = [95.94, 95.94, 78.96, 78.96, 78.96, 78.96]
#mass = [95.94, 78.96,  78.96]
#mass = [28.085,28.085,28.085,28.085,106.42,106.42,106.42,106.42]
mass = [28.085,28.085]
symprec = 1e-4
cutoff_distance = 10.0


FC_from_phonopy.read_FC(FC_symmetry, supercell, mass, symprec, cutoff_distance)
FC_from_phonopy.print_all_information()
#print FC_from_phonopy.latt_vec_prim


################################################################
FC = ForceConstant(3, FC_from_phonopy.total_num_atom_prim)
FC.get_fc_other_calculators(FC_from_phonopy, asr=True)
################################################################



################################################################

#q_path = [[0.000 ,  0.000 ,  0.000]  ,   [ 0.500 ,  0.000  , 0.000]    , [0.3333, 0.3333 , 0.000]    , [0.000,   0.000  , 0.000] ]
#q_path = [[0.000 ,  0.000 ,  0.000]  ,   [ 0.500 ,  0.000  , 0.000]    , [0.3333, 0.3333 , 0.000]    , [0.000,   0.000  , 0.000] ,    [0.000   ,0.000 ,  0.500]    , [0.500 ,  0.000   ,0.500]    , [0.3333 , 0.3333  ,0.500] ,    [0.000  , 0.000 ,  0.500]]
#q_path = [[0.0, 0.0 ,0.0] ,  [0.5 ,0.0 ,0.0] ,  [0.5 ,0.5 ,0.0] , [0.0 ,0.5 ,0.0]  , [0.0, 0.0 ,0.0] ,  [0.0, 0.0 ,0.5]  , [0.5 ,0.0 ,0.5]  , [0.5, 0.5 ,0.5] ,[0.0, 0.5, 0.5] ,[0.0 ,0.0 ,0.5]]
q_path = [[0.0 ,0.0 ,0.0] ,  [0.5, 0.0 ,0.5 ] , [0.5 ,0.25  ,0.75]   ,[0.375  ,0.375,  0.75] ,  [0.0 ,0.0 ,0.0] , [0.5,  0.5 ,0.5]  , [0.625 ,0.25 ,0.625] , [0.5 , 0.25 ,0.75] ,  [0.5 ,0.5 ,0.5]  ,[0.375 ,0.375, 0.75] ]
q_spacing = 20
alpha0 = [[0.0, 0.0, 0.0]]

DM = DynamicalMatrix('Read_FC_test', FC, alpha0)

DM.get_phonon_band(q_path,q_spacing)
DM.draw_phonon_band()