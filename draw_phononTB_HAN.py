import os.path
import time
import sys
from math import *
from math import sqrt
from numpy.linalg import *
import numpy as np
import matplotlib.pyplot as plt
import cmath
import datetime as dt
#from numeric import *

class Read_FC_from_other_calculators:
    def __init__(self, calculator):
        self.fc_information = []
        self.hopping = []
        self.calculator = calculator # phononpy and DFPT in QE
        if self.calculator == "PHONOPY" or self.calculator == "phonopy" or self.calculator == "Phonopy":
            self.version = '1.12.6.53' # latest version, FORCE_CONSTANTS file depends on the version so please check
        if self.calculator == "QE" or self.calculator == "Quantum Espresso" or self.calculator == "Quantum":
            self.version = 'Not yet'
            print "Not implemented yet"
            return 0
    
    def read_FC(self, FC_symmetry, supercell, mass, symprec, cutoff_distance):
        if self.calculator == "PHONOPY" or self.calculator == "phonopy" or self.calculator == "Phonopy":
            self.FC_symmetry = FC_symmetry ### True or False
            self.supercell = supercell ### [2, 2, 2]
            self.supercell_index = []
            for i in range(self.supercell[2]):
                for j in range(self.supercell[1]):
                    for k in range(self.supercell[0]):
                        self.supercell_index.append([k,j,i])
            self.mass = mass
            self.symprec = symprec ### default = 1e-4
            self.cutoff_distance = cutoff_distance
            self.read_FC_from_phonopy()
        else:
            # To be implemented
            return 0

    def read_POSCAR(self, POSCAR_file):
        f = open(POSCAR_file, 'r')
        tempf = f.readlines()
        f.close()

        univ_prim = float(tempf[1].split()[0])
        a_prim = np.array([float(tempf[2].split()[i]) for i in range(3)]) * univ_prim
        b_prim = np.array([float(tempf[3].split()[i]) for i in range(3)]) * univ_prim
        c_prim = np.array([float(tempf[4].split()[i]) for i in range(3)]) * univ_prim


        atom_type = np.array([str(tempf[5].split()[i]) for i in range(len(tempf[5].split()))])
        num_type = np.array([int(tempf[6].split()[i]) for i in range(len(tempf[6].split()))])
        total_num_atom = np.sum(num_type)

        #print total_num_atom

        if tempf[7][0] == 'D' or tempf[7][0] == 'd':
            # Direct coordinates
            mode = 'direct'
        elif tempf[7][0] == 'C' or tempf[7][0] == 'c':
            # Cartesian coordinates
            mode = 'cart'
        else:
            # Starnge input?
            print "POSCAR ERROR: please check direct or cartesian coordintes"
            return 0

        direct_coord = []

        if mode == 'direct':
            for i in range(total_num_atom):
                tempcoord = np.array([float(tempf[8+i].split()[j]) for j in range(3)])
                direct_coord.append(tempcoord)
        elif mode == 'cart':
            for i in range(total_num_atom):
                tempcoord = np.array([float(tempf[8+i].split()[j]) for j in range(3)])
                temp_direct_coord = from_cartesin_to_direct([a_prim, b_prim, c_prim], tempcoord)
                direct_coord.append(temp_direct_coord)
        else:
            pass

        self.latt_vec_prim = np.array([a_prim, b_prim, c_prim])
        self.atom_type_prim = atom_type
        self.num_type_prim = num_type
        self.total_num_atom_prim = total_num_atom
        self.direct_coord_prim = direct_coord

    def read_SPOSCAR(self, SPOSCAR_file):
        f = open(SPOSCAR_file, 'r')
        tempf = f.readlines()
        f.close()

        univ_super = float(tempf[1].split()[0])
        a_super = np.array([float(tempf[2].split()[i]) for i in range(3)]) * univ_super
        b_super = np.array([float(tempf[3].split()[i]) for i in range(3)]) * univ_super
        c_super = np.array([float(tempf[4].split()[i]) for i in range(3)]) * univ_super


        num_type = np.array([int(tempf[5].split()[i]) for i in range(len(tempf[5].split()))])
        total_num_atom = np.sum(num_type)

        #print total_num_atom

        if tempf[6][0] == 'D' or tempf[6][0] == 'd':
            # Direct coordinates
            mode = 'direct'
        elif tempf[6][0] == 'C' or tempf[6][0] == 'c':
            # Cartesian coordinates
            mode = 'cart'
        else:
            # Starnge input?
            print "POSCAR ERROR: please check direct or cartesian coordintes"
            return 0

        direct_coord = []

        if mode == 'direct':
            for i in range(total_num_atom):
                tempcoord = np.array([float(tempf[7+i].split()[j]) for j in range(3)])
                direct_coord.append(tempcoord)
        elif mode == 'cart':
            for i in range(total_num_atom):
                tempcoord = np.array([float(tempf[7+i].split()[j]) for j in range(3)])
                temp_direct_coord = from_cartesin_to_direct([a_super, b_super, c_super], tempcoord)
                direct_coord.append(temp_direct_coord)
        else:
            pass

        self.latt_vec_super = np.array([a_super, b_super, c_super])
        self.total_num_atom_super = total_num_atom
        self.direct_coord_super = direct_coord



    def read_FC_from_phonopy(self, POSCAR='POSCAR', SPOSCAR='SPOSCAR', FORCE_CONSTANTS='FORCE_CONSTANTS'):
        self.read_POSCAR(POSCAR)
        self.read_SPOSCAR(SPOSCAR)

        if len(self.mass) != self.total_num_atom_prim:
            print "ERROR: size error of mass or POSCAR error "
            return 0

        f = open(FORCE_CONSTANTS, 'r')
        tempf = f.readlines()
        f.close()

        total_line = 1 + self.total_num_atom_prim * (4 * self.total_num_atom_super)


        if total_line != len(tempf):
            print "ERROR: Something wrong, please check your phonopy version or FC, POSCAR, and SPOSCAR files"
            print "Version of phonopy we support is " + self.version
            return 0
        
        for i in range(self.total_num_atom_prim):
            for j in range(self.total_num_atom_super):
                templine = 1 + i*(4*self.total_num_atom_super) + (4*j)
                #print tempf[templine].split()
                a_atom, b_atom = int(tempf[templine].split()[0]) -1 , int(tempf[templine].split()[1]) -1 # 0, 1, 2, 3, 4, 5, ... atoms
                a_prim_index, a_super_pos = self.find_supercell_position(a_atom)
                b_prim_index, b_super_pos = self.find_supercell_position(b_atom)
                distance, b_min_super_pos, multi = self.find_nearest_supercell(a_prim_index, a_super_pos, b_prim_index, b_super_pos)
                #print a_prim_index, a_super_pos, b_prim_index, b_super_pos, distance, multi, b_min_super_pos
                if multi != 1:
                    #print distance, multi
                    FC_temp = np.zeros((3, 3))
                    FC_temp[0][0] = float(tempf[templine+1].split()[0]) ; FC_temp[0][1] = float(tempf[templine+1].split()[1]) ; FC_temp[0][2] = float(tempf[templine+1].split()[2])
                    FC_temp[1][0] = float(tempf[templine+2].split()[0]) ; FC_temp[1][1] = float(tempf[templine+2].split()[1]) ; FC_temp[1][2] = float(tempf[templine+2].split()[2])
                    FC_temp[2][0] = float(tempf[templine+3].split()[0]) ; FC_temp[2][1] = float(tempf[templine+3].split()[1]) ; FC_temp[2][2] = float(tempf[templine+3].split()[2])
                    #print a_prim_index, np.array(a_super_pos), a_prim_index, np.array(a_super_pos), distance, FC_temp
                    for k in range(multi): # This case: multi != 1
                        dup_check = self.check_duplicate(distance, a_prim_index, b_prim_index, b_min_super_pos[k])
                        if dup_check and distance <= self.cutoff_distance:
                            self.fc_information.append([a_prim_index, np.array(a_super_pos), b_prim_index, np.array(b_min_super_pos[k]), distance, FC_temp/multi])
                            ji_vector = np.dot(self.latt_vec_prim.transpose(), np.array(self.direct_coord_prim[b_prim_index]) + b_min_super_pos[k] - np.array(self.direct_coord_prim[a_prim_index]) - a_super_pos)
                            V_pps, V_ppp = self.find_sigma_pi_from_fc(FC_temp, ji_vector)
                            self.hopping.append([a_prim_index, b_prim_index, np.array(b_min_super_pos[k]), [V_pps, V_ppp]])
                else:
                    for k in range(multi): # This case: multi = 1
                        dup_check = self.check_duplicate(distance, a_prim_index, b_prim_index, b_min_super_pos[k])
                        if dup_check and distance <= self.cutoff_distance:
                            FC_temp = np.zeros((3, 3))
                            FC_temp[0][0] = float(tempf[templine+1].split()[0]) ; FC_temp[0][1] = float(tempf[templine+1].split()[1]) ; FC_temp[0][2] = float(tempf[templine+1].split()[2])
                            FC_temp[1][0] = float(tempf[templine+2].split()[0]) ; FC_temp[1][1] = float(tempf[templine+2].split()[1]) ; FC_temp[1][2] = float(tempf[templine+2].split()[2])
                            FC_temp[2][0] = float(tempf[templine+3].split()[0]) ; FC_temp[2][1] = float(tempf[templine+3].split()[1]) ; FC_temp[2][2] = float(tempf[templine+3].split()[2])
                            self.fc_information.append([a_prim_index, np.array(a_super_pos), b_prim_index, np.array(b_min_super_pos[k]), distance, FC_temp/multi])
                            ji_vector = np.dot(self.latt_vec_prim.transpose(), np.array(self.direct_coord_prim[b_prim_index]) + b_min_super_pos[k] - np.array(self.direct_coord_prim[a_prim_index]) - a_super_pos)
                            V_pps, V_ppp = self.find_sigma_pi_from_fc(FC_temp, ji_vector)
                            self.hopping.append([a_prim_index, b_prim_index, np.array(b_min_super_pos[k]), [V_pps, V_ppp]])
        #print self.fc_information
        #self.fc_information = np.array(self.fc_information)
        #np.save('inoformation_file', self.fc_information)
        return 0

    def check_duplicate(self, distance, a_atom_index, b_atom_index, b_atom_super_pos):
        ### This function checks whether hopping paramters are duplicated or not
        dup_check = True
        for i in range(len(self.fc_information)):
            if distance == self.fc_information[i][4]:
                if a_atom_index == b_atom_index:
                    if a_atom_index == self.fc_information[i][0] and b_atom_index == self.fc_information[i][2]:
                        if np.all(np.array(b_atom_super_pos)  == self.fc_information[i][3]):
                            if np.all(np.array(b_atom_super_pos)  == np.array([0,0,0])):
                                dup_check = False
                        if np.all(np.array(b_atom_super_pos) == -1 * self.fc_information[i][3]):
                            dup_check = False    
                else:
                    if a_atom_index == self.fc_information[i][0] and b_atom_index == self.fc_information[i][2]:
                        if np.all(np.array(b_atom_super_pos) == self.fc_information[i][3]):
                            dup_check = False
                    if a_atom_index == self.fc_information[i][2] and b_atom_index == self.fc_information[i][0]:  
                        if np.all(np.array(b_atom_super_pos)  == self.fc_information[i][3]):
                            dup_check = False
                        if np.all(np.array(b_atom_super_pos) == -1 * self.fc_information[i][3]):
                            dup_check = False                       
        return dup_check
      

    def find_supercell_position(self, atom_index):
        prim_index = int(atom_index) / (self.supercell[0] * self.supercell[1] * self.supercell[2])
        super_pos = self.supercell_index[int(atom_index) % (self.supercell[0] * self.supercell[1] * self.supercell[2])]
        return prim_index , np.array(super_pos)

    def find_nearest_supercell(self, a_atom_index, a_atom_super_pos, b_atom_index, b_atom_super_pos):
        distance = []
        supercell_pos = []
        min_super_pos = []
        s_pos = self.direct_coord_prim[a_atom_index] + np.array(a_atom_super_pos)
        p_pos = self.direct_coord_prim[b_atom_index] + np.array(b_atom_super_pos)
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for k in (-1, 0, 1):
                    diff = s_pos  - p_pos - np.array([i,j,k])*self.supercell
                    vec = np.dot(self.latt_vec_prim.transpose(), diff)
                    distance.append(np.linalg.norm(vec))
                    supercell_pos.append(np.array([i,j,k])*self.supercell + np.array(b_atom_super_pos))
        
        distance = np.array(distance)
        minimum = np.min(distance)
        for i in range(27):
            if abs(minimum - distance[i]) < self.symprec:
                min_super_pos.append(supercell_pos[i])
        multi = len(min_super_pos)

        return round(minimum,5), np.array(min_super_pos), multi
    
    def find_sigma_pi_from_fc(self, fc, ji_vector):
        if np.linalg.norm(ji_vector) == 0:
            l = 0.0 ; m = 0.0 ; n = 0.0
            final_results = [0, fc[0][0]]
        else:
            l = ji_vector[0] / np.linalg.norm(ji_vector)
            m = ji_vector[1] / np.linalg.norm(ji_vector)
            n = ji_vector[2] / np.linalg.norm(ji_vector)
            # l*l + m*m + n*n = 1.0
            if (l == 0.0 and m == 0.0): # n = 1.0
                final_results = [fc[2][2], fc[0][0]]
            elif (round(l*l,6) == round(m*m,6)): #
                if (round(l*l,6) == round(n*n,6)):
                    fc_matrix = np.array([fc[0][0], fc[0][1]])
                    transform_matrix = np.array([[l*l, 1.0-l*l],[l*m, -1.0*l*m]])
                    #print l*l, m*m, n*n, transform_matrix
                    final_results = np.dot(np.linalg.inv(transform_matrix), fc_matrix)
                else:
                    fc_matrix = np.array([fc[0][0], fc[2][2]])
                    transform_matrix = np.array([[l*l, 1.0-l*l],[n*n, 1.0-n*n]])
                    #print l*l, m*m, n*n, transform_matrix
                    final_results = np.dot(np.linalg.inv(transform_matrix), fc_matrix)
            else:
                fc_matrix = np.array([fc[0][0], fc[1][1]])
                transform_matrix = np.array([[l*l, 1.0-l*l],[m*m, 1.0-m*m]])
                final_results = np.dot(np.linalg.inv(transform_matrix), fc_matrix)
        V_pps = final_results[0] ; V_ppp = final_results[1]
        return V_pps, V_ppp

    def print_all_information(self):
        filename = 'information_file'
        f = open(filename, 'w')
        initial_line =' < ' + '\t' + 'Aatom' + '\t' + ' |  H  | ' + '\t' +  'Batom' + '\t'  + ' [ ' + '\t' + 'superx' + '\t' + 'supery' + '\t' + 'superz' + '\t' + ' ]  = ' + '\t' + 'distance' + '\t' + '===>' + '\t' + 'V_pps' + '\t' + 'V_ppp'  + '\n' 
        f.write(initial_line)
        for i in range(len(self.fc_information)):
            templine = ' < ' + '\t' + str(self.fc_information[i][0]) + '\t' + ' |  H  | ' + '\t' +  str(self.fc_information[i][2]) + '\t'  + ' [ ' + '\t' + str(self.fc_information[i][3][0]) + '\t' + str(self.fc_information[i][3][1]) + '\t' + str(self.fc_information[i][3][2]) + '\t' + ' ]  = ' + '\t' + str(self.fc_information[i][4]) + '\t' + '===>' + '\t' + str(self.hopping[i][3][0]) + '\t' + str(self.hopping[i][3][1]) 
            for j in range(3):
                for k in range(3):
                    templine += '\t' + str(self.fc_information[i][5][j][k])
            templine += '\n'
            f.write(templine)
        f.close()



class ForceConstant:
    def __init__(self, dimension, num_atom):
        self.dimension = dimension
        self.num_atom = num_atom
        self.hopping_save = []
        self.fc_info = []
        self.edge_cal = False

    def set_geometry(self, lattice_vector, atom_position, atom_mass):
        self.latt_vec = np.array(lattice_vector)
        self.atom_pos = np.array(atom_position)
        self.atom_mas = np.array(atom_mass)
        if len(self.latt_vec) != self.dimension:
            print 'Please check either dimension or lattice vector; They are differnt to each other'
        if len(self.atom_pos) != self.num_atom:
            print 'Please check either num_atom or atom_position; They are differnt to each other'
        if len(self.atom_mas) != self.num_atom:
            print 'Please check either atom_mas or atom_position; They are differnt to each other'
        
        ####reciprocal lattice####

        if self.dimension == 3:
            recip_fac =  (2*np.pi) / np.dot(self.latt_vec[0], np.cross(self.latt_vec[1], self.latt_vec[2]))
            recip_vec = np.array([np.cross(self.latt_vec[1], self.latt_vec[2]), np.cross(self.latt_vec[2], self.latt_vec[0]), np.cross(self.latt_vec[0], self.latt_vec[1])]) * recip_fac
        if self.dimension == 2:
            #print np.cross(self.latt_vec[0], self.latt_vec[1])
            recip_fac =  (2*np.pi) / np.cross(self.latt_vec[0], self.latt_vec[1])
            recip_vec = np.array([[self.latt_vec[1][1],-self.latt_vec[0][1]],[-self.latt_vec[1][0],self.latt_vec[0][0]]]) * recip_fac
            recip_vec = recip_vec.transpose()
        self.recip_vec = recip_vec


    def set_fc_direct(self,iatom, jatom, supercell_of_jatom, fc):
        fc_temp = [iatom, jatom, supercell_of_jatom, fc]
        self.fc_info.append(fc_temp)

    def set_hopping(self, iatom, jatom, supercell_of_jatom, V_info):
        #if len(force_constant) != self.dimension or len(force_constant[0]) != self.dimension:
        #    print 'Please check either dimension or the size of force constant; They are different to each other'
        #    return 0
        if len(supercell_of_jatom) != self.dimension:
            print 'Please check either dimension or the size of supercell; They are different to each other'
            return 0
        if iatom == jatom:
            print 'iatom == jatom; Do you really need it? or Do you intend to consider NNN ?'

        if not self.check_duplicate_TBmode(iatom, jatom, supercell_of_jatom):
            print "There is duplication !!!"
            print iatom, jatom, supercell_of_jatom
            sys.exit(1)

        self.hopping_save.append([iatom, jatom, supercell_of_jatom, V_info])
        fc = self.generate_fc_constant_TB(iatom, jatom, supercell_of_jatom, V_info[0], V_info[1])

        fc_temp = [iatom, jatom, supercell_of_jatom, fc]

        #print fc_temp

        self.fc_info.append(fc_temp)

    def check_duplicate_TBmode(self, iatom, jatom, supercell_of_jatom):
        for i in range(len(self.fc_info)):
            if (iatom == self.fc_info[i][0] and jatom == self.fc_info[i][1]):
                if np.all(np.array(supercell_of_jatom) - self.fc_info[i][2] == 0):
                    return False
            if (jatom == self.fc_info[i][0] and iatom == self.fc_info[i][1]):
                if np.all(np.array(supercell_of_jatom) + self.fc_info[i][2] == 0):
                    return False
        return True

    def set_acoustic_sum_rule(self):
        '''
        It determines whether on-site force constants statisfy constraints or not
        '''
        if self.dimension == 2:
            for i in range(self.num_atom):
                fc_from_asr = np.zeros((self.dimension,self.dimension))
                for j in range(len(self.fc_info)):
                    if i == self.fc_info[j][0] or i == self.fc_info[j][1]:
                        if (i == self.fc_info[j][0] and i == self.fc_info[j][1]):
                            if np.linalg.norm(self.fc_info[j][2] - np.array([0.0,0.0])) != 0:
                                fc_from_asr -= np.array(2.0*self.fc_info[j][3])  ### This case for 0, 0 [1,0] like NNN hopping
                            else:
                                fc_from_asr -= np.array(self.fc_info[j][3])
                        else:
                            fc_from_asr -= np.array(self.fc_info[j][3])
                #print fc_from_asr
                fc_temp = [i, i, [0,0], fc_from_asr]
                self.fc_info.append(fc_temp)
                hopping_from_asr = [0, fc_from_asr[0][0]]
                hopping_temp = [i, i, [0,0], hopping_from_asr]
                self.hopping_save.append(hopping_temp)

        if self.dimension == 3:
            for i in range(self.num_atom):
                fc_from_asr = np.zeros((self.dimension,self.dimension))
                for j in range(len(self.fc_info)):
                    if i == self.fc_info[j][0] or i == self.fc_info[j][1]:
                        if (i == self.fc_info[j][0] and i == self.fc_info[j][1]):
                            if np.linalg.norm(self.fc_info[j][2] - np.array([0.0,0.0,0.0])) != 0:
                                fc_from_asr -= np.array(2.0*self.fc_info[j][3])  ### This case for 0, 0 [1,0] like NNN hopping why 2.0? 0, 0, [1,0] == 0, 0 [-1, 0]
                            else:
                                fc_from_asr -= np.array(self.fc_info[j][3])  ### This case for 0, 0 [0,0] like on-site energies
                                #print i, self.fc_info[j][0], self.fc_info[j][1], self.fc_info[j][2], np.array(self.fc_info[j][3])
                        else:
                            fc_from_asr -= np.array(self.fc_info[j][3])
                fc_temp = [i, i, [0,0,0], fc_from_asr]
                #print fc_temp
                self.fc_info.append(fc_temp)
                hopping_from_asr = [0, fc_from_asr[0][0]]
                hopping_temp = [i, i, [0,0,0], hopping_from_asr]
                #self.hopping_save.append(hopping_temp)
            
        return 0

    def generate_fc_constant_TB(self, iatom, jatom, supercell_of_jatom, V_pps, V_ppp):
        '''
        Detail hopping parameters : https://en.wikipedia.org/wiki/Tight_binding
        E_xx = l^2 V_pps + (1-l^2) V_ppp
        E_xy = lm (V_pps - V_ppp)
        E_xz = ln (V_pps- V_ppp)
        l, m, n = direction cosines
        '''
        if self.dimension == 2:
            position_i = np.dot(np.array(self.latt_vec).transpose(), np.array(self.atom_pos[iatom]).transpose())
            position_j = np.dot(np.array(self.latt_vec).transpose(), np.array(self.atom_pos[jatom] + np.array(supercell_of_jatom)).transpose())
            ji_vector = position_j - position_i
            if np.linalg.norm(ji_vector) == 0:
                l, m = 0.0, 0.0
            else:
                l = ji_vector[0] / np.linalg.norm(ji_vector)
                m = ji_vector[1] / np.linalg.norm(ji_vector)
            fc_xx = l*l*V_pps + (1.0 - l*l)*V_ppp
            fc_xy = l*m*(V_pps - V_ppp)
            fc_yy = m*m*V_pps + (1.0 - m*m)*V_ppp

            fc = np.array([[fc_xx, fc_xy], [fc_xy, fc_yy]])

        if self.dimension == 3:
            position_i = np.dot(np.array(self.latt_vec).transpose(), np.array(self.atom_pos[iatom]).transpose())
            position_j = np.dot(np.array(self.latt_vec).transpose(), np.array(self.atom_pos[jatom] + np.array(supercell_of_jatom)).transpose())
            ji_vector = position_j - position_i
            if np.linalg.norm(ji_vector) == 0:
                l, m, n = 0.0, 0.0, 0.0
            else:
                l = ji_vector[0] / np.linalg.norm(ji_vector)
                m = ji_vector[1] / np.linalg.norm(ji_vector)
                n = ji_vector[2] / np.linalg.norm(ji_vector)
            fc_xx = l*l*V_pps + (1.0 - l*l)*V_ppp
            fc_xy = l*m*(V_pps - V_ppp)
            fc_yy = m*m*V_pps + (1.0 - m*m)*V_ppp
            fc_zz = n*n*V_pps + (1.0 - n*n)*V_ppp
            fc_xz = l*n*(V_pps - V_ppp)
            fc_yz = m*n*(V_pps - V_ppp)

            fc = np.array([[fc_xx, fc_xy, fc_xz], [fc_xy, fc_yy, fc_yz], [fc_xz, fc_yz, fc_zz]])

        return fc

    def get_fc_other_calculators(self, other_calculator, asr):
        self.set_geometry(other_calculator.latt_vec_prim, other_calculator.direct_coord_prim, other_calculator.mass)
        for i in range(len(other_calculator.fc_information)):
            self.set_fc_direct(other_calculator.fc_information[i][0], other_calculator.fc_information[i][2], other_calculator.fc_information[i][3], other_calculator.fc_information[i][5])
        #for i in range(len(other_calculator.hopping)):
        #    self.set_hopping(other_calculator.hopping[i][0], other_calculator.hopping[i][1], other_calculator.hopping[i][2], other_calculator.hopping[i][3])
        if asr:
            self.set_acoustic_sum_rule()

    def print_info(self):
        print 'Dimension = ' + str(self.dimension)
        print 'Number of atoms = ' + str(self.num_atom)
        if self.dimension == 2:
            print 'Lattice vector = ' + ' a = ' + str(self.latt_vec[0]) + ' b = ' + str(self.latt_vec[1])
        if self.dimension ==3:
            print 'Lattice vector = ' + ' a = ' + str(self.latt_vec[0]) + ' b = ' + str(self.latt_vec[1]) + ' c = ' + str(self.latt_vec[2])
        print 'Atomic position'
        for i in range(self.num_atom):
            line = ''
            for j in range(self.dimension):
                line += str(self.atom_pos[i][j]) + '\t'
            print line
        print 'All fc_info'
        for i in range(len(self.fc_info)):
            print str(self.fc_info[i])
        return 0

    def make_edge(self, num_repeat, direction):
        if direction == 0:
            print 'Edge direction is x '
        elif direction == 1:
            print 'Edge direction is y '
        elif direction == 2:
            print 'Edge direction is z '

        if direction + 1 > self.dimension:
            print 'ERROR: direction does not match with diemsion'
            return 0

        self.edge_cal = True
        self.edge_direction = direction

        multiple = np.zeros(self.dimension, dtype=float)
        for i in range(self.dimension):
            if i == self.edge_direction:
                multiple[i] = num_repeat
            else:
                multiple[i] = 1.0
        
        new_lattice = np.array([self.latt_vec[i]*multiple[i] for i in range(self.dimension)])
        new_atompos = []
        new_atommas = []

        atomspos_modified = self.atom_pos / multiple
        add = np.zeros(self.dimension, dtype=float)
        add[self.edge_direction] = 1.0/ num_repeat

        for i in range(self.num_atom):
            for j in range(num_repeat):
                new_atompos.append(atomspos_modified[i] + add*j)
                new_atommas.append(self.atom_mas[i])
        
        ##########

        for i in range(len(new_lattice)):
            if i == self.edge_direction:
                new_lattice[i]  = np.array(new_lattice[i]) * 1.0
        for i in range(len(new_atompos)):
            for j in range(self.dimension):
                if j == self.edge_direction:
                    new_atompos[i][j] = new_atompos[i][j] / 1.0


        ##########



        self.latt_vec = new_lattice
        self.atom_pos = new_atompos
        self.atom_mas = new_atommas
        self.num_atom = self.num_atom * num_repeat


        new_fc_info = []

        #for i in range(num_repeat):
        #    upper_limit = num_repeat - (i+1)
        #    lower_limit = -i
        #    #print upper_limit, lower_limit
        #    for j in range(len(self.hopping_save)):
        #        if self.hopping_save[j][2][self.edge_direction] > upper_limit or self.hopping_save[j][2][self.edge_direction] < lower_limit:
        #            pass
        #        else:
        #            reference_atom = int(self.hopping_save[j][0]*num_repeat + i)
        #            target_atom = int(self.hopping_save[j][1]*num_repeat + i + self.hopping_save[j][2][self.edge_direction])
        #            new_direction = [self.hopping_save[j][2][k] for k in range(self.dimension)] 
        #            new_direction[self.edge_direction] = 0.0
        #            new_fc = self.generate_fc_constant_TB(reference_atom, target_atom, new_direction, self.hopping_save[j][3][0], self.hopping_save[j][3][1])
        #            temp_new_fc_info = [reference_atom, target_atom, new_direction, new_fc]
        #            new_fc_info.append(temp_new_fc_info)
        ##print new_fc_info


        for i in range(num_repeat):
            upper_limit = num_repeat - (i+1)
            lower_limit = -i
            #print upper_limit, lower_limit
            for j in range(len(self.fc_info)):
                if self.fc_info[j][2][self.edge_direction] > upper_limit or self.fc_info[j][2][self.edge_direction] < lower_limit:
                    pass
                else:
                    reference_atom = int(self.fc_info[j][0]*num_repeat + i)
                    target_atom = int(self.fc_info[j][1]*num_repeat + i + self.fc_info[j][2][self.edge_direction])
                    new_direction = [self.fc_info[j][2][k] for k in range(self.dimension)] 
                    new_direction[self.edge_direction] = 0.0
                    new_fc = self.fc_info[j][3]
                    temp_new_fc_info = [reference_atom, target_atom, new_direction, new_fc]
                    new_fc_info.append(temp_new_fc_info)

        self.fc_info = new_fc_info
        self.num_repeat = num_repeat
        


        #for i in range

        return 0

class DynamicalMatrix:
    def __init__(self, name, force_constant, alpha):
        self.out_tag = name
        self.dimension = force_constant.dimension
        self.num_atom = force_constant.num_atom
        self.latt_vec = np.array(force_constant.latt_vec)
        self.atom_pos = force_constant.atom_pos
        self.atom_mas = force_constant.atom_mas
        self.fc_info = force_constant.fc_info
        self.TRS_broken = alpha
        self.edge_cal = force_constant.edge_cal
        self.recip_vec = force_constant.recip_vec
        if self.edge_cal:
            self.num_repeat = force_constant.num_repeat

    def obtain_qpath(self, q_path, q_spacing):
        q_vec_list = []
        q_distance_list = []
        special_q_distance_list = [0.0]
        sq_distance = 0.0
        #print np.dot(recip_vec, self.latt_vec)
        for i in range(len(q_path)-1):
            initial = np.array(q_path[i])
            final = np.array(q_path[i+1])
            #print np.dot(recip_vec.transpose(), (final-initial).transpose())
            temp_sq_distance = np.linalg.norm(np.dot(self.recip_vec.transpose(), (final-initial).transpose()))
            sq_distance += temp_sq_distance
            special_q_distance_list.append(sq_distance)
            q_distance = special_q_distance_list[i]
            for j in range(q_spacing):
                delta = (final - initial) / float(q_spacing-1)
                temp = initial + delta*j
                q_vec_list.append(temp)
                temp_q_distance = np.linalg.norm(np.dot(self.recip_vec.transpose(), np.array(delta*j).transpose()))
                #print temp_q_distance
                q_distance_list.append(q_distance+temp_q_distance)

        return special_q_distance_list, q_distance_list, q_vec_list

    def find_fc_for_pair(self, iatom, jatom):
        proper_fc_info = []  # [r_j-r_i, fc]
        if iatom == jatom:
            for i in range(len(self.fc_info)):
                if (iatom == self.fc_info[i][0] and iatom == self.fc_info[i][1]):
                    #diff_r = np.array(self.atom_pos[jatom]) + np.dot(np.array(self.latt_vec).transpose(), np.array(self.fc_info[i][2]).transpose()) - np.array(self.atom_pos[iatom])
                    diff_r = np.array(self.atom_pos[jatom]) + np.array(self.fc_info[i][2]).transpose() - np.array(self.atom_pos[iatom])
                    proper_fc_info.append([diff_r, self.fc_info[i][3]])
        else:
            for i in range(len(self.fc_info)):
                if (iatom == self.fc_info[i][0] and jatom == self.fc_info[i][1]):
                    #diff_r = np.array(self.atom_pos[jatom]) + np.dot(np.array(self.latt_vec).transpose(), np.array(self.fc_info[i][2]).transpose()) - np.array(self.atom_pos[iatom])
                    diff_r = np.array(self.atom_pos[jatom]) + np.array(self.fc_info[i][2]).transpose() - np.array(self.atom_pos[iatom])
                    proper_fc_info.append([diff_r, self.fc_info[i][3]])
                if (jatom == self.fc_info[i][0] and iatom == self.fc_info[i][1]):
                    #diff_r = np.array(self.atom_pos[jatom])  - np.array(self.atom_pos[iatom]) - np.dot(np.array(self.latt_vec).transpose(), np.array(self.fc_info[i][2]).transpose())
                    diff_r = np.array(self.atom_pos[jatom])  - np.array(self.atom_pos[iatom]) -  np.array(self.fc_info[i][2]).transpose()
                    proper_fc_info.append([diff_r, self.fc_info[i][3]])
        return proper_fc_info



    def construct_dynamicalmatrix_q(self, q_vec):
        #recip_vec = self.reciprocal()
        q = np.dot(self.recip_vec.transpose(), np.array(q_vec).transpose())
        #q = np.array(q_vec)
        dm = np.zeros((self.dimension * self.num_atom, self.dimension * self.num_atom), dtype=complex)
        for i in range(self.num_atom):
            for j in range(self.num_atom):
                if i < j:
                    mass_factor = np.sqrt(float(self.atom_mas[i]) * float(self.atom_mas[j]))
                    dm_local = np.zeros((self.dimension, self.dimension), dtype=complex)
                    proper_fc_info = self.find_fc_for_pair(i,j)
                    for k in range(len(proper_fc_info)):
                        #r = np.array(proper_fc_info[k][0])
                        r = np.dot(self.latt_vec.transpose(), np.array(proper_fc_info[k][0]).transpose())
                        #print r
                        phase_factor = np.exp(np.vdot(r, q) * 1j)
                        dm_local += proper_fc_info[k][1] * phase_factor / mass_factor
                    dm[(i*self.dimension):(i*self.dimension+self.dimension), (j*self.dimension):(j*self.dimension+self.dimension)] += dm_local
                    dm[(j*self.dimension):(j*self.dimension+self.dimension), (i*self.dimension):(i*self.dimension+self.dimension)] += dm_local.conj().transpose()
                if i == j:
                    mass_factor = np.sqrt(float(self.atom_mas[i]) * float(self.atom_mas[j]))
                    dm_local = np.zeros((self.dimension, self.dimension), dtype=complex)
                    proper_fc_info = self.find_fc_for_pair(i,j)
                    for k in range(len(proper_fc_info)):
                        r = np.dot(self.latt_vec.transpose(), np.array(proper_fc_info[k][0]).transpose())
                        phase_factor = np.exp(np.vdot(r, q) * 1j)
                        if np.linalg.norm(r) == 0.0:
                            dm_local += proper_fc_info[k][1] * phase_factor / mass_factor
                        else:
                            dm_local += proper_fc_info[k][1] * phase_factor / mass_factor + (proper_fc_info[k][1] * phase_factor / mass_factor).conj().transpose()
                    dm[(i*self.dimension):(i*self.dimension+self.dimension), (j*self.dimension):(j*self.dimension+self.dimension)] += dm_local
        
        dynamical_matrix = (dm + dm.conj().transpose()) / 2.0 
        #dynamical_matrix = dm

        return dynamical_matrix

    def DM_spectral_decomposition(self, initial_dm):
        '''
        If A = VDV-1 --> sqrt(A) = Vsqrt(D)V-1 where D is diagnoal matrix which consists of eigenvalues
        V = eigenvectors
        '''
        eigval = np.linalg.eigh(initial_dm)[0]
        eigval = np.abs(eigval) #due to small imaginary at gamma
        sqrt_eigval = np.sqrt(eigval)

        V = np.linalg.eigh(initial_dm)[1]
        V_inv = np.linalg.inv(V)

        size = len(sqrt_eigval)
        sqrt_D = np.zeros((size, size), dtype=complex)

        for i in range(size):
            sqrt_D[i,i] = sqrt_eigval[i]

        decomp_dm = np.dot(np.dot(V, sqrt_D), V_inv)


        decomp_dm = (decomp_dm + decomp_dm.conj().transpose()) / 2
        return decomp_dm

    def make_phTB_H_ver2(self, dm):
        dm_new = np.zeros((len(dm)*2,len(dm)*2), dtype=complex)
        gamma = np.zeros((len(dm),len(dm)), dtype=complex)
    
        if len(self.TRS_broken) == self.num_atom:
            alpha_set = self.TRS_broken
        elif len(self.TRS_broken) == 1:
            #print 'constant alpha'
            alpha_set = []
            for i in range(self.num_atom):
                alpha_set.append(self.TRS_broken[0]) 
        #print alpha_set
        else:
            print 'alpha error'

        #print alpha_set

        for i in range(self.num_atom):
            alpha = alpha_set[i]
            #print alpha
            if self.dimension == 2:
                gamma[i*self.dimension+0,i*self.dimension+1] += alpha[2]
                gamma[i*self.dimension+1,i*self.dimension+0] += alpha[2] *-1
            if self.dimension == 3:
                gamma[i*self.dimension+0,i*self.dimension+1] += alpha[2] *-1
                gamma[i*self.dimension+0,i*self.dimension+2] += alpha[1] 
                gamma[i*self.dimension+1,i*self.dimension+2] += alpha[0] *-1 
                gamma[i*self.dimension+1,i*self.dimension+0] += alpha[2]
                gamma[i*self.dimension+2,i*self.dimension+0] += alpha[1] *-1
                gamma[i*self.dimension+2,i*self.dimension+1] += alpha[0]     

        for i in range(2):
            for j in range(2):
                for k in range(len(dm)):
                    for l in range(len(dm)):
                        position_x = i*len(dm) + k
                        position_y = j*len(dm) + l
                        if i == 0 and j == 0:
                            dm_new[position_x, position_y] = 0.0
                        elif i == 0 and j ==1:
                            #dm_new[position_x, position_y] = dm.transpose()[k,l]
                            dm_new[position_x,position_y] = 1j * dm[k,l]
                        elif i ==1 and j == 0:
                            #dm_new[position_x, position_y] = dm[k,l]
                            dm_new[position_x,position_y] = -1j * dm[k,l]
                        elif i ==1 and j == 1:
                            dm_new[position_x, position_y] = -2j* gamma[k,l]
                        else:
                            pass
        return dm_new



    def get_phonon_band(self, q_path, q_spacing):
        vasp2THZ = 15.633302
        '''
        First, obtain q_path
        '''
        special_q_distance_list, q_distance_list, q_vec_list = self.obtain_qpath(q_path, q_spacing)


        '''
        Second, solve dynamical matrix at a q point
        '''
        
        band_structure = []
        atom_projected = []
        pos_expectation = []

        pos_expectation_nk = []
        if self.edge_cal:
            position_operator = np.zeros((2*self.dimension*self.num_atom, 2*self.dimension*self.num_atom), dtype=complex)
            for i in range(self.num_atom):
                cell_number = i % self.num_repeat
                for j in range(self.dimension):
                    position_operator[self.dimension*i + j, self.dimension*i + j] = cell_number
                    position_operator[self.dimension*i + j + self.dimension*self.num_atom, self.dimension*i + j + self.dimension*self.num_atom] = cell_number

        for i in range(len(q_vec_list)):
            print 'Process: ' + str(i+1) +'/' + str(len(q_vec_list))
            q_vec = q_vec_list[i]
            dyn = self.construct_dynamicalmatrix_q(q_vec)
            modified_dyn = self.DM_spectral_decomposition(dyn)
            modified_dyn = self.make_phTB_H_ver2(modified_dyn)
            #modified_dyn = transform_H_for_u(modified_dyn,num_atom, latt_vec, atom_pos, recip_vec, q_vec, dimension)
            #print modified_dyn[6][6]
            #w1 = (np.linalg.eigvalsh(modified_dyn).real) *vasp2THZ
            w1, v1 = np.linalg.eigh(modified_dyn)
            band_num = len(w1)
            band_structure.append(w1)

            """atom projection module"""

            atom_projected_nk = []
            for j in range(band_num):
                wave_square = np.multiply(v1[:,j].transpose(), np.conjugate(v1[:,j].transpose()))
                atom_temp = []
                for k in range(self.num_atom):
                    initial = (band_num/2)  + self.dimension*(k)
                    end = (band_num/2)  + self.dimension*(k+1)
                    #print initial, end
                    projection = np.sum(wave_square[initial:end]) + np.sum(wave_square[initial-(band_num/2):end-(band_num/2)])
                    atom_temp.append(np.real(projection))
                #print np.sum(atom_temp)
                atom_projected_nk.append(atom_temp)

            atom_projected.append(atom_projected_nk)

            """atom projection module"""

            """position expectation module"""
            if self.edge_cal:
                pos_expectation_nk = []
                for j in range(band_num):
                    expectation = np.dot(np.conjugate(v1[:,j]), np.dot(position_operator, v1[:,j].transpose()))
                    #expectation = np.matmul(np.conjugate(v1[j]), v1[j].transpose())
                    pos_expectation_nk.append(expectation)
                pos_expectation.append(pos_expectation_nk)


        print 'total number of band is ' + str(band_num)


        '''
        Third, write phonon frequency file
        '''

        out_name = 'ph_frequecny_'+self.out_tag+'.out'
        g = open(out_name,'w')
        out_name2 = 'ph_frequecny_'+self.out_tag+'_projected.out'
        g2 = open(out_name2,'w')
        if self.edge_cal:
            print "Edge calculations: position_expectation calcualtion"
            out_name3 = 'ph_frequecny_'+self.out_tag+'_pos_expectation.out'
            g3 = open(out_name3, 'w')

        templine = ''

        for i in range(len(special_q_distance_list)):
            templine += str(special_q_distance_list[i])+ ' '
        templine += '\n'
        g.write(templine)
        g2.write(templine)
        if self.edge_cal:
            g3.write(templine)

        for i in range(len(q_vec_list)):
            templine = str(q_distance_list[i])
            for j in range(band_num):
                templine += ' '+ str(band_structure[i][j]*vasp2THZ)
            templine += '\n'    
            g.write(templine)
        g.close()

        for i in range(len(q_vec_list)):
            for j in range(band_num):
                templine = str(q_distance_list[i]) + ' ' + str(band_structure[i][j])
                for k in range(self.num_atom):
                    templine += ' ' + str(atom_projected[i][j][k])
                templine += '\n'
                g2.write(templine)
            g2.write('\n')
        g2.close()

        if self.edge_cal:
            for i in range(len(q_vec_list)):
                for j in range(band_num):
                    templine = str(q_distance_list[i]) + ' ' + str(band_structure[i][j]*vasp2THZ) + ' ' + str(pos_expectation[i][j])
                    templine += '\n'
                    g3.write(templine)
                g3.write('\n')
            g3.close()

        return 0

    def draw_phonon_band(self):
        '''
        Fourth, draw phonon band along q path
        '''
        file_name = 'ph_frequecny_'+self.out_tag+'.out'
        f = open(file_name,'r')
        tempf = f.readlines()
        f.close()

        totalline = len(tempf)
        #print totalline
        band_num = len(tempf[1].split()) - 1


        sqx = [float(tempf[0].split()[i]) for i in range(len(tempf[0].split()))]

        #print sqx
        fig = plt.figure()
        plt.axhline(y=0, color='black', linewidth=2)
        for i in range(len(sqx)):
            plt.axvline(x=sqx[i], color='black', linewidth=2)

        qx = []
        eigenval = np.zeros((band_num, totalline-1))


        for i in range(1, totalline):
            temp_val = tempf[i].split()
            qx.append(float(temp_val[0]))
            for j in range(band_num):
                eigenval[j][i-1] = float(temp_val[1+j])


        for i in range(band_num):
            plt.plot(qx, eigenval[i], color='black')

        #plt.ylim(0, max(eigenval[:][:]))
        plt.xlim(min(sqx)-0.1, max(sqx)+0.1)
        fig.savefig('phband.png')
        plt.show()
        
        
        return 0

    def draw_phonon_band_edge_projection(self):
        if not self.edge_cal:
            print "ERROR: You didn't perform edge calculations"
            return 0

        f = open('ph_frequecny_'+self.out_tag+'_pos_expectation.out' , 'r')
        #f = open('ph_frequecny_'+self.out_tag+'_projected.out' , 'r')
        tempf = f.readlines()
        f.close()

        band_num = 2* self.dimension * self.num_atom
        kpoint_num = (len(tempf)-1) / (band_num+1)
        repeat_info = np.array([float(i) % self.num_repeat for i in range(self.num_atom)])
        print band_num, kpoint_num

        sqx = [float(tempf[0].split()[i]) for i in range(len(tempf[0].split()))]
        fig = plt.figure()
        plt.axhline(y=0, color='black', linewidth=2)
        for i in range(len(sqx)):
            plt.axvline(x=sqx[i], color='black', linewidth=2)

        qx = []
        eigenval = np.zeros((band_num, kpoint_num))
        projected = np.zeros((band_num, kpoint_num))

        for i in range(kpoint_num):
            line = (band_num+1) * i  + 1
            qx.append(float(tempf[line].split()[0]))
            for j in range(band_num):
                temp_val = tempf[line+j].split()
                frequency = float(temp_val[1])
                #atominfo = np.array([float(temp_val[k+2]) for k in range(self.num_atom)])
                #edgeinfo = np.sum(atominfo * repeat_info)
                #print np.real(complex(temp_val[2]))
                edgeinfo = np.real(complex(temp_val[2]))
                

                eigenval[j][i] = frequency
                projected[j][i] = edgeinfo

        bubble_size = 40
        for i in range(band_num/2, band_num):
            plt.plot(qx, eigenval[i], linewidth=0.3, color='black')
            plt.scatter(qx, eigenval[i], bubble_size, c=projected[i], cmap='RdBu', vmin=0, vmax=self.num_repeat-1, edgecolors='face')


        plt.xlim(min(sqx)-0.1, max(sqx)+0.1)
        fig.savefig('phband_edge.png')
        plt.show()

        return 0

    def draw_phonon_projected_band(self, PartA, PartB):
        '''
        draw phonon band projected on two parts along q-path
        '''
        return 0


    def make_phband_PROCAR_format(self, q_grid):
        vasp2THZ = 15.633302
        #recip_vec = self.reciprocal()
        '''
        First, generate q points depending on a selected mode
        Input examples
        #q_grid = ['slice',[31, 31, 1], 0.5]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
        #q_grid = ['line', q_path7, 20, 2, 51 ]  #### [[fixed_qpoint_direction, fixed_qpoint_value], [mxn for 2D planes]]
        #q_grid = ['node', [0.5,0.5,0.5], 0.005, 21, 21]  #### [q_node mode, center_direct_coord, radius from the node. theta_spacing, phi_spacing]
        #q_grid = ['berrycurv', [0, 1.8262], [0, 1.111], [2, 0.5], 101, 101] ###
        #q_grid = ['berryphase', q_path_berry_MoSe2, 10]  ####
        '''
        q_vec_list = []
        #print type(q_grid[0])

        if q_grid[0] == 'slice':
            mode_name = 'slicemode'
            print 'Dimension of your system is ' + str(self.dimension)
            fixed_direction = np.argwhere(np.array(q_grid[1])==1)[0][0]
            for k in range(q_grid[1][2]):
                for j in range(q_grid[1][1]):
                    for i in range(q_grid[1][0]):
                        #kx, ky, kz range --> -0.5 < kx, ky, kz < 0.5  !!! -0.5 = 0.5 for periodicity
                        deltax = 1.0 /(q_grid[1][0])
                        deltay = 1.0 /(q_grid[1][1])
                        deltaz = 1.0 /(q_grid[1][2])
                        temp = [-0.5 + deltax*i, -0.5 + deltay*j, -0.5 + deltaz*k] # For 2D case
                        #temp = [-0.5 + deltay*j, -0.5 + deltax*i, -0.5 + deltaz*k] # For 2D case
                        temp[fixed_direction] = float(q_grid[2])
                        if self.dimension == 2:
                            temp = temp[:2]
                        #temp = temp - np.round(temp)
                        #print temp
                        q_vec_list.append(temp)
        elif q_grid[0] == 'line':
            if self.dimension !=3:
                print 'Dimension ERROR:'
                return 0
            mode_name = 'linemode'
            fixed_direction = int(q_grid[3])
            fixed_spacing = int(q_grid[4])
            q_path = np.array(q_grid[1])
            q_spacing = int(q_grid[2])

            q_temp_list = []
            for i in range(len(q_path)-1):
                initial = np.array(q_path[i])
                final = np.array(q_path[i+1])
                for j in range(q_spacing):
                    delta = (final - initial) / float(q_spacing-1)
                    temp = initial + delta*j
                    q_temp_list.append(temp)
        
            for i in range(len(q_temp_list)):
                for j in range(fixed_spacing):
                    tempp = q_temp_list[i]
                    delta = 1.0 /(fixed_spacing)
                    tempp[fixed_direction] = -0.5 + delta*j
                    #print temp
                    q_vec_list.append(tempp)
                    #fixed part doesn't change!!!!!!!!!! I don't know!!!!!!!!!!!!!!
            print 'Not implemented yet'
            return 0
        elif q_grid[0] == 'node':
            if self.dimension !=3:
                print 'Dimension ERROR:, currently support only 3D cases'
                return 0
            mode_name = 'nodemode'
            center_direct = np.array(q_grid[1]) ; r = float(q_grid[2]) ; theta_spacing = int(q_grid[3]) ; phi_spacing = int(q_grid[4])
            #print center_direct, r, theta_spacing, phi_spacing
            center_cart = np.dot(self.recip_vec.transpose(), np.array(center_direct).transpose())
            inv_recip_vec = np.linalg.inv(self.recip_vec.transpose())
            for i in range(theta_spacing):
                for j in range(phi_spacing):
                    # phi --> periodic from 0 to 2pi,  theta --> not periodic from 0 to pi
                    delta_theta = np.pi / (theta_spacing-1)
                    delta_phi = 2*np.pi / phi_spacing
                    #temp_cart = np.array([r*np.sin(delta_theta*i)*np.cos(delta_phi*j), r*np.sin(delta_theta*i)*np.sin(delta_phi*j), r*np.cos(delta_theta*i)])
                    temp_cart = np.array([r*np.sin(delta_theta*i)*np.cos(-np.pi + delta_phi*j), r*np.sin(delta_theta*i)*np.sin(-np.pi + delta_phi*j), r*np.cos(delta_theta*i)])
                    temp_cart += center_cart
                    #temp_cart = temp_cart  /  (2*np.pi)
                    temp_direct = np.dot(inv_recip_vec, temp_cart.transpose()).transpose()
                    #temp_direct[1] = temp_direct[1] - np.round(temp_direct[1])
                    q_vec_list.append(temp_direct)
            #print len(q_vec_list)
        elif q_grid[0] == 'berrycurv_slice':
            mode_name = 'berrycurv_slice'
            qx_range = np.array(q_grid[1])
            qy_range = np.array(q_grid[2])
            qz_infor = np.array(q_grid[3])
            Nqx = int(q_grid[4])
            Nqy = int(q_grid[5])
            deltax = (qx_range[1]-qx_range[0]) / float(Nqx-1)
            deltay = (qy_range[1]-qy_range[0]) / float(Nqy-1)
            inv_recip_vec = np.linalg.inv(self.recip_vec.transpose())
            for i in range(Nqx+1):
                for j in range(Nqy+1):
                    temp_cart = np.array([qx_range[0]+deltax*i, qy_range[0]+deltay*j, 0])
                    temp_direct = np.dot(inv_recip_vec, temp_cart.transpose()).transpose()
                    temp_direct[int(qz_infor[0])] = float(qz_infor[1])
                    if self.dimension == 2:
                        temp_direct = temp_direct[:2]
                    q_vec_list.append(temp_direct)
            self.prepare_3Dplot_data([q_grid[4]+1, q_grid[1][0], (q_grid[1][1]-q_grid[1][0])/(q_grid[4]-1)], [q_grid[5]+1, q_grid[2][0],  (q_grid[2][1]-q_grid[2][0])/(q_grid[5]-1)])
        elif q_grid[0] == 'berryphase':
            mode_name = 'berryphase'
            q_path = np.array(q_grid[1])
            q_spacing = int(q_grid[2])

            q_vec_list = []
            for i in range(len(q_path)-1):
                initial = np.array(q_path[i])
                final = np.array(q_path[i+1])
                for j in range(q_spacing):
                    delta = (final - initial) / float(q_spacing)
                    temp = initial + delta*j
                    q_vec_list.append(temp)
        else:
            print 'There is no mode I can calculate'
            print 'Check your mode'
            print 'Mode name should belongs to the slice, line, node, berrycurv, berryphase'
            return 0

        '''
        Second, set the file name
        '''

        bandfilename = 'phband_PROCAR_'+self.out_tag+'_'+mode_name+'.out'
        bf = open(bandfilename, 'w')
        if q_grid[0] == 'slice':
            bf.write(str(q_grid[1][0])+'\t'+str(q_grid[1][1])+'\t'+str(q_grid[1][2])+'\t'+str(self.num_atom*self.dimension*2)+'\n') # For q_slice mode
        if q_grid[0] == 'line':
            bf.write(str(q_grid[4])+'\t'+str(len(q_grid[1])*q_grid[2])+'\t'+'1'+'\t'+str(self.num_atom*self.dimension*2)+'\n') # For q_line mode
        if q_grid[0] == 'node':
            bf.write(str(q_grid[4])+'\t'+str(q_grid[3])+'\t'+'1'+'\t'+str(self.num_atom*self.dimension*2)+'\n') # For q_node mode
        if q_grid[0] == 'berrycurv_slice':
            bf.write(str(q_grid[4])+'\t'+str(q_grid[5])+'\t'+'1'+'\t'+str(self.num_atom*self.dimension*2)+'\n') # For berry curvature mode
        if q_grid[0] == 'berryphase':
            bf.write(str(len(q_vec_list))+'\t'+str(self.num_atom*self.dimension*2)+'\n') # For berry curvature mode
        bf.write('\n')

        '''
        Third, solve dynamical matrix
        '''

        band_structure = []

        for i in range(len(q_vec_list)):
            print 'Process: ' + str(i+1) +'/' + str(len(q_vec_list))
            q_vec = q_vec_list[i]
            G = np.round(q_vec) 
            print G
            q_vec_refined = q_vec - G 
            dyn = self.construct_dynamicalmatrix_q(q_vec_refined)
            modified_dyn = self.DM_spectral_decomposition(dyn)
            modified_dyn = self.make_phTB_H_ver2(modified_dyn)
            w1, v1 = np.linalg.eigh(modified_dyn)
            q = np.dot(self.recip_vec.transpose(), np.array(q_vec).transpose())
            ##################################################################
            refined_v1 = np.zeros((len(w1),len(w1)), dtype=complex)
            G_cart = np.dot(self.recip_vec.transpose(), np.array(G).transpose())
            x = np.dot(self.latt_vec.transpose(), np.array(self.atom_pos).transpose())
            product = np.dot(G_cart,x)
            for cc in range(len(w1)):
                for bb in range(self.num_atom):
                    for aa in range(self.dimension):
                        refined_v1[bb*self.dimension + aa][cc] = v1[bb*self.dimension + aa][cc] * np.exp(-1.j*product[bb])
                        refined_v1[bb*self.dimension + aa + len(w1)/2][cc] = v1[bb*self.dimension + aa + len(w1)/2][cc] * np.exp(-1.j*product[bb])
            ##################################################################
            band_num = len(w1)
            if self.dimension ==3:
                line = 'kx' + '\t' + 'ky' + '\t' + 'kz' + '\t' + str(q_vec[0]) + '\t' + str(q_vec[1]) + '\t' + str(q_vec[2]) + '\t' + str(q[0]) + '\t' + str(q[1]) + '\t' + str(q[2]) + '\n'
            if self.dimension ==2:
                line = 'kx' + '\t' + 'ky' + '\t' + 'kz' + '\t' + str(q_vec[0]) + '\t' + str(q_vec[1]) + '\t' + '0' + '\t' + str(q[0]) + '\t' + str(q[1]) + '\t' + '0' + '\n'
            bf.write(line)
            for j in range(band_num):
                line = str(j+1) + '\t' + str(w1[j]) + '\t'
                for k in range(band_num):
                    #line += str(v1[k][j]) + '\t'
                    line += str(refined_v1[k][j]) + '\t'
                line += '\n'
                bf.write(line)
            bf.write('\n')
        bf.close()

        return 0    

    def get_3Dplot_data(self, band):
        vasp2THZ = 15.633302
        filename ='phband_PROCAR_'+self.out_tag+'_berrycurv'+'.out'
        f = open(filename, 'r')
        tempf = f.readlines()
        f.close()

        nx, ny, nz, nb = int(tempf[0].split()[0])+1, int(tempf[0].split()[1])+1, int(tempf[0].split()[2]), int(tempf[0].split()[3])

        eigenval = []

        for i in range(ny):
            for j in range(nx):
                position = (nb+2)*nx*i + (nb+2)*j +(1+band) + 2
                #print position
                eigenval.append(float(tempf[position].split()[1]))

        eigenval = np.array(eigenval) * vasp2THZ
        eigenval_rsh = eigenval.reshape((nx, ny))

        #print eigenval_rsh[1][1]
    
        outname = '3Dplot_'+self.out_tag+'_band_'+str(band)+'.out'
        g = open(outname, 'w')

        for i in range(ny):
            temp_line = ''
            for j in range(nx):
                temp_line += str(eigenval_rsh[i][j]) + '\t'
            temp_line += '\n'
            g.write(temp_line)

        g.close()



        return 0

    def prepare_3Dplot_data(self, x_info, y_info):
        fx = open('3Dplot_xgrid', 'w')
        fy = open('3Dplot_ygrid', 'w')
    
        num_x, start_x, diff_x = x_info[0], x_info[1], x_info[2]
        num_y, start_y, diff_y = y_info[0], y_info[1], y_info[2]
    
        for i in range(num_y):
            linex = ''
            liney = ''
            for j in range(num_x):
                linex += str(start_x + j*diff_x) + '\t'
                liney += str(start_y + i*diff_y) + '\t'
            linex += '\n'
            liney += '\n'
            fx.write(linex)
            fy.write(liney)

        fx.close()
        fy.close()
    
        return 0

    def _set_cell_oriented(self):
        self._angles = get_angles(self.latt_vec)
        self._cell_params = get_cell_parameters(self.latt_vec)
        a, b, c = self._cell_params
        alpha, beta, gamma = self._angles
        #print alpha, beta, gamma, a, b, c
        self._lattice_oriented = get_cell_matrix(a, b, c, alpha, beta, gamma)
        self._positions_oriented = self._get_oriented_displacements(np.dot(self.atom_pos, self.latt_vec))
        return 0
    
    def _get_oriented_displacements(self, vec_cartesian):
        return np.dot(np.dot(vec_cartesian, np.linalg.inv(self.latt_vec)), self._lattice_oriented)
    

    def make_anime_file_for_vsim(self, q_point):
        if self.dimension == 2:
            print 'Error: current version only supports 3-dim case'
            return 0
        self._set_cell_oriented()
        lat_oriented = self._lattice_oriented
        pos_oriented = self._positions_oriented
        vasp2THZ = 15.633302

        filename ='anime_'+self.out_tag+'.ascii'
        f = open(filename, 'w')
        f.write('# PhononTB generated file for v_sim 3.6'+'\n')

        templine = '   '+str(lat_oriented[0][0]) + '\t' +str(lat_oriented[1][0]) + '\t' +str(lat_oriented[1][1]) + '\n'
        f.write(templine)
        templine = '   '+str(lat_oriented[2][0]) + '\t' +str(lat_oriented[2][1]) + '\t' +str(lat_oriented[2][2]) + '\n'
        f.write(templine)

        real_atom_pos = np.dot(self.latt_vec.transpose(), np.array(self.atom_pos).transpose()).transpose()
        #print real_atom_pos

        for i in range(self.num_atom):
            templine = '   '+str(pos_oriented[i][0]) + '\t' +str(pos_oriented[i][1]) + '\t' +str(pos_oriented[i][2]) + '\t' + 'A' + '\n'
            f.write(templine)    
        

        q_vec = q_point
        dyn = self.construct_dynamicalmatrix_q(q_vec)
        modified_dyn = self.DM_spectral_decomposition(dyn)
        modified_dyn = self.make_phTB_H_ver2(modified_dyn)
        w1, v1 = np.linalg.eigh(modified_dyn)
        #print np.linalg.norm(v1[:,17])
    
        #w2 = np.linalg.eigvalsh(dyn).real
        #w = refine_frequency(w2)*vasp2THZ
        band_num = len(w1)

        for i in range(band_num/2, band_num):
            templine_start = '#metaData:  qpt=[' + str(q_vec[0]) + ';' + str(q_vec[1]) + ';' + str(q_vec[2]) + ';' + str(w1[i]*vasp2THZ) + '\\' + '\n'
            f.write(templine_start)
            for j in range(self.num_atom):
                index = j 
                real_part_x, imag_part_x = np.real(v1[:,i][3*index + 0]), np.imag(v1[:,i][3*index + 0])
                real_part_y, imag_part_y = np.real(v1[:,i][3*index + 1]), np.imag(v1[:,i][3*index + 1])
                real_part_z, imag_part_z = np.real(v1[:,i][3*index + 2]), np.imag(v1[:,i][3*index + 2])
                templine = '#; ' + str(real_part_x) + '; ' + str(real_part_y) + '; ' + str(real_part_z) + '; ' + str(imag_part_x) + '; ' + str(imag_part_y) + '; ' + str(imag_part_z) + '; ' +' \\' + '\n'
                f.write(templine)
            templine_end = '# ]' + '\n'
            f.write(templine_end)

        f.close()
        return 0


class ComputeTopologicalInvariants:
    def __init__(self, out_tag, band_range, q_grid):
        self.out_tag = out_tag
        self.q_grid = q_grid
        if self.q_grid[0] == 'slice':
            self.phband_file = 'phband_PROCAR_'+self.out_tag+'_slicemode.out'
        elif self.q_grid[0] == 'line':
            self.phband_file = 'phband_PROCAR_'+self.out_tag+'_linemode.out'
        elif self.q_grid[0] == 'node':
            self.phband_file = 'phband_PROCAR_'+self.out_tag+'_nodemode.out'
        elif self.q_grid[0] == 'berrycurv_line':
            self.phband_file = 'phband_PROCAR_'+self.out_tag+'_berrycurve_line.out'
        elif self.q_grid[0] == 'berrycurv_slice':
            self.phband_file = 'phband_PROCAR_'+self.out_tag+'_berrycurve_slice.out'
        elif self.q_grid[0] == 'berryphase':
            self.phband_file = 'phband_PROCAR_'+self.out_tag+'_berryphase.out'
        self.band_range = band_range
        self.q_grid = q_grid
        self.br = len(self.band_range)
        print str(self.phband_file)


    def get_Willsons_loop(self):
        if self.q_grid[0] == 'slice':
            print 'Slice mode calculation'
            print '======================'
            print 'Information'
            print 'Nx, Ny, Nz = ' +str(self.q_grid[1][0]) + ' ' +str(self.q_grid[1][1]) + ' ' +str(self.q_grid[1][2]) + ' '
            print 'The q-point for fixed values = ' + str(self.q_grid[2])
            print '======================'
            self.calculate_all_theta_slicemode()
            return 0
        if self.q_grid[0] == 'line':
            print 'Line mode calculation (Not implemented yet)'
            #print '======================'
            #print 'Information'
            #print 'Nx, Ny, Nz = ' +str(q_grid[1][0]) + ' ' +str(q_grid[1][1]) + ' ' +str(q_grid[1][2]) + ' '
            #print 'The q-point for fixed values = ' + str(q_grid[2])
            #print '======================'
            # self.calculate_all_theta_linemode()
            return 0
        if self.q_grid[0] == 'node':
            print 'Node mode calculation'
            print '======================'
            print 'Information'
            print 'Center of node = ' +str(self.q_grid[1][0]) + ' ' +str(self.q_grid[1][1]) + ' ' +str(self.q_grid[1][2]) + ' '
            print 'Radius = ' + str(self.q_grid[2])
            print 'N_theta, N_phi = ' +str(self.q_grid[3]) + ' ' +str(self.q_grid[4]) + ' '
            print '======================' 
            # self.calculate_all_theta_nodemode()
            return 0


    def calculate_theta_at_fixed_ky(self, ky):
        # band range: bands for calculting theta below an universal gap
        f = open(self.phband_file,'r')
        tempf = f.readlines()
        f.close()
    
        mx, my, mz, nb = int(tempf[0].split()[0]), int(tempf[0].split()[1]), int(tempf[0].split()[2]), int(tempf[0].split()[3])

        if ky > my:
            print 'Your ky point is out of range'
            print ky, my
            return 0
        else:
            ky_pos = -0.5 + (1.0/my)*(ky-1)
            #print 'ky information: ' + str(ky) +'th = ' + str(ky_pos)


        ''' First, read wavefunction and make proper array '''


        wf_data = [] # wf_data[ith kx point][nth ph band]

        for i in range(mx):
            wf_data.append([])
            start_line = 2 + ((1+nb+1)*mx)*(ky-1) + (1+nb+1)*i
            for j in range(nb):
                tempf_wf = []
                for k in range(nb):
                    tempf_wf.append(complex(tempf[start_line+j+1].split()[2+k]))
                wf_data[i].append(tempf_wf)

        ''' Second, calculate F and D matrix ''' 

        D = np.identity(self.br, dtype=complex)
        H = np.identity(self.br, dtype=complex)
  
        for i in range(mx):
            if i == mx-1:
                A_point = wf_data[i]
                B_point = wf_data[0]
            else:
                A_point = wf_data[i]
                B_point = wf_data[i+1]

            temp_F_matrix = self.get_F_matrix(A_point,B_point)
            #print temp_F_matrix
            temp_G_matrix = self.get_G_matrix(A_point,B_point)
            U, s, V = np.linalg.svd(temp_F_matrix, full_matrices=True)
            temp_F_matrix = np.dot(U, V)
            D = np.dot(D, temp_F_matrix)
            H = np.dot(H, temp_G_matrix)

        ''' Third, diagonalize D matrix and calculate theta at fixed ky'''

        Eigen = np.linalg.eigvals(D)
        A = np.imag(np.log(Eigen)) / (-2.0*np.pi)
        A = np.sort(A)
    
        return ky_pos, A


    def get_F_matrix(self, A_point, B_point):
        # A_point : kx,i   nb bands and wavefunctions
        # B_point : kx,i+1 nb bands and wavefunctions
        F = np.zeros((self.br, self.br), dtype=complex) 
        for i in range(self.br):
            for j in range(self.br):
                x = np.conjugate(A_point[self.band_range[i]])
                y = B_point[self.band_range[j]]
                F[i,j] = np.dot(x,np.transpose(y))
        return F

    def get_G_matrix(self, A_point, B_point):
        # A_point : kx,i   nb bands and wavefunctions
        # B_point : kx,i+1 nb bands and wavefunctions
        G = np.zeros((self.br, self.br), dtype=complex) 
        for i in range(self.br):
            for j in range(self.br):
                x = np.conjugate(A_point[self.band_range[i]])
                y = B_point[self.band_range[j]]
                if i == j:
                    G[i,j] = np.dot(x,np.transpose(y))
                else:
                    G[i,j] = 0.0
        return G



    def calculate_all_theta_slicemode(self):
        x = [] ; y = []

        fig = plt.figure()

        for i in range(self.br):
            y.append([])

        for i in range(self.q_grid[1][1]):
            ky_pos, A = self.calculate_theta_at_fixed_ky(i+1)
            x.append(ky_pos)
            for j in range(self.br):
                y[j].append(A[j])
            #line = str(ky_pos) + '\t'
            #for j in range(len(band_range)):
            #    line += str(A[j]) + '\t'
            #print line

        ky_pos, A = self.calculate_theta_at_fixed_ky(1)
        #print ky_pos
        #print A
        x.append(0.5)
        for j in range(self.br):
            y[j].append(A[j])
        #line = str(abs(ky_pos)) + '\t'
        #for j in range(len(band_range)):
        #    line += str(A[j]) + '\t'
        #print line

        for i in range(self.br):
            plt.plot(x, y[i], linestyle="none", marker="o", ms=5, markeredgecolor='black', markerfacecolor='black')
        plt.title(str(self.band_range))
        plt.xlabel(r'$k$-point')
        plt.ylabel(r'Wannier Charge Center')
        plt.axis([-0.5 -0.001 , 0.5 + 0.001, -0.5, 0.5])
        fig.savefig('WCC.png')

        plt.show()
    
        return 0

    def calculate_Berry_phase(self):
    # band range: bands for calculting theta below an universal gap
        f = open(self.phband_file,'r')
        tempf = f.readlines()
        f.close()

        nkx, nb = int(tempf[0].split()[0]), int(tempf[0].split()[1])

        ''' First, read wavefunction and make proper array '''

        wf_data = [] # wf_data[ith kx point][nth ph band]

        for i in range(nkx):
            wf_data.append([])
            start_line = 2  + (1+nb+1)*i
            #print start_line
            for j in range(nb):
                tempf_wf = []
                for k in range(nb):
                    #print start_line+j+1
                    tempf_wf.append(complex(tempf[start_line+j+1].split()[2+k]))
                wf_data[i].append(tempf_wf)

        ''' Second, calculate F and D matrix ''' 
   
        F =  np.identity(self.br, dtype=complex)
        phi =  np.linalg.det(F)

        for i in range(nkx):
            temp = np.zeros((self.br, self.br), dtype=complex)
            if i == nkx-1:
                for k in range(self.br):
                    for j in range(self.br):
                        A = np.conjugate(wf_data[i][self.band_range[k]])
                        B = np.array(wf_data[0][self.band_range[j]])
                        temp[k,j] = np.dot(A,B)
            else:
                for k in range(self.br):
                    for j in range(self.br):
                        A = np.conjugate(wf_data[i][self.band_range[k]])
                        B = np.array(wf_data[i+1][self.band_range[j]])
                        temp[k,j] = np.dot(A,B)
            phi = phi * np.linalg.det(temp)

        berry = -1 * np.imag(np.log(phi))
        print 'Berry phase is ' + str(berry)
        return 0

###########################################################################
def get_angles(lattice):
    a, b, c = get_cell_parameters(lattice)
    alpha = np.arccos(np.vdot(lattice[1], lattice[2]) / b / c) / np.pi * 180
    beta  = np.arccos(np.vdot(lattice[2], lattice[0]) / c / a) / np.pi * 180
    gamma = np.arccos(np.vdot(lattice[0], lattice[1]) / a / b) / np.pi * 180
    return alpha, beta, gamma

def get_cell_parameters(lattice):
    return np.sqrt(np.dot (lattice, lattice.transpose()).diagonal())

def get_cell_matrix(a, b, c, alpha, beta, gamma):
    # These follow 'matrix_lattice_init' in matrix.c of GDIS
    alpha *= np.pi / 180
    beta *= np.pi / 180
    gamma *= np.pi / 180
    a1 = a
    a2 = 0.0
    a3 = 0.0
    b1 = np.cos(gamma)
    b2 = np.sin(gamma)
    b3 = 0.0
    c1 = np.cos(beta)
    c2 = (2 * np.cos(alpha) + b1**2 + b2**2 - 2 * b1 * c1 - 1) / (2 * b2)
    c3 = np.sqrt(1 - c1**2 - c2**2)
    lattice = np.zeros((3, 3), dtype='double')
    lattice[0, 0] = a
    lattice[1] = np.array([b1, b2, b3]) * b
    lattice[2] = np.array([c1, c2, c3]) * c
    return lattice    

def from_cartesin_to_direct(lattice, cart_coord):
    # This will transform Cartesian coordinates to direct coordinates
    lattice_inv = np.linalg.inv(np.array(lattice).T)
    direct = np.matmul(lattice_inv, np.array(cart_coord).T)
    return direct.T

###########################################################################


if __name__ == "__main__":
    FC_test =  ForceConstant(3, 3)
