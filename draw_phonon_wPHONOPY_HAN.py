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
#from numeric import *

def select_specific_atoms_from_phonopy(filename,start_i, end_i):
    si = int(start_i)
    ei = int(end_i)

    f = open(filename,'r')
    tempf = f.readlines()
    f.close()

    g = open('convert_FC', 'w')

    total = int(tempf[0].split()[0])

    num = end_i - start_i + 1

    g.write(str(num)+'\n')

    for i in range(num):
        initial = si + i
        for j in range(num):
            final = si + j
            ij_line = (initial-1)*total*4 + (final-1)*4 + 1

            templine1 = str(i+1) + ' ' + str(j+1) + '\n'
            g.write(templine1)
            for k in range(1,4):
                g.write(tempf[ij_line+k])
    g.close()


def read_info_cell(primitive_file):
    f = open(primitive_file,'r')
    temp = f.readlines()
    f.close()

    univ = float(temp[1].split()[0])

    lat_a = np.array([float(temp[2].split()[i]) for i in range(3)]) * univ
    lat_b = np.array([float(temp[3].split()[i]) for i in range(3)]) * univ
    lat_c = np.array([float(temp[4].split()[i]) for i in range(3)]) * univ

    latt_vec = [lat_a, lat_b, lat_c]

    num_atom = np.sum(np.array([int(temp[6].split()[i]) for i in range(len(temp[6].split()))]))

    atom_pos = []

    for i in range(num_atom):
        atom_pos.append(np.array([float(temp[8+i].split()[j]) for j in range(3)]))

    #print atom_pos
    #print latt_vec

    return np.array(latt_vec), np.array(atom_pos)

def read_force_constants(force_file):
    f = open(force_file,'r')
    temp = f.readlines()
    f.close()

    data_num = int(temp[0].split()[0])

    fc2 = np.zeros((data_num, data_num,3,3))

    '''
    fc2 format: fc2[i][j][a][b]
      i: Atom index of finitely displaced atom.
      j: Atom index at which force on the atom is measured.
      a, b: Cartesian direction indices = (0, 1, 2) for i and j, respectively
    '''

    for i in range(data_num):
        for j in range(data_num):
            ij_line = (i)*data_num*4 + (j)*4 + 1
            fc2[i][j][0][0] = float(temp[ij_line+1].split()[0]) ; fc2[i][j][0][1] = float(temp[ij_line+1].split()[1]) ; fc2[i][j][0][2] = float(temp[ij_line+1].split()[2])
            fc2[i][j][1][0] = float(temp[ij_line+2].split()[0]) ; fc2[i][j][1][1] = float(temp[ij_line+2].split()[1]) ; fc2[i][j][1][2] = float(temp[ij_line+2].split()[2])
            fc2[i][j][2][0] = float(temp[ij_line+3].split()[0]) ; fc2[i][j][2][1] = float(temp[ij_line+3].split()[1]) ; fc2[i][j][2][2] = float(temp[ij_line+3].split()[2])


    return fc2

def manipulate_forceconstants(fc2, super_size, input_data):
    a_num = len(input_data[1])
    mag_x, mag_y, mag_z = float(input_data[2][0]), float(input_data[2][1]), float(input_data[2][2])
    print mag_x, mag_y, mag_z
    #print fc_size
    for i in range(a_num):
        for j in range(a_num):
            if input_data[1][i] != input_data[1][j]:
                for k in range(super_size):
                    for l in range(super_size):
                        fc2[i*super_size + k][j*super_size + l][0][0] = fc2[i*super_size + k][j*super_size + l][0][0] * mag_x
                        fc2[i*super_size + k][j*super_size + l][1][1] = fc2[i*super_size + k][j*super_size + l][1][1] * mag_y
                        fc2[i*super_size + k][j*super_size + l][2][2] = fc2[i*super_size + k][j*super_size + l][2][2] * mag_z
    return fc2




def get_equivalent_smallest_vectors(atom_number_supercell, atom_number_primitive, super_size, supercell_lattice, supercell_pos, symprec):
    distances = []
    differences = []
    positions = supercell_pos

    # Atomic positions are confined into the lattice made of reduced bases.
    for pos in positions:
        pos -= np.rint(pos)

    p_pos = positions[atom_number_primitive*super_size]
    s_pos = positions[atom_number_supercell]
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                # The vector arrow is from the atom in primitive to
                # the atom in supercell cell plus a supercell lattice
                # point. This is related to determine the phase
                # convension when building dynamical matrix.
                diff = s_pos + [i, j, k] - p_pos
                differences.append(diff)
                vec = np.dot(diff, supercell_lattice)
                distances.append(np.linalg.norm(vec))

    minimum = min(distances)
    smallest_vectors = []
    #print minimum
    for i in range(27):
        if abs(minimum - distances[i]) < symprec:
            smallest_vectors.append(np.dot(differences[i], supercell_lattice))
    #print smallest_vectors
    return smallest_vectors



def get_shortest_vectors(super_latt_vec, super_atom_pos, atom_pos, super_size, symprec):
    """
    shortest_vectors:

      Shortest vectors from an atom in primitive cell to an atom in
      supercell in the fractional coordinates. If an atom in supercell
      is on the border centered at an atom in primitive and there are
      multiple vectors that have the same distance and different
      directions, several shortest vectors are stored. The
      multiplicity is stored in another array, "multiplicity".
      [atom_super, atom_primitive, multiple-vectors, 3]
      
    multiplicity:
      Number of multiple shortest vectors (third index of "shortest_vectors")
      [atom_super, atom_primitive]
    """
    shortest_vectors = np.zeros((super_size*len(atom_pos), len(atom_pos), 27, 3), dtype='double')
    multiplicity = np.zeros((super_size*len(atom_pos), len(atom_pos)), dtype='int')

    for i in range(super_size*len(atom_pos)):
        for j in range(len(atom_pos)):
            vectors = get_equivalent_smallest_vectors(i,j,super_size,super_latt_vec,super_atom_pos,symprec)
            multiplicity[i][j] = len(vectors)
            for k in range(len(vectors)):
                shortest_vectors[i][j][k] = vectors[k]

    return shortest_vectors, multiplicity


def find_atom_number_in_supercell(atom_pos, ith_in_primitive, multiple):
    a_num = len(atom_pos)
    super_size = multiple[0]*multiple[1]*multiple[2]

    if ith_in_primitive > a_num-1:
        print 'your atom is out of range of primitive cell'


    atom_in_supercell = []
    for i in range(super_size):
        atom_in_supercell.append(ith_in_primitive*super_size+i)

    return atom_in_supercell




def construct_dynamical_matrix_q(fc2, super_latt_vec, super_atom_pos, atom_pos, latt_vec, q_vec, multiple, mass, symprec):
    a_num = len(atom_pos)

    super_size = multiple[0]*multiple[1]*multiple[2]

    latt_vec = np.array(latt_vec)
    recip_fac =  (2*np.pi) / np.dot(latt_vec[0], np.cross(latt_vec[1], latt_vec[2]))
    recip_vec = np.array([np.cross(latt_vec[1], latt_vec[2]), np.cross(latt_vec[2], latt_vec[0]), np.cross(latt_vec[0], latt_vec[1])]) * recip_fac

    q = np.dot(recip_vec.transpose(), np.array(q_vec).transpose())
    #print (2*np.pi) / recip_fac
    #print latt_vec
    #print recip_vec / (2*np.pi)
    #print np.dot(recip_vec, q_vec)

    vecs, multiplicity = get_shortest_vectors(super_latt_vec, super_atom_pos, atom_pos, super_size, symprec)

    dm = np.zeros((3 * a_num, 3 * a_num), dtype=complex)

    '''
    D^(a,b)_(k,k') = sum mass(k, k') fc2[k][k']
    '''

    for i in range(a_num):
        for j in range(a_num):
            mass_factor = np.sqrt(mass[i] * mass[j]) 
            dm_local = np.zeros((3, 3), dtype=complex)
            for k in range(super_size*len(atom_pos)):
                if k in find_atom_number_in_supercell(atom_pos,j,multiple):
                    multi = multiplicity[k][i]
                    #multi=1
                    #print multi
                    phase = []
                    for l in range(multi):
                        vec = vecs[k][i][l]
                        phase.append(np.vdot(vec, q) * 1j)
                    phase_factor = np.exp(phase).sum()
                    dm_local += fc2[i*super_size, k] * phase_factor / mass_factor / multi
            dm[(i*3):(i*3+3), (j*3):(j*3+3)] += dm_local


    dynamical_matrix = (dm + dm.conj().transpose()) / 2



    return dynamical_matrix


def DM_spectral_decomposition(initial_dm):
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

def reduce_dimension_to_2D(dm, num_atom):
    dm_reduced = np.zeros((num_atom*2,num_atom*2), dtype=complex)

    for x in range(num_atom*3):
        for y in range(num_atom*3):
            if x == 2 or y == 2:
                dm[x,y] = 0
            if x == 5 or y == 5:
                dm[x,y] = 0
            if x == 8 or y == 8:
                dm[x,y] = 0

    for x in range(num_atom):
        for y in range(num_atom):
            dm_reduced[x*2+0,y*2+0] = dm[x*3+0,y*3+0]
            dm_reduced[x*2+1,y*2+0] = dm[x*3+1,y*3+0]
            dm_reduced[x*2+0,y*2+1] = dm[x*3+0,y*3+1]
            dm_reduced[x*2+1,y*2+1] = dm[x*3+1,y*3+1]

    return dm_reduced

def add_TRS_broken_term(dm, num_atom, alpha):
    dm_reduced = np.zeros((num_atom*2,num_atom*2), dtype=complex)
    for i in range(num_atom):
        dm[i*2+0,i*2+1] += alpha * 1j
        dm[i*2+1,i*2+0] += alpha * -1j


    return dm

def make_special_H(dm,num_atom,alpha0, dimension):
    dm_new = np.zeros((len(dm)*2,len(dm)*2), dtype=complex)
    gamma = np.zeros((len(dm),len(dm)), dtype=complex)
    #dimen = 2
    for i in range(num_atom):
        if i == 0 or i == 3 or i == 4 or i == 7:
            alpha =  +1 * np.array(alpha0)
        if i == 1 or i == 2 or i == 5 or i == 6:
            alpha =  -1 * np.array(alpha0)
        if dimension == 2:
            gamma[i*dimension+0,i*dimension+1] += alpha[2]
            
            gamma[i*dimension+1,i*dimension+0] += alpha[2] *-1
        if dimension ==3:
            gamma[i*dimension+0,i*dimension+1] += alpha[2] *-1
            gamma[i*dimension+0,i*dimension+2] += alpha[1] 
            gamma[i*dimension+1,i*dimension+2] += alpha[0] *-1 
            gamma[i*dimension+1,i*dimension+0] += alpha[2]
            gamma[i*dimension+2,i*dimension+0] += alpha[1] *-1
            gamma[i*dimension+2,i*dimension+1] += alpha[0]     

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


def make_special_H_ver2(dm,num_atom,alpha0, dimension):
    dm_new = np.zeros((len(dm)*2,len(dm)*2), dtype=complex)
    gamma = np.zeros((len(dm),len(dm)), dtype=complex)
    #dimen = 2

    if len(alpha0) == num_atom:
        alpha_set = alpha0
    elif len(alpha0) == 1:
        #print 'constant alpha'
        alpha_set = []
        for i in range(num_atom):
            alpha_set.append(alpha0[0]) 
        #print alpha_set
    else:
        print 'alpha error'

    #print alpha_set

    for i in range(num_atom):
        alpha = alpha_set[i]
        #print alpha
        if dimension == 2:
            gamma[i*dimension+0,i*dimension+1] += alpha[2]
            gamma[i*dimension+1,i*dimension+0] += alpha[2] *-1
        if dimension ==3:
            gamma[i*dimension+0,i*dimension+1] += alpha[2] *-1
            gamma[i*dimension+0,i*dimension+2] += alpha[1] 
            gamma[i*dimension+1,i*dimension+2] += alpha[0] *-1 
            gamma[i*dimension+1,i*dimension+0] += alpha[2]
            gamma[i*dimension+2,i*dimension+0] += alpha[1] *-1
            gamma[i*dimension+2,i*dimension+1] += alpha[0]     

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

def transform_H_for_u(dm,num_atom,latt_vec,atom_pos,recip_vec,q_vec, dimension):
    '''
    This case is for 2D
    '''
    #dimension = 2
    dm_new = np.zeros((len(dm),len(dm)), dtype=complex)

    #q = np.dot(recip_vec.transpose(), np.array(q_vec).transpose())
    q = np.array(q_vec)


    for i in range(2):
        for j in range(2):
            for k in range(num_atom):
                for l in range(num_atom):
                    for ii in range(dimension):
                        for jj in range(dimension):
                            position_i = i*(len(dm)/2) + k*(dimension) + ii
                            position_j = j*(len(dm)/2) + l*(dimension) + jj
                            #x_i = np.dot(latt_vec.transpose(), np.array(atom_pos[k]).transpose())
                            #x_j = np.dot(latt_vec.transpose(), np.array(atom_pos[l]).transpose())
                            x_i = np.array(atom_pos[k])
                            x_j = np.array(atom_pos[l])
                            dm_new[position_i,position_j] = np.exp(-1.0j*q[ii]*x_i[ii]) * dm[position_i,position_j] * np.exp(1.0j*q[jj]*x_j[jj])
    return dm_new


def refine_frequency(matrix):
    find_negative = np.where(matrix<0.0)[0]
    #print find_negative
    sqrt_matrix = np.sqrt(np.abs(matrix))
    final_matrx = sqrt_matrix
    #print final_matrx
    for i in range(len(find_negative)):
        temp = find_negative[i]
        final_matrx[temp] = -final_matrx[temp]
    return final_matrx


def get_phonon_band(out_tag, q_path, primitive_file, supercell_file, FORCE_CONSTANTS_file, multiple, mass, q_spacing, symprec, broken_strength, dimension, manipulate_fc_info):
    super_latt_vec, super_atom_pos = read_info_cell(supercell_file)
    latt_vec, atom_pos = read_info_cell(primitive_file)
    fc2 = read_force_constants(FORCE_CONSTANTS_file)

    super_size = multiple[0] * multiple[1] * multiple[2]
    if manipulate_fc_info[0] == True:  
        fc2 = manipulate_forceconstants(fc2, super_size, manipulate_fc_info)

    num_atom = len(atom_pos)
    vasp2THZ = 15.633302

    latt_vec = np.array(latt_vec)
    recip_fac =  (2*np.pi) / np.dot(latt_vec[0], np.cross(latt_vec[1], latt_vec[2]))
    recip_vec = np.array([np.cross(latt_vec[1], latt_vec[2]), np.cross(latt_vec[2], latt_vec[0]), np.cross(latt_vec[0], latt_vec[1])]) * recip_fac

    #q_spacing = 10 # default

    '''
    First, find q_path and q_distance
    '''

    q_vec_list = []
    q_distance_list = []
    special_q_distance_list = [0.0]
    sq_distance = 0.0

    for i in range(len(q_path)-1):
        initial = np.array(q_path[i])
        final = np.array(q_path[i+1])
        #print np.dot(recip_vec.transpose(), (final-initial).transpose())
        temp_sq_distance = np.linalg.norm(np.dot(recip_vec.transpose(), (final-initial).transpose()))
        sq_distance += temp_sq_distance
        special_q_distance_list.append(sq_distance)
        q_distance = special_q_distance_list[i]
        #print 'Process: ' + str(i+1) +'/' + str(len(q_path)-1)
        for j in range(q_spacing):
            delta = (final - initial) / float(q_spacing-1)
            temp = initial + delta*j
            q_vec_list.append(temp)
            temp_q_distance = np.linalg.norm(np.dot(recip_vec.transpose(), np.array(delta*j).transpose()))
            #print temp_q_distance
            q_distance_list.append(q_distance+temp_q_distance)
    #print len(q_vec_list)
    #print len(special_q_distance_list)
    #print len(q_distance_list)

    '''
    Second, solve dynamcial matrix (only get eigenvalues)
    '''

    band_structure = []
    atom_projected = []
    #print q_vec_list

    for i in range(len(q_vec_list)):
        print 'Process: ' + str(i+1) +'/' + str(len(q_vec_list))
        q_vec = q_vec_list[i]
        dyn = construct_dynamical_matrix_q(fc2, super_latt_vec, super_atom_pos, atom_pos, latt_vec, q_vec, multiple, mass, symprec)
        modified_dyn = DM_spectral_decomposition(dyn)
        if dimension == 2:
            modified_dyn = reduce_dimension_to_2D(modified_dyn, num_atom)
        #modified_dyn = add_TRS_broken_term(modified_dyn, num_atom, broken_strength/2.0)
        modified_dyn = make_special_H_ver2(modified_dyn,num_atom,broken_strength, dimension)
        #modified_dyn = transform_H_for_u(modified_dyn,num_atom, latt_vec, atom_pos, recip_vec, q_vec, dimension)
        #print modified_dyn[6][6]
        #w1 = (np.linalg.eigvalsh(modified_dyn).real) *vasp2THZ
        w1, v1 = np.linalg.eigh(modified_dyn)
    
        #w2 = np.linalg.eigvalsh(dyn).real
        #w = refine_frequency(w2)*vasp2THZ
        band_num = len(w1)
        band_structure.append(w1)

        """atom projection module"""

        atom_projected_nk = []
        for j in range(band_num):
            wave_square = np.multiply(v1.transpose()[j], np.conjugate(v1.transpose()[j]))
            atom_temp = []
            for k in range(num_atom):
                initial = (band_num/2)  + dimension*(k)
                end = (band_num/2)  + dimension*(k+1)
                #print initial, end
                projection = np.sum(wave_square[initial:end]) + np.sum(wave_square[initial-(band_num/2):end-(band_num/2)])
                atom_temp.append(np.real(projection))
            #print np.sum(atom_temp)
            atom_projected_nk.append(atom_temp)

        atom_projected.append(atom_projected_nk)
         


    print 'total number of band is ' + str(band_num)


    '''
    Third, write phonon frequency file
    '''

    out_name = 'ph_frequecny_'+out_tag+'.out'
    g = open(out_name,'w')
    out_name2 = 'ph_frequecny_'+out_tag+'_projected.out'
    g2 = open(out_name2,'w')

    templine = ''

    for i in range(len(special_q_distance_list)):
        templine += str(special_q_distance_list[i])+ ' '
    templine += '\n'
    g.write(templine)
    g2.write(templine)

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
            for k in range(num_atom):
                templine += ' ' + str(atom_projected[i][j][k])
            templine += '\n'
            g2.write(templine)
        g2.write('\n')
    g2.close()


    return 0




def draw_phonon_band(file_name):
    '''
    Fourth, draw phonon band along q path
    '''

    f = open(file_name,'r')
    tempf = f.readlines()
    f.close()

    totalline = len(tempf)
    #print totalline
    band_num = len(tempf[1].split()) - 1


    sqx = [float(tempf[0].split()[i]) for i in range(len(tempf[0].split()))]

    #print sqx
    
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
        plt.plot(qx, eigenval[i])

    plt.xlim(min(sqx)-0.1, max(sqx)+0.1)
    plt.show()
    return 0

def draw_phonon_projected_band(file_name, PartA, PartB):
    f = open(file_name, 'r')
    tempf = f.readlines()
    f.close()

    num_atom = len(PartA) + len(PartB)
    dimension = 3
    band_num = 2* dimension * num_atom
    kpoint_num = (len(tempf)-1) / (band_num+1)
    print band_num, kpoint_num

    sqx = [float(tempf[0].split()[i]) for i in range(len(tempf[0].split()))]

    plt.axhline(y=0, color='black', linewidth=2)
    for i in range(len(sqx)):
        plt.axvline(x=sqx[i], color='black', linewidth=2)

    qx = []
    eigenval = np.zeros((band_num, kpoint_num))
    projected_A = np.zeros((band_num, kpoint_num))
    projected_B = np.zeros((band_num, kpoint_num))

    for i in range(kpoint_num):
        line = (band_num+1) * i  + 1
        qx.append(float(tempf[line].split()[0]))
        for j in range(band_num):
            temp_val = tempf[line+j].split()
            frequency = float(temp_val[1])
            tempa = np.array([float(temp_val[k+2]) for k in PartA])
            tempb = np.array([float(temp_val[k+2]) for k in PartB])
            parta_temp = np.sum(tempa) / (np.sum(tempa) + np.sum(tempb))
            partb_temp = np.sum(tempb) / (np.sum(tempa) + np.sum(tempb))
            
            eigenval[j][i] = frequency
            projected_A[j][i] = parta_temp
            projected_B[j][i] = partb_temp

    #print projected_A - projected_B

    bubble_size = 10
    for i in range(band_num/2, band_num):
        plt.plot(qx, eigenval[i], linewidth=0.3, color='black')
        #plt.scatter(qx, eigenval[i], bubble_size, c=projected_A[i], cmap='Greens')
        plt.scatter(qx, eigenval[i], bubble_size, c=2*np.absolute(projected_A[i]-0.5), cmap='seismic')
        #print projected_A[i]
        #olor_norm = plt.Normalize(projected_A[i].min(), projected_A[i].max())
        #points = np.array([qx, eigenval[i]]).T.reshape(-1, 1, 2)
        #segments = np.concatenate([points[:-1], points[1:]], axis=1)
        #lc = LineCollection(segments, cmap='viridis', norm=norm)
    

    plt.xlim(min(sqx)-0.1, max(sqx)+0.1)
    plt.show()


def make_phband_PROCAR_format(out_tag, q_grid, primitive_file, supercell_file, FORCE_CONSTANTS_file, multiple, mass, symprec, dimension, broken_strength):
    super_latt_vec, super_atom_pos = read_info_cell(supercell_file)
    latt_vec, atom_pos = read_info_cell(primitive_file)
    fc2 = read_force_constants(FORCE_CONSTANTS_file)
    super_size = multiple[0] * multiple[1] * multiple[2]

    vasp2THZ = 15.633302

    latt_vec = np.array(latt_vec)
    recip_fac =  (2*np.pi) / np.dot(latt_vec[0], np.cross(latt_vec[1], latt_vec[2]))
    recip_vec = np.array([np.cross(latt_vec[1], latt_vec[2]), np.cross(latt_vec[2], latt_vec[0]), np.cross(latt_vec[0], latt_vec[1])]) * recip_fac



    '''
    First, make q-points grid for calculating the chern number ex) q_grid = [10,10,1]
    qx first, qy second, qz third
    '''

    q_vec_list = []

    if q_grid[0] == 'slice':
        mode_name = 'slicemode'
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
                    #temp = temp - np.round(temp)
                    #print temp
                    q_vec_list.append(temp)

    
    if q_grid[0] == 'line':
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
        #print q_temp_list

        for i in range(len(q_temp_list)):
            for j in range(fixed_spacing):
                tempp = q_temp_list[i]
                delta = 1.0 /(fixed_spacing)
                tempp[fixed_direction] = -0.5 + delta*j
                #print temp
                q_vec_list.append(tempp)
                #print tempp
        #print q_vec_list # 
        #fixed part doesn't change!!!!!!!!!! I don't know!!!!!!!!!!!!!!
        return 0



    if q_grid[0] == 'node':
        mode_name = 'nodemode'
        center_direct = np.array(q_grid[1]) ; r = float(q_grid[2]) ; theta_spacing = int(q_grid[3]) ; phi_spacing = int(q_grid[4])
        #print center_direct, r, theta_spacing, phi_spacing
        center_cart = np.dot(recip_vec.transpose(), np.array(center_direct).transpose())
        inv_recip_vec = np.linalg.inv(recip_vec.transpose())
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
        print len(q_vec_list)


    if q_grid[0] == 'berry':
        mode_name = 'berrycurv'
        qx_range = np.array(q_grid[1])
        qy_range = np.array(q_grid[2])
        qz_infor = np.array(q_grid[3])
        Nqx = int(q_grid[4])
        Nqy = int(q_grid[5])
        deltax = (qx_range[1]-qx_range[0]) / float(Nqx-1)
        deltay = (qy_range[1]-qy_range[0]) / float(Nqy-1)
        inv_recip_vec = np.linalg.inv(recip_vec.transpose())
        for i in range(Nqx+1):
            for j in range(Nqy+1):
                temp_cart = np.array([qx_range[0]+deltax*i, qy_range[0]+deltay*j, 0])
                temp_direct = np.dot(inv_recip_vec, temp_cart.transpose()).transpose()
                temp_direct[int(qz_infor[0])] = float(qz_infor[1])
                q_vec_list.append(temp_direct)
        print len(q_vec_list)


    if q_grid[0] == 'berryphase':
        mode_name = 'berryphase'
        q_path = np.array(q_grid[1])
        q_spacing = int(q_grid[2])

        q_vec_list = []
        for i in range(len(q_path)-1):
            initial = np.array(q_path[i])
            final = np.array(q_path[i+1])
            for j in range(q_spacing):
                delta = (final - initial) / float(q_spacing-1)
                temp = initial + delta*j
                q_vec_list.append(temp)

        #print q_vec_list
        #return 0



    #q_vec_list = [[0.262,0.4512,0.0], [-0.262,-0.4512,0.0]]
    

    #print q_vec_list
    num_atom = len(atom_pos)
    bandfilename = 'phband_PROCAR_'+out_tag+'_'+mode_name+'.out'
    bf = open(bandfilename, 'w')
    if q_grid[0] == 'slice':
        bf.write(str(q_grid[1][0])+'\t'+str(q_grid[1][1])+'\t'+str(q_grid[1][2])+'\t'+str(num_atom*dimension*2)+'\n') # For q_slice mode
    if q_grid[0] == 'line':
        bf.write(str(q_grid[4])+'\t'+str(len(q_grid[1])*q_grid[2])+'\t'+'1'+'\t'+str(num_atom*dimension*2)+'\n') # For q_line mode
    if q_grid[0] == 'node':
        bf.write(str(q_grid[4])+'\t'+str(q_grid[3])+'\t'+'1'+'\t'+str(num_atom*dimension*2)+'\n') # For q_node mode
    if q_grid[0] == 'berry':
        bf.write(str(q_grid[4])+'\t'+str(q_grid[5])+'\t'+'1'+'\t'+str(num_atom*dimension*2)+'\n') # For berry curvature mode
    if q_grid[0] == 'berryphase':
        bf.write(str(len(q_vec_list))+'\t'+str(num_atom*dimension*2)+'\n') # For berry curvature mode
    bf.write('\n')


    '''
    Second, solve dynamcial matrix for getting eigenvectors
    '''

    #q_vec_list = [[0.6,0.6,0.0], [-0.4, -0.4, 0.0]]

    band_structure = []

    for i in range(len(q_vec_list)):
        print 'Process: ' + str(i+1) +'/' + str(len(q_vec_list))
        q_vec = q_vec_list[i]
        dyn = construct_dynamical_matrix_q(fc2, super_latt_vec, super_atom_pos, atom_pos, latt_vec, q_vec, multiple, mass, symprec)
        modified_dyn = DM_spectral_decomposition(dyn)
        if dimension ==2:
            modified_dyn = reduce_dimension_to_2D(modified_dyn, num_atom)
        #modified_dyn = add_TRS_broken_term(modified_dyn, num_atom, broken_strength/2.0)
        modified_dyn = make_special_H_ver2(modified_dyn,num_atom,broken_strength,dimension)
        #modified_dyn = transform_H_for_u(modified_dyn,num_atom, latt_vec, atom_pos, recip_vec, q_vec, dimension)
        w1, v1 = np.linalg.eigh(modified_dyn)
        #for k in range(dimension*num_atom*2):
        ########################################################################################
        # from psi_nk to u_nk
        #A = np.dot(recip_vec.transpose(), np.array([1.0,0.0,0.0]).transpose())
        #B = np.dot(latt_vec.transpose(), np.array([0.0,1.0,0.0]).transpose())
        #print A[0]*B[0]
        q = np.dot(recip_vec.transpose(), np.array(q_vec).transpose())
        G = np.zeros(3)
        for j in range(3):
            if q_vec_list[i][j] >= 0.5:
                G[j] = 1.0
        #print G
        #G_cart = np.dot(recip_vec.transpose(), np.array(G).transpose())
        #x = np.dot(latt_vec.transpose(), np.array(atom_pos).transpose()).transpose()
        #G_cart = G
        #x = atom_pos
        ##print x
        #for aa in range(dimension):
        #    for bb in range(num_atom):
        #        for cc in range(len(w1)):
        #            v1[bb*dimension + aa][cc] = v1[bb*dimension + aa][cc] * np.exp(1.j*G_cart[aa]*x[bb][aa])
        #            v1[bb*dimension + aa + len(w1)/2][cc] = v1[bb*dimension + aa + len(w1)/2][cc] * np.exp(1.j*G_cart[aa]*x[bb][aa])
        ########################################################################################      
        band_num = len(w1)
        #print band_num
        line = 'kx' + '\t' + 'ky' + '\t' + 'kz' + '\t' + str(q_vec[0]) + '\t' + str(q_vec[1]) + '\t' + str(q_vec[2]) + '\t' + str(q[0]) + '\t' + str(q[1]) + '\t' + str(q[2]) + '\n'
        bf.write(line)
        for j in range(band_num):
            line = str(j+1) + '\t' + str(w1[j]) + '\t'
            for k in range(band_num):
                line += str(v1[k][j]) + '\t'
            line += '\n'
            bf.write(line)
        bf.write('\n')

    bf.close()

    return 0    


def calculate_theta_at_fixed_ky(phband_file, ky, band_range):
    # band range: bands for calculting theta below an universal gap
    f = open(phband_file,'r')
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

    #print wf_data[24][4]


    ''' Second, calculate F and D matrix ''' 

    br = len(band_range)
    #print br

    D = np.identity(br, dtype=complex)
    H = np.identity(br, dtype=complex)

    #print D

    for i in range(mx):
        if i == mx-1:
            A_point = wf_data[i]
            B_point = wf_data[0]
        else:
            A_point = wf_data[i]
            B_point = wf_data[i+1]
        #print A_point
        #print B_point

        temp_F_matrix = get_F_matrix(A_point,B_point, band_range)
        #print temp_F_matrix
        temp_G_matrix = get_G_matrix(A_point,B_point, band_range)
        U, s, V = np.linalg.svd(temp_F_matrix, full_matrices=True)
        temp_F_matrix = np.dot(U, V)
        D = np.dot(D, temp_F_matrix)
        H = np.dot(H, temp_G_matrix)
        #print D

    #print D

    


    ''' Third, diagonalize D matrix and calculate theta at fixed ky'''

    Eigen = np.linalg.eigvals(D)
    #print Eigen

    #print np.imag(np.log(Eigen)) / (2*np.pi)
    

    A = np.imag(np.log(Eigen)) / (-2.0*np.pi)

    A = np.sort(A)
    #print A

    return ky_pos, A

def get_F_matrix(A_point, B_point, band_range):
    # A_point : kx,i   nb bands and wavefunctions
    # B_point : kx,i+1 nb bands and wavefunctions
    br = len(band_range)
    F = np.zeros((br, br), dtype=complex) 
    for i in range(br):
        for j in range(br):
            x = np.conjugate(A_point[band_range[i]])
            y = B_point[band_range[j]]
            F[i,j] = np.dot(x,np.transpose(y))
            #if i == j:
            #    print i, j
            #    print F[i,j]

    return F

def get_G_matrix(A_point, B_point, band_range):
    # A_point : kx,i   nb bands and wavefunctions
    # B_point : kx,i+1 nb bands and wavefunctions
    br = len(band_range)
    G = np.zeros((br, br), dtype=complex) 
    for i in range(br):
        for j in range(br):
            x = np.conjugate(A_point[band_range[i]])
            y = B_point[band_range[j]]
            if i == j:
                G[i,j] = np.dot(x,np.transpose(y))
            else:
                G[i,j] = 0.0
            #if i == j:
            #    print i, j
            #    print F[i,j]

    return G

def get_Willsons_loop(phband_file, band_range,q_grid):
    if q_grid[0] == 'slice':
        print 'Slice mode calculation'
        print '======================'
        print 'Information'
        print 'Nx, Ny, Nz = ' +str(q_grid[1][0]) + ' ' +str(q_grid[1][1]) + ' ' +str(q_grid[1][2]) + ' '
        print 'The q-point for fixed values = ' + str(q_grid[2])
        print '======================'
        calculate_all_theta_slicemode(phband_file, band_range,q_grid)
        return 0
    if q_grid[0] == 'line':
        print 'Line mode calculation (Not implemented yet)'
        #print '======================'
        #print 'Information'
        #print 'Nx, Ny, Nz = ' +str(q_grid[1][0]) + ' ' +str(q_grid[1][1]) + ' ' +str(q_grid[1][2]) + ' '
        #print 'The q-point for fixed values = ' + str(q_grid[2])
        #print '======================'
        #calculate_all_theta_linemode(phband_file, band_range,q_grid)
        return 0
    if q_grid[0] == 'node':
        print 'Node mode calculation'
        print '======================'
        print 'Information'
        print 'Center of node = ' +str(q_grid[1][0]) + ' ' +str(q_grid[1][1]) + ' ' +str(q_grid[1][2]) + ' '
        print 'Radius = ' + str(q_grid[2])
        print 'N_theta, N_phi = ' +str(q_grid[3]) + ' ' +str(q_grid[4]) + ' '
        print '======================' 
        calculate_all_theta_nodemode(phband_file, band_range,q_grid)
        return 0

    return 0    


def calculate_all_theta_nodemode(phband_file, band_range,q_grid):
    x = []
    y = []
    br = len(band_range)

    fig = plt.figure()

    for i in range(br):
        y.append([])

    theta_spacing = int(q_grid[3])
    delta_theta = np.pi / (theta_spacing-1)
    for i in range(theta_spacing):
        ky_pos, A = calculate_theta_at_fixed_ky(phband_file, i+1, band_range)
        x.append(delta_theta*i)
        for j in range(br):
            y[j].append(A[j])


    for i in range(br):
        plt.plot(x, y[i], linestyle="none", marker="o", ms=5, markeredgecolor='black', markerfacecolor='black')
    plt.title(str(band_range))
    plt.xlabel(r'$k$-point')
    plt.ylabel(r'Wannier Charge Center')
    plt.axis([-0.001 , np.pi + 0.001, -0.5, 0.5])
    fig.savefig('WCC.png')

    plt.show()
    
    return 0


def calculate_all_theta_slicemode(phband_file, band_range,q_grid):
    x = []
    y = []
    br = len(band_range)

    fig = plt.figure()

    for i in range(br):
        y.append([])

    for i in range(q_grid[1][1]):
        ky_pos, A = calculate_theta_at_fixed_ky(phband_file, i+1, band_range)
        x.append(ky_pos)
        for j in range(br):
            y[j].append(A[j])
        #line = str(ky_pos) + '\t'
        #for j in range(len(band_range)):
        #    line += str(A[j]) + '\t'
        #print line

    ky_pos, A = calculate_theta_at_fixed_ky(phband_file, 1, band_range)
    #print ky_pos
    #print A
    x.append(0.5)
    for j in range(br):
        y[j].append(A[j])
    #line = str(abs(ky_pos)) + '\t'
    #for j in range(len(band_range)):
    #    line += str(A[j]) + '\t'
    #print line

    for i in range(br):
        plt.plot(x, y[i], linestyle="none", marker="o", ms=5, markeredgecolor='black', markerfacecolor='black')
    plt.title(str(band_range))
    plt.xlabel(r'$k$-point')
    plt.ylabel(r'Wannier Charge Center')
    plt.axis([-0.5 -0.001 , 0.5 + 0.001, -0.5, 0.5])
    fig.savefig('WCC.png')

    plt.show()
    
    return 0



def find_nearest_supercell(super_latt_vec, super_atom_pos, ref_atom, target_atom):
    index_array = []
    dist_array = []
    for i in (-1, 0 , 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                index = np.array([i,j,k])
                ref = super_atom_pos[ref_atom]
                tar = super_atom_pos[target_atom] + index
                dist = tar - ref
                #dist = dist - np.round(tar-ref)
                dist = np.transpose(np.dot(np.transpose(super_latt_vec), np.transpose(dist)))
                dist = np.linalg.norm(dist)
                dist_array.append(dist)
                index_array.append(index)
    #print np.argmin(dist_array)
    #print index_array
    #print dist_array
    min_dist = dist_array[np.argmin(dist_array)]
    min_index = index_array[np.argmin(dist_array)]

    return min_dist, min_index

def extract_tb_parameter(primitive_file, supercell_file, FORCE_CONSTANTS_file, multiple, mass, dimension, symprec):
    super_latt_vec, super_atom_pos = read_info_cell(supercell_file)
    latt_vec, atom_pos = read_info_cell(primitive_file)
    fc2 = read_force_constants(FORCE_CONSTANTS_file)
    super_size = multiple[0] * multiple[1] * multiple[2]

    num_atom = len(atom_pos)
    vasp2THZ = 15.633302

    latt_vec = np.array(latt_vec)
    recip_fac =  (2*np.pi) / np.dot(latt_vec[0], np.cross(latt_vec[1], latt_vec[2]))
    recip_vec = np.array([np.cross(latt_vec[1], latt_vec[2]), np.cross(latt_vec[2], latt_vec[0]), np.cross(latt_vec[0], latt_vec[1])]) * recip_fac

    #print latt_vec

    #print atom_pos

    #print fc2[0][0][0][2]

    fc_tb_info_from_supercell = []

    g = open('fc_parameter','w')

    initial_line = 'primitive number' + '\t' + 'supercell_info' + '\t' + 'supercell number' + '\t' + 'distance' + '\t' + 'i' + '\t' + 'j' + '\t' + 'fc/mass' + '\n'
    #print initial_line
    g.write(initial_line)

    for i in range(num_atom):
        for j in range(num_atom):
            for x in range(multiple[0]):
                for y in range(multiple[1]):
                    for z in range(multiple[2]):
                        supercell_position = [x,y,z]
                        primitive_number = i*super_size
                        supercell_number = j*super_size + x  + y*multiple[0] + z*multiple[0]*multiple[1]
                        mass_term = sqrt(mass[i]*mass[j])
                        temp_info = [i, supercell_position, j, fc2[primitive_number][supercell_number]]
                        fc_tb_info_from_supercell.append(temp_info)
                        #force = []
                        #for k in range(3):
                        #    for l in range(3):
                        #        force.append(fc2[primitive_number][supercell_number][k][l])
                        #force = np.array(force) / mass_term
                        distance = np.dot(super_latt_vec.transpose(), (super_atom_pos[primitive_number]-super_atom_pos[supercell_number]).transpose())
                        distance = np.linalg.norm(distance)
                        #temp_line = str(i) + '\t' + str(supercell_position) + '\t' + str(j) + '\t' + str(distance) + '\t'
                        for k in range(dimension):
                            for l in range(dimension):
                                #temp_line = str(i) + '\t' + str(supercell_position) + '\t' + str(j) + '\t' + str(distance) + '\t' + str(k) + '\t' + str(l) + '\t' + str(fc2[primitive_number][supercell_number][k][l]/mass_term) + '\n'
                                temp_line = str(i) + '(' + str(k) + ')' + '\t' + str(supercell_position) + '\t' + str(j) + '(' + str(l) + ')' + '\t' + '==>' + '\t' + str(i*dimension+k) + '\t' + str(j*dimension+l) + '\t' + str(distance) + '\t' + str(fc2[primitive_number][supercell_number][k][l]/mass_term) + '\n'
                                g.write(temp_line)
                        #print temp_line 
                        g.write(temp_line)
    g.close()
    #print len(fc_tb_info_from_supercell)

    '''
    This subroutine will find other tb parameters from primitive cell to other adjacent cells
    '''

    #for i in range(len(fc_tb_info_from_supercell)):
    #for i in range(len(fc_tb_info_from_supercell)):
    #    primitive_atom = atom_pos[fc_tb_info_from_supercell[i][0]]
    #    supercell_atom = atom_pos[fc_tb_info_from_supercell[i][2]] + np.array(fc_tb_info_from_supercell[i][1])
    #    diff0 = supercell_atom - primitive_atom
    #    distance0 = np.linalg.norm(np.dot(diff0, latt_vec))
    #    #print distance0
    #    for x in range(-(multiple[0]-1), multiple[0]):
    #        for y in range(-(multiple[1]-1), multiple[1]):
    #            for z in range(-(multiple[2]-1), multiple[2]):
    #                if x < 0 or y < 0 or z < 0:
    #                    #print x, y, z
    #                    supercell_atom2 = atom_pos[fc_tb_info_from_supercell[i][2]] + np.array([x,y,z])
    #                    diff = supercell_atom2 - primitive_atom
    #                    distance = np.linalg.norm(np.dot(diff, latt_vec))
    #                    if abs(distance0 - distance) < symprec:
    #                        print x, y, z, fc_tb_info_from_supercell[i][0], fc_tb_info_from_supercell[i][2], fc_tb_info_from_supercell[i][1]#, fc_tb_info_from_supercell[i][3]

    return 0

def calculate_Berry_phase(phband_file, band_range):
    # band range: bands for calculting theta below an universal gap
    f = open(phband_file,'r')
    tempf = f.readlines()
    f.close()


    nkx, nb = int(tempf[0].split()[0]), int(tempf[0].split()[1])
    #print nkx, nb

    ''' First, read wavefunction and make proper array '''


    wf_data = [] # wf_data[ith kx point][nth ph band]

    for i in range(nkx):
        wf_data.append([])
        start_line = 2  + (1+nb+1)*i
        #print start_line
        for j in range(nb):
            tempf_wf = []
            for k in range(nb):
                tempf_wf.append(complex(tempf[start_line+j+1].split()[2+k]))
            wf_data[i].append(tempf_wf)

    #print wf_data[0][3]

    ''' Second, calculate F and D matrix ''' 

    br = len(band_range)
    F =  np.identity(br, dtype=complex)
    phi =  np.linalg.det(F)

    for i in range(nkx-1):
        temp = np.zeros((br, br), dtype=complex)
        for k in range(br):
            for j in range(br):
                A = np.conjugate(wf_data[i][band_range[k]])
                B = np.array(wf_data[i+1][band_range[j]])
                temp[k,j] = np.dot(A,B)
        phi = phi * np.linalg.det(temp)
 

    berry = -1 * np.imag(np.log(phi))

    print 'Berry phase is ' + str(berry)
    return 0

def calculate_Berry_curvature(filename, band_range, qx_range, qy_range):
    f = open(filename, 'r')
    tempf = f.readlines()
    f.close()

    g = open('berry.out', 'w')


    nqx = int(tempf[0].split()[0])
    nqy = int(tempf[0].split()[1])
    nb = int(tempf[0].split()[3])
    br = len(band_range)

    xvec = np.linspace(qx_range[0],qx_range[1],nqx)  
    yvec = np.linspace(qy_range[0],qy_range[1],nqy)  
    x, y = np.meshgrid(xvec, yvec)
    #print x
    F = np.zeros((nqx,nqy))

    Berry = 0


    for i in range(nqx):
        for j in range(nqy):
            line = (nb+2) * j + (nb+2) * (nqy+1) * i + 2
            q_cart = np.array([tempf[line].split()[k] for k in range(6,9)])
            #print q_cart
            line_xx = (nb+2) * j + (nb+2) * (nqy+1) * (i+1) + 2
            line_yy = (nb+2) * (j+1) + (nb+2) * (nqy+1) * i + 2
            line_xy = (nb+2) * (j+1) + (nb+2) * (nqy+1) * (i+1) + 2

            psi = [] ; psi_xx = [] ; psi_yy = [] ; psi_xy = []

            for m in range(br):
                #print tempf[line+band_range[m]+1]
                #print np.array([complex(tempf[line+band_range[m]+1].split()[k]) for k in range(2,2+nb)])
                psi.append(np.array([complex(tempf[line+band_range[m]+1].split()[k]) for k in range(2,2+nb)]))
                psi_xx.append(np.array([complex(tempf[line_xx+band_range[m]+1].split()[k]) for k in range(2,2+nb)]))
                psi_yy.append(np.array([complex(tempf[line_yy+band_range[m]+1].split()[k]) for k in range(2,2+nb)]))
                psi_xy.append(np.array([complex(tempf[line_xy+band_range[m]+1].split()[k]) for k in range(2,2+nb)]))
            #print len(psi)
            psi = np.array(psi).transpose() ; psi_xx = np.array(psi_xx).transpose() ; psi_yy = np.array(psi_yy).transpose() ; psi_xy = np.array(psi_xy).transpose()
            U_x_00 = np.dot(np.conjugate(psi).transpose(), psi_xx) ; U_x_00 = np.linalg.det(U_x_00) / np.linalg.norm(np.linalg.det(U_x_00))
            U_y_xx = np.dot(np.conjugate(psi_xx).transpose(), psi_xy) ; U_y_xx = np.linalg.det(U_y_xx) / np.linalg.norm(np.linalg.det(U_y_xx)) 
            U_x_yy = np.dot(np.conjugate(psi_yy).transpose(), psi_xy) ; U_x_yy = np.linalg.det(U_x_yy) / np.linalg.norm(np.linalg.det(U_x_yy)) 
            U_y_00 = np.dot(np.conjugate(psi).transpose(), psi_yy) ; U_y_00 = np.linalg.det(U_y_00) / np.linalg.norm(np.linalg.det(U_y_00))
            

            temp_F = np.imag(np.log((U_x_00*U_y_xx)/(U_x_yy*U_y_00)))
            Berry += temp_F
            #print np.imag(np.log(U_x_00)) + np.imag(np.log(U_y_xx)) - np.imag(np.log(U_x_yy)) - np.imag(np.log(U_y_00)), temp_F
            #print temp_F
            #print np.real(temp_F)
            temp_line = str(q_cart[0]) + '\t' + str(q_cart[1]) + '\t' + str(q_cart[2]) + '\t' + str(temp_F) + '\n'
            g.write(temp_line)
            #print temp_line
            F[i][j] = temp_F

    print Berry
    g.close()

    fig, ax = plt.subplots(figsize=(4*qx_range[1], 4*qy_range[1]))                     
    #ax.imshow(F.T)
    plt.contourf(x,y,F.T)
    #print F[0][29]                          
    plt.show()  

    return 0

def get_3Dplot_data(filename, band):
    vasp2THZ = 15.633302
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
    
    outname = '3Dplot_band_'+str(band)+'.out'
    g = open(outname, 'w')

    for i in range(ny):
        temp_line = ''
        for j in range(nx):
            temp_line += str(eigenval_rsh[i][j]) + '\t'
        temp_line += '\n'
        g.write(temp_line)

    g.close()



    return 0

def prepare_3Dplot_data(x_info, y_info):
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

def make_finite_along_one_direction(finite_direction, finite_length):
    '''
    current limitations: length > 2N +1 and odd number
    '''
    return 0

def read_band_yaml(filename, target_band_1, target_band_2):
    f = open(filename, 'r')
    tempf = f.readlines()
    f.close()

    nq = int(tempf[0].split()[1])

    for i in range(len(tempf)):
        if len(tempf[i].split()) != 0:
            if str(tempf[i].split()[0]) == 'natom:':
                atom = i
                break

    nb = 3*int(tempf[atom].split()[1])
    #print nb

    for i in range(len(tempf)):
        if len(tempf[i].split()) != 0:
            if str(tempf[i].split()[0]) == 'phonon:':
                start = i
                break
    
    q_pos = []
    target_1 = []
    target_2 = []
    diff = []

    for i in range(nq):
        templine = start + 1 + (4+nb*2)*i
        #print templine
        temp_q = np.array([str(tempf[templine].split()[j+3]) for j in range(3)])
        freq1 = float(tempf[templine+2+(target_band_1-1)*2+2].split()[1])
        freq2 = float(tempf[templine+2+(target_band_2-1)*2+2].split()[1])
        #print freq1, freq2
        q_pos.append(temp_q) ; target_1.append(freq1) ; target_2.append(freq2) ; diff.append(freq2-freq1)

    for i in range(nq):
        templine =  str(q_pos[i][0]) + ' ' + str(q_pos[i][1]) + ' ' + str(q_pos[i][2]) + ' ' + str(diff[i])
        print templine 
    #print nq

    return 0




mass1 = [40.078, 40.078, 40.078, 40.078, 14.007, 14.007, 10.81, 10.81, 10.81]
mass2 = [10.81, 10.81, 10.81]
mass3 = [10.81, 10.81, 10.81,10.81, 10.81, 10.81,10.81, 10.81, 10.81]
mass4 = [12.011,12.011,12.011,12.011,12.011,12.011,12.011,12.011,12.011,14.007,14.007,14.007,14.007,14.007,14.007]
mass5 = [12.011,12.011]
mass55 =[11.9, 12.1]
mass555 = [12.011, 12.011, 12.011, 12.011]
mass5555 = [14, 10, 10, 14]
mass6 = [30.973, 30.973, 30.973, 30.973]
mass7 = [28.085, 28.085]
mass77 = [28.085, 30.085]
mass8 = [28.085, 28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,28.085,1.0,1.0,1.0,1.0]
mass9 = [55.845, 55.845, 55.845, 55.845, 28.085, 28.085, 28.085, 28.085]
mass10 = [96.96, 127.6, 127.6]
mass11 = [58.933, 58.933, 58.933, 58.933, 28.085, 28.085, 28.085, 28.085]
mass12 = [96.96, 96.96, 127.6, 127.6, 127.6, 127.6]
mass1212 = [96.96, 960000.96, 127.6, 1270000.6, 127.6, 1270000.6]
mass13 = [96.96, 32.06, 32.06]
mass14 = [51.9961, 51.9961, 126.90, 126.90, 126.90, 126.90, 126.90, 126.90]
mass15 = [208.980,208.980,127.6,127.6,127.6]
mass16 = [96.96, 96.96, 78.96, 78.96, 78.96, 78.96]
mass1616 = [960000.96, 96.96, 78.96, 780000.96, 780000.96, 78.96]
mass17 = [96.96, 78.96, 78.96]
mass18 = [28.085, 28.085,28.085, 28.085, 195.084,195.084,195.084,195.084]
mass19 = [28.085, 28.085,28.085, 28.085, 106.42, 106.42, 106.42, 106.42]
mass20 = [47.867, 47.867, 16.0, 16.0, 16.0, 16.0]
mass21 = [69.723, 69.723, 14.0, 14.0]
mass22 = [65.38, 65.38, 65.38, 65.38, 65.38, 65.38, 65.38, 65.38, 121.76, 121.76, 121.76, 121.76, 121.76, 121.76, 121.76, 121.76]
mass23 = [58.6934, 58.6934, 58.6934, 58.6934, 30.973, 30.973, 30.973, 30.973, 30.973, 30.973, 30.973, 30.973]
mass24 = [96.96, 96.96, 32.06, 32.06, 32.06, 32.06]
mass25 = [26.981, 26.981, 26.981, 26.981, 26.981, 26.981, 26.981, 26.981, 26.981, 26.981, 26.981, 26.981, 55.845, 55.845, 55.845, 55.845, 28.085, 28.085, 28.085, 28.085, 28.085, 28.085, 28.085, 28.085]
mass26 = [40.078, 10.81, 10.81, 10.81, 10.81, 10.81, 10.81]
mass27 = [55.845, 55.845, 55.845, 55.845, 78.96, 78.96, 78.96, 78.96]
q_path1 = [[0, 0, 0], [0.333, 0.3333, 0.0], [0.0,0.0,0.0], [-0.3333, -0.3333, 0.0], [0, 0, 0]]
q_path2 = [[0, 0, 0], [0.5,0.0,0.0], [0.3333, 0.3333, 0.0], [0, 0, 0]]
q_path3 = [[-0.3333, -0.3333, 0], [0,0,0], [0.3333,0.3333,0]]
q_path4 = [[0, 0, 0], [0.5,0.0,0.0], [0.333, 0.3333, 0.0], [0.0,0.0,0.0], [-0.3333, -0.3333, 0.0], [-0.5,0.0,0.0], [0, 0, 0]]
q_path5 = [[0, 0, 0], [0.3333, 0.3333, 0.0], [0.6666, -0.3333, 0.0], [0, 0, 0]]
q_path6 = [[0, 0 ,0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0, 0, 0]]
q_path7 = [[0.333, 0.3333, 0.0],[-0.333, -0.3333, 0.0]]
q_path8 = [[0.0, 0.0 ,0.0 ], [0.5, 0.0 ,0.0] , [0.5, 0.5 ,0.0],  [0.0 ,0.5 ,0.0] , [0.0, 0.0 ,0.0]]
q_path9 = [[0.0 ,0.0, 0.0]  ,[0.000  , 0.500  , 0.000],  [0.500  , 0.500   ,0.000] ,   [0.0, 0.0 ,0.0], [0.500  ,0.500  , 0.500], [0.500 ,  0.500,   0.000]]
q_path10 = [[0.3333, 0.3333, 0.0], [0.2222, 0.2222, 0]]
q_path11 = [[0, 0, 0], [0.5,0.0,0.0], [0.3333, 0.3333, 0.0], [0, 0, 0], [0, 0, 0.5], [0.5,0.0,0.5], [0.3333, 0.3333, 0.5], [0, 0, 0.5]]
q_path12 = [[0.0, 0.0, 0.0], [0.5,0.0,0.0], [1.0,0.0,0.0]]
q_path13 = [[0.3333, 0.3333, 0.5], [0.5,0.0,0.5], [0.6666, -0.3333, 0.5]]
q_path14 = [[0.0, 0.0 ,0.0 ], [0.5, 0.0 ,0.0] , [0.5, 0.5 ,0.0],  [0.0 ,0.5 ,0.0] , [0.0, 0.0 ,0.0], [0.0, 0.0 ,0.5 ], [0.5, 0.0 ,0.5] , [0.5, 0.5 ,0.5],  [0.0 ,0.5 ,0.5] , [0.0, 0.0 ,0.5]]
q_path15 = [[0.0, 0.0 ,0.5 ], [0.5, 0. ,0.5] , [0.5, 0.5 ,0.5], [0.0 ,0.5 ,0.5], [0.0 ,0.0 ,0.5] ]
q_path16 = [[0.0, 0.0 ,0.5 ], [0., 0.5 ,0.5] , [0.5, 0.5 ,0.5],  [0.0 ,0.0 ,0.5] ]
q_path17 = [[0, 0, 0.5], [0.3333, 0.3333, 0.5], [0.6666, -0.3333, 0.5], [0, 0, 0.5]]
q_path18 = [[0.5, 0.0 ,0.5 ], [0.5, 0.5 ,0.5] , [0.0, 0.5 ,0.5],  [0.0 ,0.0 ,0.5]]
q_path19 = [[0.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.5,0.5,0.0] , [0.0,0.0,0.0]]
multiple = [4, 4, 1]
multiple3 = [4, 4, 1]
multiple4 = [2, 2, 1]
multiple5 = [4, 2, 1]
multiple6 = [2, 3, 1]
multiple7 = [2, 2, 2]
multiple8 = [3, 3, 1]
multiple9 = [4, 4, 2]
multiple10 = [3, 3, 3]

symprec = 1e-4
q_spacing = 20
dimension = 3

alpha0 = np.array([[0.0,0.0,0.0]])
alpha = np.array([0.0, 0.0, 0.01])
alpha_xy = np.array([0.03, 0.03, 0.0])
alpha_z = np.array([0.0, 0.0, 0.00])
alpha_MoSe2_bulk_Inv = np.array([alpha, -alpha, -alpha, alpha, alpha, -alpha])
alpha_MoSe2_bulk_noInv = np.array([alpha, alpha, alpha, alpha, alpha, alpha])
alpha_MoSe2_bulk_Inv_xy = np.array([alpha_xy, -alpha_xy, -alpha_xy, alpha_xy, alpha_xy, -alpha_xy])
alpha_MoSe2_bulk_noInv_xy = np.array([alpha_xy, alpha_xy, alpha_xy, alpha_xy, alpha_xy, alpha_xy])
alpha_MoSe2_bulk_Inv_z = np.array([alpha_z, -alpha_z, -alpha_z, alpha_z, alpha_z, -alpha_z])
alpha_MoSe2_bulk_noInv_z = np.array([alpha_z, alpha_z, alpha_z, alpha_z, alpha_z, alpha_z])
alpha_MoSe2_mono = np.array([alpha, alpha, alpha])
alpha_SiPd_bulk_Inv = np.array([-alpha, -alpha, alpha, alpha, alpha, alpha, -alpha, -alpha])
alpha_SiPd_bulk_noInv = np.array([alpha, alpha, alpha, alpha, alpha, alpha, alpha, alpha])
alpha_GaN_bulk_Inv = np.array([alpha, -alpha, -alpha, alpha])
alpha_graphene = np.array([alpha, alpha])

manipulate_fc_info0 = [False]
manipulate_fc_info1 = [True, [0,0,1,1],[1.0, 1.0, 1.0]]
manipulate_fc_info2 = [True, [0,1,0,1,0,1],[1.0, 1.0, 1.0]]
#extract_tb_parameter('POSCAR_Gra', 'SPOSCAR_Gra', 'FORCE_CONSTANTS_Gra', multiple4, mass5, 3, symprec)
#extract_tb_parameter('POSCAR_BP', 'SPOSCAR_BP', 'FORCE_CONSTANTS_BP', multiple5, mass6, 2)
#get_phonon_band(q_path2, 'POSCAR_B', 'SPOSCAR_B', 'FORCE_CONSTANTS_B', multiple, mass2, q_spacing, symprec, -0.2) 
#get_phonon_band(q_path2, 'POSCAR_SK', 'SPOSCAR_SK', 'FORCE_CONSTANTS_SK', multiple3, mass3, q_spacing, symprec, 0.2) 
#get_phonon_band(q_path2, 'POSCAR_C4N3', 'SPOSCAR_C4N3', 'FORCE_CONSTANTS_C4N3', multiple4, mass4, q_spacing, symprec, 0.2) 
#get_phonon_band('Graphene_noSPI', q_path5, 'POSCAR_Gra', 'SPOSCAR_Gra', 'FORCE_CONSTANTS_Gra', multiple4, mass5, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band(q_path5, 'POSCAR_2DSi', 'SPOSCAR_2DSi', 'FORCE_CONSTANTS_2DSi', multiple4, mass7, q_spacing, symprec, 0.1, dimension)
#get_phonon_band(q_path6, 'POSCAR_BP', 'SPOSCAR_BP', 'FORCE_CONSTANTS_BP', multiple5, mass6, q_spacing, symprec, 0.1, dimension)
#get_phonon_band(q_path8, 'POSCAR_GM', 'SPOSCAR_GM', 'FORCE_CONSTANTS_GM', multiple6, mass8, q_spacing, symprec, 0.0)
#get_phonon_band(q_path9, 'POSCAR_FeSi', 'SPOSCAR_FeSi', 'FORCE_CONSTANTS_FeSi', multiple7, mass9, q_spacing, symprec, 0.01, dimension)
#get_phonon_band(q_path9, 'POSCAR_CoSi', 'SPOSCAR_CoSi', 'FORCE_CONSTANTS_CoSi', multiple7, mass11, q_spacing, symprec, 0.00, dimension)
#get_phonon_band(q_path2, 'POSCAR_2H_MoS2', 'SPOSCAR_2H_MoS2', 'FORCE_CONSTANTS_2H_MoS2', multiple, mass13, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#get_phonon_band(q_path2, 'POSCAR_2H_MoTe2', 'SPOSCAR_2H_MoTe2', 'FORCE_CONSTANTS_2H_MoTe2', multiple, mass10, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_mono_SPI_z', q_path2, 'POSCAR_MoSe2', 'SPOSCAR_MoSe2', 'FORCE_CONSTANTS_MoSe2', multiple, mass17, q_spacing, symprec, alpha_MoSe2_mono, dimension, manipulate_fc_info0)
#get_phonon_band(q_path3, 'POSCAR_2H_MoTe2_Bi', 'SPOSCAR_2H_MoTe2_Bi', 'FORCE_CONSTANTS_2H_MoTe2_Bi', multiple, mass12, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#get_phonon_band(q_path2, 'POSCAR_BLG', 'SPOSCAR_BLG', 'FORCE_CONSTANTS_BLG', multiple4, mass555, q_spacing, symprec, alpha, dimension, manipulate_fc_info1)
#get_phonon_band(q_path2, 'POSCAR_CrI3', 'SPOSCAR_CrI3', 'FORCE_CONSTANTS_CrI3', multiple4, mass14, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#get_phonon_band(q_path2, 'POSCAR_CrI3_331', 'SPOSCAR_CrI3_331', 'FORCE_CONSTANTS_CrI3_331', multiple8, mass14, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#get_phonon_band(q_path2, 'POSCAR_Bi2Te3', 'SPOSCAR_Bi2Te3', 'FORCE_CONSTANTS_Bi2Te3', multiple, mass15, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#get_phonon_band('GaN_bulk_noSPI', q_path11, 'POSCAR_GaN', 'SPOSCAR_GaN', 'FORCE_CONSTANTS_GaN', multiple7, mass16, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band('GaN_bulk_SPI_Inv', q_path11, 'POSCAR_GaN', 'SPOSCAR_GaN', 'FORCE_CONSTANTS_GaN', multiple7, mass16, q_spacing, symprec, alpha_GaN_bulk_Inv, dimension, manipulate_fc_info0)
get_phonon_band('MoSe2_bulk_noSPI', q_path11, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulkAA_noSPI', q_path11, 'POSCAR_MoSe2_bulkAA', 'SPOSCAR_MoSe2_bulkAA', 'FORCE_CONSTANTS_MoSe2_bulkAA', multiple9, mass16, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulk_SPI_Inv', q_path11, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, q_spacing, symprec, alpha_MoSe2_bulk_Inv, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulk_SPI_noInv_003z', q_path11, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, q_spacing, symprec, alpha_MoSe2_bulk_noInv_z, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulk_SPI_Inv_003z', q_path11, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, q_spacing, symprec, alpha_MoSe2_bulk_Inv_z, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulkAA_SPI_Inv', q_path11, 'POSCAR_MoSe2_bulkAA', 'SPOSCAR_MoSe2_bulkAA', 'FORCE_CONSTANTS_MoSe2_bulkAA', multiple9, mass16, q_spacing, symprec, alpha_MoSe2_bulk_Inv, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulk_SPI_noInv', q_path11, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, q_spacing, symprec, alpha_MoSe2_bulk_noInv, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulk_SPI_noInv_003xy', q_path11, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, q_spacing, symprec, alpha_MoSe2_bulk_noInv_xy, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulk_SPI_Inv_005xy', q_path11, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, q_spacing, symprec, alpha_MoSe2_bulk_Inv_xy, dimension, manipulate_fc_info0)
#get_phonon_band('MoSe2_bulkAA_SPI_noInv', q_path11, 'POSCAR_MoSe2_bulkAA', 'SPOSCAR_MoSe2_bulkAA', 'FORCE_CONSTANTS_MoSe2_bulkAA', multiple9, mass16, q_spacing, symprec, alpha_MoSe2_bulk_noInv, dimension, manipulate_fc_info0)
#get_phonon_band(q_path13, 'POSCAR_MoSe2', 'SPOSCAR_MoSe2', 'FORCE_CONSTANTS_MoSe2', multiple, mass17, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#get_phonon_band(q_path14, 'POSCAR_SiPt', 'SPOSCAR_SiPt', 'FORCE_CONSTANTS_SiPt', multiple7, mass18, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#get_phonon_band('AlFeSi_bulk_noSPI', q_path14, 'POSCAR_ASF', 'SPOSCAR_ASF', 'FORCE_CONSTANTS_ASF', multiple7, mass25, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band('SiPt_lowsym_bulk_noSPI', q_path14, 'POSCAR_SiPt_lowsym', 'SPOSCAR_SiPt_lowsym', 'FORCE_CONSTANTS_SiPt_lowsym', multiple7, mass18, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band('SiPd_bulk_noSPI', q_path18, 'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band('SiPd_bulk_noSPI', q_path14, 'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band('SiPd_bulk_SPI_Inv', q_path14,   'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, q_spacing, symprec, alpha_SiPd_bulk_Inv, dimension, manipulate_fc_info0)
#get_phonon_band('SiPd_bulk_SPI_noInv', q_path14, 'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, q_spacing, symprec, alpha_SiPd_bulk_noInv, dimension, manipulate_fc_info0)
#get_phonon_band('CaB6_bulk_noSPI', q_path19, 'POSCAR_CaB', 'SPOSCAR_CaB', 'FORCE_CONSTANTS_CaB', multiple7, mass26, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band('FeSe_mono_AFMcol_noSPI', q_path19, 'POSCAR_FeSe_mono_AFMcol', 'SPOSCAR_FeSe_mono_AFMcol', 'FORCE_CONSTANTS_FeSe_mono_AFMcol', multiple4, mass27, q_spacing, symprec, alpha0, dimension, manipulate_fc_info0)
#get_phonon_band(q_path16, 'POSCAR_TiO2', 'SPOSCAR_TiO2', 'FORCE_CONSTANTS_TiO2', multiple7, mass20, q_spacing, symprec, alpha, dimension, manipulate_fc_info0)
#draw_phonon_band('ph_frequecny_SiPd_bulk_noSPI.out')
#draw_phonon_band('ph_frequecny_MoSe2_bulk_SPI_noInv_005z.out')
#draw_phonon_projected_band('ph_frequecny_MoSe2_bulk_SPI_Inv_003z_projected.out', [0,3,4], [1,2,5])
#draw_phonon_projected_band('ph_frequecny_GaN_bulk_noSPI_projected.out', [0,3], [1,2])

#draw_phonon_projected_band('ph_frequecny_MoSe2_bulk_SPI_Inv_projected.out', [0,3,4], [1,2,5])
#draw_phonon_projected_band('ph_frequecny_GaN_bulk_SPI_Inv_projected.out', [0,3], [1,2])

#draw_phonon_projected_band('ph_frequecny_SiPd_bulk_noSPI_projected.out', [0,1,6,7], [2,3,4,5])
#draw_phonon_band('ph_frequecny_SiPd_bulk_SPI_Inv.out')

#read_band_yaml('band.yaml', 9,10)

q_path_berry_CaB6_nl_91011 = [[0.0+0.001,0.3921508,0.0-0.001],[0.0+0.001,0.3921508,0.0+0.001],[0.0-0.001,0.3921508,0.0+0.001],[0.0-0.001,0.3921508,0.0-0.001],[0.5+0.001,0.3921508, 0.0-0.001]]
q_path_berry_CaB6_nl_678 = [[0.0+0.001,0.4378863,0.0-0.001],[0.0+0.001,0.4378863,0.0+0.001],[0.0-0.001,0.4378863,0.0+0.001],[0.0-0.001,0.4378863,0.0-0.001],[0.5+0.001,0.4378863, 0.0-0.001]]
q_path_berry_CaB6_nl_910 = [[0.2385661+0.001,0.2385661,0.0-0.001],[0.2385661+0.001,0.2385661,0.0+0.001],[0.2385661-0.001,0.2385661,0.0+0.001],[0.2385661-0.001,0.2385661,0.0-0.001],[0.2385661+0.001,0.2385661, 0.0-0.001]]
q_path_berry_SiPd = [[0.5+0.001,0.5-0.001,0.5],[0.5+0.001,0.5+0.001,0.5],[0.5-0.001,0.5+0.001,0.5],[0.5-0.001,0.5-0.001,0.5],[0.5+0.001,0.5-0.001, 0.5]]
q_path_berry_Graphene = [[0.333+0.001,0.333-0.001,0.0],[0.333+0.001,0.333+0.001,0.0],[0.333-0.001,0.333+0.001,0.0],[0.333-0.001,0.333-0.001,0.0],[0.333+0.001,0.333-0.001,0.0]]
q_path_berry_MoSe2 = [[0.0+0.001,0.0-0.001,0.5],[0.0+0.001,0.0+0.001,0.5],[0.0-0.001,0.0+0.001,0.5],[0.0-0.001,0.0-0.001,0.5],[0.0+0.001,0.0-0.001,0.5]]
q_path_berry_ZnSb_z = [[0.5+0.001,0.5-0.001,0.3],[0.5+0.001,0.5+0.001,0.3],[0.5-0.001,0.5+0.001,0.3],[0.5-0.001,0.5-0.001,0.3],[0.5+0.001,0.5-0.001,0.3]]
q_path_berry_ZnSb_x = [[0.3,0.5-0.001,0.5+0.001],[0.3,0.5+0.001,0.5+0.001],[0.3,0.5+0.001,0.5-0.001],[0.3,0.5-0.001,0.5-0.001],[0.3,0.5-0.001,0.5+0.001]]
q_path_berry_ZnSb_y = [[0.5+0.001,0.5-0.001,0.3],[0.5+0.001,0.5+0.001,0.3],[0.5-0.001,0.5+0.001,0.3],[0.5-0.001,0.5-0.001,0.3],[0.5+0.001,0.5-0.001,0.3]]
q_path_berry_SiPt = [[0.5+0.001,0.5-0.001,0.5],[0.5+0.001,0.5+0.001,0.5],[0.5-0.001,0.5+0.001,0.5],[0.5-0.001,0.5-0.001,0.5],[0.5+0.001,0.5-0.001,0.5]]
q_path_berry_NiP2 = [[0.5+0.001,0.5-0.001,0.3],[0.5+0.001,0.5+0.001,0.3],[0.5-0.001,0.5+0.001,0.3],[0.5-0.001,0.5-0.001,0.3],[0.5+0.001,0.5-0.001,0.3]]
q_path_berry_MoS2_mono_nl = [[0.1766027+0.005,0.1766027,0.0-0.005],[0.1766027+0.005,0.1766027,0.0+0.005],[0.1766027-0.005,0.1766027,0.0+0.005],[0.1766027-0.005,0.1766027,0.0-0.005],[0.1766027+0.005,0.1766027,0.0-0.005]]
q_path_berry_MoS2_mono_nl2 = [[0.2738924+0.005,0.2738924,0.0-0.005],[0.2738924+0.005,0.2738924,0.0+0.005],[0.2738924-0.005,0.2738924,0.0+0.005],[0.2738924-0.005,0.2738924,0.0-0.005],[0.2738924+0.005,0.2738924,0.0-0.005]]
q_path_berry_MoSe2_mono_nl = [[0.2493426+0.001,0.2493426,0.0-0.001],[0.2493426+0.001,0.2493426,0.0+0.001],[0.2493426-0.001,0.2493426,0.0+0.001],[0.2493426-0.001,0.2493426,0.0-0.001],[0.2493426+0.001,0.2493426,0.0-0.001]]
q_path_berry_MoSe2_mono = [[0.0+0.001,0.0-0.001,0.0],[0.0+0.001,0.0+0.001,0.0],[0.0-0.001,0.0+0.001,0.0],[0.0-0.001,0.0-0.001,0.0],[0.0+0.001,0.0-0.001,0.0]]
q_path_berry_MoSe2_nl_14_17 = [[0.2473151+0.001,0.2473151,0.0-0.001],[0.2473151+0.001,0.2473151,0.0+0.001],[0.2473151-0.001,0.2473151,0.0+0.001],[0.2473151-0.001,0.2473151,0.0-0.001],[0.2473151+0.001,0.2473151,0.0-0.001]]
q_path_berry_MoSe2_nl_15_17 = [[0.2455780+0.001,0.2455780,0.0-0.001],[0.2455780+0.001,0.2455780,0.0+0.001],[0.2455780-0.001,0.2455780,0.0+0.001],[0.2455780-0.001,0.2455780,0.0-0.001],[0.2455780+0.001,0.2455780,0.0-0.001]]
q_path_berry_MoSe2_nl_14_16 = [[0.2494864+0.001,0.2494864,0.0-0.001],[0.2494864+0.001,0.2494864,0.0+0.001],[0.2494864-0.001,0.2494864,0.0+0.001],[0.2494864-0.001,0.2494864,0.0-0.001],[0.2494864+0.001,0.2494864,0.0-0.001]]
q_path_berry_MoSe2_nl_15_16 = [[0.2452523+0.001,0.2452523,0.0-0.001],[0.2452523+0.001,0.2452523,0.0+0.001],[0.2452523-0.001,0.2452523,0.0+0.001],[0.2452523-0.001,0.2452523,0.0-0.001],[0.2452523+0.001,0.2452523,0.0-0.001]]
#q_grid = ['slice',[31, 31, 1], 0.5]  #### [q_slice mode, [nx, ny, nz], fixed_qpoints]
#q_grid = ['line', q_path7, 20, 2, 51 ]  #### [[fixed_qpoint_direction, fixed_qpoint_value], [mxn for 2D planes]]
#q_grid = ['node', [0.5,0.5,0.5], 0.005, 21, 21]  #### [q_node mode, center_direct_coord, radius from the node. theta_spacing, phi_spacing]
q_grid = ['berryphase', q_path_berry_MoSe2, 10]
dimension = 3
#make_phband_PROCAR_format(q_grid, 'POSCAR_BP', 'SPOSCAR_BP', 'FORCE_CONSTANTS_BP', multiple5, mass6, symprec, dimension, 0.04)
#make_phband_PROCAR_format('Graphene_noSPI', q_grid, 'POSCAR_Gra', 'SPOSCAR_Gra', 'FORCE_CONSTANTS_Gra', multiple4, mass5, symprec, dimension, alpha0)
#make_phband_PROCAR_format('Graphene_SPI_0.2z', q_grid, 'POSCAR_Gra', 'SPOSCAR_Gra', 'FORCE_CONSTANTS_Gra', multiple4, mass55, symprec, dimension, alpha_graphene)
#make_phband_PROCAR_format(q_grid, 'POSCAR_2DSi', 'SPOSCAR_2DSi', 'FORCE_CONSTANTS_2DSi', multiple4, mass7, symprec, dimension, 0.1)
#make_phband_PROCAR_format(q_grid, 'POSCAR_SK', 'SPOSCAR_SK', 'FORCE_CONSTANTS_SK', multiple3, mass3, symprec, dimension, 0.2)
#make_phband_PROCAR_format(q_grid, 'POSCAR_C4N3', 'SPOSCAR_C4N3', 'FORCE_CONSTANTS_C4N3', multiple4, mass4, symprec, dimension, 0.2)
#make_phband_PROCAR_format(q_grid, 'POSCAR_FeSi', 'SPOSCAR_FeSi', 'FORCE_CONSTANTS_FeSi', multiple7, mass9, symprec, dimension, 0.0)
#make_phband_PROCAR_format(q_grid, 'POSCAR_CoSi', 'SPOSCAR_CoSi', 'FORCE_CONSTANTS_CoSi', multiple7, mass11, symprec, dimension, 0.0)
#make_phband_PROCAR_format(q_grid, 'POSCAR_2H_MoTe2', 'SPOSCAR_2H_MoTe2', 'FORCE_CONSTANTS_2H_MoTe2', multiple, mass10, symprec, dimension, alpha)
#make_phband_PROCAR_format(q_grid, 'POSCAR_2H_MoTe2_Bi', 'SPOSCAR_2H_MoTe2_Bi', 'FORCE_CONSTANTS_2H_MoTe2_Bi', multiple, mass12, symprec, dimension, alpha)
#make_phband_PROCAR_format(q_grid, 'POSCAR_BLG', 'SPOSCAR_BLG', 'FORCE_CONSTANTS_BLG', multiple4, mass555, symprec, dimension, [0.0,0.0,0.2])
#make_phband_PROCAR_format(q_grid, 'POSCAR_CrI3', 'SPOSCAR_CrI3', 'FORCE_CONSTANTS_CrI3', multiple4, mass14, symprec, dimension, alpha)
#make_phband_PROCAR_format(q_grid, 'POSCAR_Bi2Te3', 'SPOSCAR_Bi2Te3', 'FORCE_CONSTANTS_Bi2Te3', multiple, mass15, symprec, dimension, alpha)
#make_phband_PROCAR_format(q_grid, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, symprec, dimension, alpha)
#make_phband_PROCAR_format('MoS2_mono_noSPI', q_grid, 'POSCAR_2H_MoS2', 'SPOSCAR_2H_MoS2', 'FORCE_CONSTANTS_2H_MoS2', multiple, mass13, symprec, dimension, alpha0)
#make_phband_PROCAR_format('MoSe2_mono_noSPI', q_grid, 'POSCAR_MoSe2', 'SPOSCAR_MoSe2', 'FORCE_CONSTANTS_MoSe2', multiple, mass17, symprec, dimension, alpha0)
#make_phband_PROCAR_format('MoSe2_mono_SPI_0015z', q_grid, 'POSCAR_MoSe2', 'SPOSCAR_MoSe2', 'FORCE_CONSTANTS_MoSe2', multiple, mass17, symprec, dimension, alpha_MoSe2_mono)
#make_phband_PROCAR_format('MoSe2_mono_SPI_xy', q_grid, 'POSCAR_MoSe2', 'SPOSCAR_MoSe2', 'FORCE_CONSTANTS_MoSe2', multiple, mass17, symprec, dimension, alpha_MoSe2_mono)
#make_phband_PROCAR_format('MoSe2_bulkAA_SPI_noInv', q_grid, 'POSCAR_MoSe2_bulkAA', 'SPOSCAR_MoSe2_bulkAA', 'FORCE_CONSTANTS_MoSe2_bulkAA', multiple9, mass16, symprec, dimension, alpha_MoSe2_bulk_noInv)
#make_phband_PROCAR_format('MoSe2_bulk_SPI_Inv', q_grid, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, symprec, dimension, alpha_MoSe2_bulk_Inv)
#make_phband_PROCAR_format('MoSe2_bulk_SPI_noInv_003xy_kzpi', q_grid, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, symprec, dimension, alpha_MoSe2_bulk_noInv_xy)
#make_phband_PROCAR_format('MoSe2_bulk_SPI_Inv_005xy_kzpi', q_grid, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, symprec, dimension, alpha_MoSe2_bulk_Inv_xy)
#make_phband_PROCAR_format('MoSe2_bulk_SPI_Inv_003z_kzpi', q_grid, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, symprec, dimension, alpha_MoSe2_bulk_Inv_z)
#make_phband_PROCAR_format('MoSe2_bulk_SPI_noInv_kz0', q_grid, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, symprec, dimension, alpha_MoSe2_bulk_noInv)
#make_phband_PROCAR_format('MoS2_bulk_noSPI_kz5', q_grid, 'POSCAR_2H_MoS2_bulk', 'SPOSCAR_2H_MoS2_bulk', 'FORCE_CONSTANTS_2H_MoS2_bulk', multiple9, mass24, symprec, dimension, alpha0)
#make_phband_PROCAR_format('MoSe2_bulk_noSPI', q_grid, 'POSCAR_MoSe2_bulk', 'SPOSCAR_MoSe2_bulk', 'FORCE_CONSTANTS_MoSe2_bulk', multiple9, mass16, symprec, dimension, alpha0)
#make_phband_PROCAR_format('SiPt_lowsym_bulk_noSPI', q_grid, 'POSCAR_SiPt_lowsym', 'SPOSCAR_SiPt_lowsym', 'FORCE_CONSTANTS_SiPt_lowsym', multiple7, mass18, symprec, dimension, alpha0)
#make_phband_PROCAR_format('NiP2_bulk_noSPI', q_grid, 'POSCAR_NiP2', 'SPOSCAR_NiP2', 'FORCE_CONSTANTS_NiP2', multiple7, mass23, symprec, dimension, alpha0)
#make_phband_PROCAR_format('ZnSb_bulk_noSPI', q_grid, 'POSCAR_ZnSb', 'SPOSCAR_ZnSb', 'FORCE_CONSTANTS_ZnSb', multiple7, mass22, symprec, dimension, alpha0)
#make_phband_PROCAR_format('SiPd_bulk_noSPI', q_grid, 'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, symprec, dimension, alpha0)
#make_phband_PROCAR_format('SiPd_bulk_SPI_Inv', q_grid, 'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, symprec, dimension, alpha_SiPd_bulk_Inv)
#make_phband_PROCAR_format('SiPd_bulk_SPI_noInv', q_grid, 'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, symprec, dimension, alpha_SiPd_bulk_noInv)
#make_phband_PROCAR_format('FeSi_bulk_SPI_noInv', q_grid, 'POSCAR_FeSi', 'SPOSCAR_FeSi', 'FORCE_CONSTANTS_FeSi', multiple7, mass9, symprec, dimension, alpha0)
#make_phband_PROCAR_format('FeSi_bulk_333_SPI_noInv', q_grid, 'POSCAR_FeSi_333', 'SPOSCAR_FeSi_333', 'FORCE_CONSTANTS_FeSi_333', multiple10, mass9, symprec, dimension, alpha0)
#make_phband_PROCAR_format('CaB6_bulk_noSPI', q_grid, 'POSCAR_CaB', 'SPOSCAR_CaB', 'FORCE_CONSTANTS_CaB', multiple7, mass19, symprec, dimension, alpha0)
#make_phband_PROCAR_format(q_grid, 'POSCAR_SiPt', 'SPOSCAR_SiPt', 'FORCE_CONSTANTS_SiPt', multiple7, mass18, symprec, dimension, alpha)
#make_phband_PROCAR_format(q_grid, 'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, symprec, dimension, alpha)
#print calculate_theta_at_fixed_ky('phband_PROCAR.out', 1, band_range)


band_range = [int(i) for i in range(18,30)] ; print '# of bands = ' + str(len(band_range)) + ' Detail bands = ' + str(band_range)
#band_range = [44,4]
#get_Willsons_loop('phband_PROCAR_Graphene_SPI_0.2z_slicemode.out', band_range, q_grid)
#get_Willsons_loop('phband_PROCAR_MoS2_bulk_noSPI_kz5_slicemode.out', band_range, q_grid)
#get_Willsons_loop('phband_PROCAR_MoSe2_bulk_SPI_Inv_003z_kzpi_slicemode.out', band_range, q_grid)
#get_Willsons_loop('phband_PROCAR_MoSe2_bulk_SPI_noInv_001xy_kzpi_slicemode.out', band_range, q_grid)
#get_Willsons_loop('phband_PROCAR_MoSe2_mono_SPI_0015z_slicemode.out', band_range, q_grid)
#get_Willsons_loop('phband_PROCAR_SiPd_bulk_noSPI_nodemode.out', band_range, q_grid)
#get_Willsons_loop('phband_PROCAR_FeSi_bulk_333_SPI_noInv_nodemode.out', band_range, q_grid)
#get_Willsons_loop('phband_PROCAR_slicemode.out', band_range, q_grid)
#get_Willsons_loop('phband_PROCAR_SiPt_lowsym_bulk_noSPI_nodemode.out', band_range, q_grid)
#calculate_Berry_phase('phband_PROCAR_SiPd_bulk_noSPI_berryphase.out', band_range)
#calculate_Berry_phase('phband_PROCAR_Graphene_noSPI_berryphase.out', band_range)
#calculate_Berry_phase('phband_PROCAR_MoSe2_bulk_noSPI_berryphase.out', band_range)
#calculate_Berry_phase('phband_PROCAR_ZnSb_bulk_noSPI_berryphase.out', band_range)
#calculate_Berry_phase('phband_PROCAR_SiPt_lowsym_bulk_noSPI_berryphase.out', band_range)
#calculate_Berry_phase('phband_PROCAR_NiP2_bulk_noSPI_berryphase.out', band_range)
#calculate_Berry_phase('phband_PROCAR_MoSe2_mono_noSPI_berryphase.out', band_range)
#calculate_Berry_phase('phband_PROCAR_MoS2_mono_noSPI_berryphase.out', band_range)
#calculate_Berry_phase('phband_PROCAR_CaB6_bulk_noSPI_berryphase.out', band_range)
#calculate_all_theta_nodemode('phband_PROCAR_nodemode.out', band_range, q_grid)



latt_vec, atom_pos = read_info_cell('POSCAR_SiPd')
recip_fac =  (2*np.pi) / np.dot(latt_vec[0], np.cross(latt_vec[1], latt_vec[2]))
recip_vec = np.array([np.cross(latt_vec[1], latt_vec[2]), np.cross(latt_vec[2], latt_vec[0]), np.cross(latt_vec[0], latt_vec[1])]) * recip_fac
print recip_vec
q_grid = ['berry', [0, 1.8262], [0, 1.111], [2, 0.5], 101, 101]
dimension =3
alpha = [0.0,0.0,0.0]
manipulate_fc_info0 = [False]
band_range = [int(i) for i in range(15,16)] ; print '# of bands = ' + str(len(band_range)) + ' Detail bands = ' + str(band_range)
#make_phband_PROCAR_format(q_grid, 'POSCAR_Gra', 'SPOSCAR_Gra', 'FORCE_CONSTANTS_Gra', multiple4, mass55, symprec, dimension, [0.0,0.0,0,0])
#make_phband_PROCAR_format(q_grid, 'POSCAR_2DSi', 'SPOSCAR_2DSi', 'FORCE_CONSTANTS_2DSi', multiple4, mass7, symprec, dimension, 0.1)
#make_phband_PROCAR_format(q_grid, 'POSCAR_2H_MoS2', 'SPOSCAR_2H_MoS2', 'FORCE_CONSTANTS_2H_MoS2', multiple, mass13, symprec, dimension, alpha)
#make_phband_PROCAR_format(q_grid, 'POSCAR_2H_MoTe2', 'SPOSCAR_2H_MoTe2', 'FORCE_CONSTANTS_2H_MoTe2', multiple, mass10, symprec, dimension, alpha)
#make_phband_PROCAR_format(q_grid, 'POSCAR_MoSe2', 'SPOSCAR_MoSe2', 'FORCE_CONSTANTS_MoSe2', multiple, mass17, symprec, dimension, alpha)
#make_phband_PROCAR_format('SiPd_bulk_noSPI_kzpi', q_grid, 'POSCAR_SiPd', 'SPOSCAR_SiPd', 'FORCE_CONSTANTS_SiPd', multiple7, mass19, symprec, dimension, alpha0)
#calculate_Berry_curvature('phband_PROCAR_berrycurv.out', band_range, q_grid[1], q_grid[2])
#get_3Dplot_data('phband_PROCAR_SiPd_bulk_noSPI_kzpi_berrycurv.out', 44)
#prepare_3Dplot_data([q_grid[4]+1, q_grid[1][0], (q_grid[1][1]-q_grid[1][0])/(q_grid[4]-1)], [q_grid[5]+1, q_grid[2][0],  (q_grid[2][1]-q_grid[2][0])/(q_grid[5]-1)])

#super_latt_vec, super_atom_pos = read_info_cell('SPOSCAR')

#print len(atom_pos)
#fc2 = read_force_constants('FORCE_CONSTANTS')
#mass = [40.078, 40.078, 40.078, 40.078, 14.007, 14.007, 10.81, 10.81, 10.81]

#super_size = multiple[0] * multiple[1] * multiple[2]

#super_size=16
#q_vec = [0.3333,0.3333,0.0]



#vecs, multiplicity = get_shortest_vectors(super_latt_vec, super_atom_pos, atom_pos, super_size, symprec)
#print multiplicity[15][0], np.linalg.norm(vecs[15][0])
#print find_atom_number_in_supercell(atom_pos,1,multiple)


#dyn = construct_dynamical_matrix_q(fc2, super_latt_vec, super_atom_pos, atom_pos, latt_vec, q_vec, multiple, mass, symprec)
#w = np.linalg.eigvalsh(dyn).real
#print np.sqrt(np.abs(w)) * 15.633302
#print np.sort(np.sqrt(w))
#print len(dyn)


#print refine_frequency(np.array([-1.0, -4.0, 4.0, 9.0, -9.0,-2.0]))


#print fc2[1][10][0][2]
#print find_atom_number_in_supercell(atom_pos, 2, multiple)
#get_equivalent_smallest_vectors(3, 0, 16, super_latt_vec, super_atom_pos, latt_vec, 1e-5)
#print len(get_shortest_vectors(super_latt_vec, super_atom_pos, atom_pos, super_size, symprec)[0])


#select_specific_atoms_from_phonopy('FORCE_CONSTANTS',97, 144)

