#!/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Gechuanqi Pan'

import os
import time
import glob
import pymatgen as pmg
from pymatgen.io import cif

global Periodic_table_of_elements
Periodic_table_of_elements = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm','Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo')


def get_info_from_PIM_file(molten_salt_system):
	"""
		Description:
			get info from PIM_paramters_file, the result is stored in list

		Args:
			molten_salt_system: the simulated system

		Returns:
			return a list
	"""

	parameters_filename = '/WORK/nscc-gz_material_1/pangchq/parameter_files/PIM_parameters_file'
	with open(parameters_filename, 'r', encoding='UTF-8') as f_PIM_para:
		# firstly, judge if system is in PIM_parameters_file, if false, report an error and exit
		flag_exist = False
		for num, line in enumerate(f_PIM_para):
			if line.strip() == ('#'+molten_salt_system):
				flag_exist = True
				break
		if flag_exist == False:
			print(molten_salt_system+" is not in PIM_parameters_file, please check!")
			exit()
		f_PIM_para.seek(0, 0)

		PIM_para_list = []
		flag = False
		for num, line in enumerate(f_PIM_para):
			if line.strip() == ('#'+molten_salt_system):
				flag = True
				print('begin reading force field parameters...')
				continue
			if line.strip() == ('#END '+molten_salt_system):
				print('end reading force field parameters...')
				break
			if flag and (not line.startswith('#')):
				PIM_para_list.append(list(filter(None,line.strip().split(' '))))

	return PIM_para_list


def gen_directory(molten_salt_system, PIM_para_list):
	"""
		Description:
			create subdirectories by target temperature for molten_salt_system and the result is stored in list, e.g. ['900K', '950K', '1000K', ..., '1650K', 'gen_liquid_model']

		Args:
			molten_salt_system: the simulated system

		Returns:
			return a list
	"""

	Ncount = 0
	Lower_temperature_limit = 0
	Higher_temperature_limit = 0
	Temperature_interval = 0
	for i in PIM_para_list:
		if i[0] == 'Lower_temperature_limit':
			Lower_temperature_limit = int(i[1])
		if i[0] == 'Higher_temperature_limit':
			Higher_temperature_limit = int(i[1])
		if i[0] == 'Temperature_interval':
			Temperature_interval = int(i[1])
			break
	Ncount = 1 + (Higher_temperature_limit - Lower_temperature_limit)//Temperature_interval
	Temperature_list = []
	for n in range(Ncount):
		Temperature_list.append(Lower_temperature_limit + Temperature_interval*n)
	Temperature_list = [str(T)+'K' for T in Temperature_list]
	Temperature_list.append('gen_liquid_model')
	# create subdirectories
	dir_current = os.getcwd()
	for directory in Temperature_list:
		dir_next = os.path.join(dir_current, directory)
		if not os.path.exists(dir_next):
			os.mkdir(dir_next)

	return Temperature_list


def get_info_from_ICSD_file(molten_salt_system):
	"""
		Description:
			get info from ICSD_paramters_file, the result is stored in a list

		Args:
			molten_salt_system: the simulated system

		Returns:
			return a list
	"""

	parameters_filename = '/WORK/nscc-gz_material_1/pangchq/parameter_files/ICSD_parameters_file'
	with open(parameters_filename, 'r', encoding='UTF-8') as f_ICSD_para:
		# firstly, judge if system is in ICSD_parameters_file, if false, report an error and exit
		flag_exist = False
		for num, line in enumerate(f_ICSD_para):
			if line.strip() == ('#'+molten_salt_system): 
				flag_exist = True
				break
		if flag_exist == False:
			print(molten_salt_system+" is not in ICSD_parameters_file, please check!")
			exit()
		f_ICSD_para.seek(0, 0)

		ICSD_para_list = []
		flag = False
		for num, line in enumerate(f_ICSD_para):
			if line.strip() == ('#'+molten_salt_system):
				flag = True
				print('begin reading structure parameters...')
				continue
			if line.strip() == ('#END '+molten_salt_system):
				print('end reading structure parameters...')
				break
			if flag and (not line.startswith('#')):
				ICSD_para_list.append(list(filter(None,line.strip('\n').split(' '))))

	return ICSD_para_list


def get_structure_parameters(ICSD_para_list):
	"""
		Description:
			get spacegroup、unit_cell_list、element_list、fraction_coords_list from ICSD_para_list for building model by pymatgen

		Args:
			ICSD_para_list: the list contains needed info exacted from ICSD_paramters_file

		Returns:
			spacegroup: str
			unit_cell_list: list, 1d
			element_list: list, 1d
			fraction_coords_list: list, 2d
			Z_number: int
	"""

	atoms_list = []
	unit_cell_list = []
	for i in ICSD_para_list:
		if i[0] == 'SG':
			spacegroup = i[2]
		if i[0] == 'Unit':
			for n in range(2, len(i)):
				unit_cell_list.append(i[n])
		if i[0] == 'Z':
			Z_number = int(i[1])
		if i[0] in Periodic_table_of_elements:
			atoms_list.append(i[0:1]+i[2:3]+i[5:8]) # Atom, OX, x, y, z
	
	# delete (x) and ., for example, 3.816(1) 11.815(2) 8.507(1) 90. 90. 90. to 3.816 11.815 8.507 90 90 90
	for i in range(len(unit_cell_list)):
		if unit_cell_list[i].find('(') != -1:
			unit_cell_list[i] = unit_cell_list[i][0:unit_cell_list[i].find('(')]
		if unit_cell_list[i].endswith('.'):
			unit_cell_list[i] = unit_cell_list[i].strip('.')
		unit_cell_list[i] = float(unit_cell_list[i])

	element_list = []
	fraction_coords_list = []
	for i in atoms_list:
		element_list.append(i[0]+i[1][::-1])
		fraction_coords_list.append(i[2:])
	for i in range(len(fraction_coords_list)):
		for j in range(len(fraction_coords_list[i])):
			index = fraction_coords_list[i][j].find('(')
			if index > 0:
				fraction_coords_list[i][j] = fraction_coords_list[i][j][0:index] # delete (x)
			fraction_coords_list[i][j] = float(fraction_coords_list[i][j])

	return spacegroup, unit_cell_list, element_list, fraction_coords_list, Z_number


def modify_pdb(molten_salt_system):
	"""
		Description:
			make necessary modify in pdb file, e.g. delete bonds、angles and impropers

		Args:
			molten_salt_system: the simulated system

		Returns:
			no return
	"""

	with open(molten_salt_system+'.pdb', 'r') as f_old_pdb:
		lines = f_old_pdb.readlines()
	with open(molten_salt_system+'.pdb', 'w') as f_new_pdb:
		for line in lines:
			if line.find('CONECT') != -1: break
			f_new_pdb.write(line)
		f_new_pdb.write('END'+'\n')


def gen_model_for_cp2k_simple(molten_salt_system):
	"""
		Description:
			generate initial configuration for pure halide salt by pymatgen module, and transfer cif format to pdb format by Openbabel module

		Args:
			molten_salt_system: the simulated system

		Returns:
			no return
	"""

	global dir_molten_salt_system
	global total_atoms

	os.chdir(dir_molten_salt_system)
	ICSD_para_list = get_info_from_ICSD_file(molten_salt_system)
	spacegroup, unit_cell_list, element_list, fraction_coords_list, Z_number = get_structure_parameters(ICSD_para_list)
	print(spacegroup, unit_cell_list, element_list, fraction_coords_list, Z_number)
	# use pymatgen module to build initial configuration
	lattice = pmg.Lattice.from_parameters(a=unit_cell_list[0], b=unit_cell_list[1], c=unit_cell_list[2], alpha=unit_cell_list[3], beta=unit_cell_list[4], gamma=unit_cell_list[5])
	struct = pmg.Structure.from_spacegroup(spacegroup, lattice, element_list, fraction_coords_list)
	# determine supercell according to total_atoms
	Nx = Ny = Nz = round(pow(total_atoms/struct.num_sites, 1/3))
	# make supercell, cif format
	struct.make_supercell([Nx, Ny, Nz])
	# write structure into a cif file
	w = cif.CifWriter(struct)
	w.write_file(molten_salt_system+'.cif')
	# transfer cif format to pdb format by Openbabel module
	os.system('obabel -icif '+molten_salt_system+'.cif -opdb -O '+molten_salt_system+'.pdb')
	# make necessary modify, e.g. delete excess bonds、angles and impropers
	modify_pdb(molten_salt_system)
	# if initial configuration is nonorthogonal, then make a transformation by put the nonorthogonal system into a big orthogonal box using packmol module and get a orthogonal pdb file which is used by cp2k
	if unit_cell_list[3] != 90 or unit_cell_list[4] != 90 or unit_cell_list[5] != 90:
		print('nonorthogonal system! need transfored to orthogonal system!')
		change2orthogonal(molten_salt_system, unit_cell_list, Nx, Ny, Nz)


def change2orthogonal(molten_salt_system, unit_cell_list, Nx, Ny, Nz):
	"""
		Description:
			if initial configuration is nonorthogonal, then make a transformation by putting the nonorthogonal system into a bigger orthogonal box using packmol module and get a orthogonal pdb file which is afterward used by cp2k

		Args:
			molten_salt_system: the simulated system

		Returns:
			no return
	"""

	os.system('mv '+molten_salt_system+'.pdb tmp.pdb')
	with open(os.path.join(os.getcwd(), molten_salt_system+'.inp'), 'w') as f_inp:
		f_inp.write('tolerance 2.0'+'\n')
		f_inp.write('filetype pdb'+'\n')
		f_inp.write('output '+molten_salt_system+'.pdb'+'\n')
		f_inp.write('\n')
		f_inp.write('structure tmp.pdb'+'\n')
		f_inp.write('number 1'+'\n')
		a = round(unit_cell_list[0]*Nx*1.5, 2)
		b = round(unit_cell_list[1]*Ny*1.5, 2)
		c = round(unit_cell_list[2]*Nz*1.5, 2)
		f_inp.write('inside box 0 0 0 '+str(a)+' '+str(b)+' '+str(c)+'\n')
		f_inp.write('end structure'+'\n')
	os.system('yhrun -N 1 -n 1 -p work packmol < '+molten_salt_system+'.inp')
	while True:
		if glob.glob(molten_salt_system+'.pdb') == []:
			time.sleep(3)
			continue
		else:
			break
	# modify header and write length of simulated box into pdb file, for convience in function get_liquid_model()
	with open(molten_salt_system+'.pdb', 'r') as f_old_pdb:
		lines = f_old_pdb.readlines()
	tmp_list = []
	with open(molten_salt_system+'.pdb', 'w') as f_new_pdb:
		for line in lines:
			if line.startswith('REMARK'): continue
			tmp_list.append(line)
		tmp_str = "CRYST1"+str(a).rjust(9)+str(b).rjust(9)+str(c).rjust(9)+str(90.00).rjust(7)+str(90.00).rjust(7)+str(90.00).rjust(7)+" P 1         1\n"
		tmp_list.insert(2, tmp_str)
		f_new_pdb.writelines(tmp_list)


def write_input_file(molten_salt_system, input_pdb_file, PIM_para_list, Tstr, ensemble, steps, length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z, timestep = 1.0):
	"""
		Description:
			generate input file of cp2k named xx.in

		Args:
			molten_salt_system: the simulated system
			input_pdb_file: the pdb file is the initial configuration 
			Tstr: a string end with 'K', e.g. '1500K'
			ensemble: the simulated condition, i.e. 'NPT' or 'NVT'
			steps: the total simulated time, unit in fs

		Returns:
			no return
	"""

	T = int(Tstr.strip('K'))
	inputfilename = molten_salt_system+'-'+str(T)+'K'+'-'+ensemble+'.in'
	with open(inputfilename, 'w') as f_input:
		f_input.write('# PIMD of '+molten_salt_system+'\n')
		for i in PIM_para_list:
			if i[0] == 'Reference':
				f_input.write('# Reference '+' '.join(i[1:])+'\n')
				break
		f_input.write('\n')
		f_input.write('@SET SYSTEM '+str(T)+'K'+'-'+ensemble+'\n')
		f_input.write('\n')
		#GLOBAL
		f_input.write('&GLOBAL'+'\n')
		f_input.write('    '+'PROJECT ${SYSTEM}'+'\n')
		f_input.write('    '+'RUN_TYPE MD'+'\n')
		f_input.write('    '+'IOLEVEL MEDIUM'+'\n')
		f_input.write('    '+'&PRINT MEDIUM'+'\n')
		f_input.write('    '*2+'PHYSCON TRUE'+'\n')
		f_input.write('    '+'&END PRINT'+'\n')
		f_input.write('    '+'&PROGRAM_RUN_INFO MEDIUM'+'\n')
		f_input.write('    '+'&END PROGRAM_RUN_INFO'+'\n')
		f_input.write('    '+'&REFERENCES MEDIUM'+'\n')
		f_input.write('    '+'&END REFERENCES'+'\n')
		f_input.write('&END GLOBAL'+'\n')
		f_input.write('\n')
		#FORCE_EVAL
		f_input.write('&FORCE_EVAL'+'\n')
		f_input.write('    '+'METHOD FIST'+'\n')
		f_input.write('    '+'STRESS_TENSOR ANALYTICAL'+'\n')
		f_input.write('    '+'&MM'+'\n')
		f_input.write('    '*2+'&FORCEFIELD'+'\n')
		f_input.write('    '*3+'&SPLINE'+'\n')
		f_input.write('    '*4+'EMAX_ACCURACY [hartree] 2E-2'+'\n')
		f_input.write('    '*4+'EMAX_SPLINE [hartree] 5E1'+'\n')
		f_input.write('    '*4+'EPS_SPLINE [hartree] 1E-7'+'\n')
		f_input.write('    '*3+'&END SPLINE'+'\n')
		f_input.write('    '*3+'DO_NONBONDED TRUE'+'\n')
		#CHARGE section
		for i in PIM_para_list:
			if i[0] == 'num_of_elements': 
				num_of_elements = int(i[1])
				break
		charge_list = []
		for i in PIM_para_list:
			if len(charge_list) > num_of_elements:
				break
			elif i[0].startswith('element_'): 
				charge_list.append(i[1:])
		for n in range(len(charge_list)):
			f_input.write('    '*3+'&CHARGE'+'\n')
			f_input.write('    '*4+'ATOM'+' '+charge_list[n][0]+'\n')
			f_input.write('    '*4+'CHARGE'+' '+charge_list[n][1]+'\n')
			f_input.write('    '*3+'&END CHARGE'+'\n')
		f_input.write('    '*3+'&NONBONDED'+'\n')
		#BMH section
		Ncount_BMH = 0
		for i in PIM_para_list:
			if i[0] == 'BMH': 
				Ncount_BMH += 1
		BMH_list = []
		for i in PIM_para_list:
			if len(BMH_list) > Ncount_BMH:
				break
			elif i[0] == 'BMH': 
				BMH_list.append(i[1:])
		for n in range(len(BMH_list)):
			if BMH_list[n][0] == 'no':
				f_input.write('    '*4+'&BMHFT'+'\n')
			elif BMH_list[n][0] == 'yes':
				f_input.write('    '*4+'&BMHFTD'+'\n')
			f_input.write('    '*5+'ATOMS'+' '+BMH_list[n][1]+' '+BMH_list[n][2]+'\n')
			f_input.write('    '*5+'A [hartree]'+' '+BMH_list[n][3]+'\n')
			f_input.write('    '*5+'B [bohr^-1]'+' '+BMH_list[n][4]+'\n')
			f_input.write('    '*5+'C [bohr^6*hartree]'+' '+BMH_list[n][5]+'\n')
			f_input.write('    '*5+'D [bohr^8*hartree]'+' '+BMH_list[n][6]+'\n')
			if BMH_list[n][0] == 'yes':
				f_input.write('    '*5+'BD [bohr^-1]'+' '+BMH_list[n][7]+'\n')
				f_input.write('    '*5+'ORDER 3'+'\n')
			f_input.write('    '*5+'RCUT'+' '+str(rcut)+'\n')
			if BMH_list[n][0] == 'no':
				f_input.write('    '*4+'&END BMHFT'+'\n')
			elif BMH_list[n][0] == 'yes':
				f_input.write('    '*4+'&END BMHFTD'+'\n')
		f_input.write('    '*3+'&END NONBONDED'+'\n')
		#DIPOLE section
		Ncount_DIPOLE = 0
		for i in PIM_para_list:
			if i[0] == 'DIPOLE': 
				Ncount_DIPOLE += 1
		DIPOLE_list = []
		for i in PIM_para_list:
			if len(DIPOLE_list) > Ncount_DIPOLE:
				break
			elif i[0] == 'DIPOLE': 
				DIPOLE_list.append(i[1:])
		for n in range(len(DIPOLE_list)):
			f_input.write('    '*3+'&DIPOLE'+'\n')
			f_input.write('    '*4+'ATOM'+' '+DIPOLE_list[n][0]+'\n')
			f_input.write('    '*4+'APOL [bohr^3]'+' '+DIPOLE_list[n][2]+'\n')
			f_input.write('    '*4+'&DAMPING'+'\n')
			f_input.write('    '*5+'ATOM'+' '+DIPOLE_list[n][1]+'\n')
			f_input.write('    '*5+'BIJ [bohr^-1]'+' '+DIPOLE_list[n][3]+'\n')
			f_input.write('    '*5+'CIJ'+' '+DIPOLE_list[n][4]+'\n')
			f_input.write('    '*5+'ORDER 4'+'\n')
			f_input.write('    '*5+'TYPE TANG-TOENNIES'+'\n')
			f_input.write('    '*4+'&END DAMPING'+'\n')
			f_input.write('    '*3+'&END DIPOLE'+'\n')
		f_input.write('    '*2+'&END FORCEFIELD'+'\n')
		f_input.write('    '*2+'&POISSON'+'\n')
		f_input.write('    '*3+'PERIODIC XYZ'+'\n')
		f_input.write('    '*3+'POISSON_SOLVER PERIODIC'+'\n')
		f_input.write('    '*3+'&EWALD'+'\n')
		f_input.write('    '*4+'EWALD_ACCURACY 1E-6'+'\n')
		f_input.write('    '*4+'EWALD_TYPE EWALD'+'\n')
		f_input.write('    '*4+'GMAX'+' '+str(gmax_x)+' '+str(gmax_y)+' '+str(gmax_z)+'\n')
		f_input.write('    '*4+'&MULTIPOLES TRUE'+'\n')
		f_input.write('    '*5+'EPS_POL 5E-8'+'\n')
		f_input.write('    '*5+'MAX_IPOL_ITER 500'+'\n')
		f_input.write('    '*5+'MAX_MULTIPOLE_EXPANSION DIPOLE'+'\n')
		f_input.write('    '*5+'POL_SCF CONJUGATE_GRADIENT'+'\n') # Linear conjugate-gradient optimization of the sum of the electrostatic and induction energy. This method does not support non-linear polarization but is sometimes faster
		f_input.write('    '*4+'&END MULTIPOLES'+'\n')
		f_input.write('    '*3+'&END EWALD'+'\n')
		f_input.write('    '*2+'&END POISSON'+'\n')
		f_input.write('    '+'&END MM'+'\n')
		f_input.write('    '+'&SUBSYS'+'\n')
		f_input.write('    '*2+'&CELL'+'\n')
		f_input.write('    '*3+'ABC'+' '+str(length_x)+' '+str(length_y)+' '+str(length_z)+'\n')
		f_input.write('    '*3+'ALPHA_BETA_GAMMA 90 90 90'+'\n')
		f_input.write('    '*2+'&END CELL'+'\n')
		f_input.write('    '*2+'&TOPOLOGY'+'\n')
		f_input.write('    '*3+'CONN_FILE_FORMAT OFF'+'\n')
		f_input.write('    '*3+'COORD_FILE_FORMAT PDB'+'\n')
		f_input.write('    '*3+'COORD_FILE_NAME'+' ' +input_pdb_file+'\n')
		f_input.write('    '*2+'&END TOPOLOGY'+'\n')
		f_input.write('    '+'&END SUBSYS'+'\n')
		f_input.write('&END FORCE_EVAL'+'\n')
		f_input.write('\n')
		# MOTION
		f_input.write('&MOTION'+'\n')
		f_input.write('    '+'&PRINT'+'\n')
		f_input.write('    '*2+'&CELL MEDIUM'+'\n')
		f_input.write('    '*3+'ADD_LAST SYMBOLIC'+'\n')
		f_input.write('    '*3+'&EACH'+'\n')
		f_input.write('    '*4+'MD 10'+'\n')
		f_input.write('    '*3+'&END EACH'+'\n')
		f_input.write('    '*2+'&END CELL'+'\n')
		f_input.write('    '*2+'&FORCES MEDIUM'+'\n')
		f_input.write('    '*3+'ADD_LAST SYMBOLIC'+'\n')
		f_input.write('    '*3+'&EACH'+'\n')
		f_input.write('    '*4+'MD 100'+'\n')
		f_input.write('    '*3+'&END EACH'+'\n')
		f_input.write('    '*2+'&END FORCES'+'\n')
		f_input.write('    '*2+'&STRESS MEDIUM'+'\n')
		f_input.write('    '*3+'ADD_LAST SYMBOLIC'+'\n')
		f_input.write('    '*3+'&EACH'+'\n')
		f_input.write('    '*4+'MD 10'+'\n')
		f_input.write('    '*3+'&END EACH'+'\n')
		f_input.write('    '*2+'&END STRESS'+'\n')
		f_input.write('    '*2+'&TRAJECTORY MEDIUM'+'\n')
		f_input.write('    '*3+'LOG_PRINT_KEY TRUE'+'\n')
		f_input.write('    '*3+'FORMAT PDB'+'\n')
		f_input.write('    '*3+'FILENAME =${SYSTEM}.pdb'+'\n')
		f_input.write('    '*3+'ADD_LAST SYMBOLIC'+'\n')
		f_input.write('    '*3+'&EACH'+'\n')
		f_input.write('    '*4+'MD 100'+'\n')
		f_input.write('    '*3+'&END EACH'+'\n')
		f_input.write('    '*2+'&END TRAJECTORY'+'\n')
		f_input.write('    '*2+'&VELOCITIES MEDIUM'+'\n')
		f_input.write('    '*3+'ADD_LAST SYMBOLIC'+'\n')
		f_input.write('    '*3+'&EACH'+'\n')
		f_input.write('    '*4+'MD 10'+'\n')
		f_input.write('    '*3+'&END EACH'+'\n')
		f_input.write('    '*2+'&END VELOCITIES'+'\n')
		f_input.write('    '*2+'&RESTART_HISTORY OFF'+'\n')
		f_input.write('    '*2+'&END RESTART_HISTORY'+'\n')
		f_input.write('    '*2+'&RESTART OFF'+'\n')
		f_input.write('    '*2+'&END RESTART'+'\n')
		f_input.write('    '+'&END PRINT'+'\n')
		f_input.write('    '+'&MD'+'\n')
		f_input.write('    '*2+'&PRINT'+'\n')
		f_input.write('    '*3+'&ENERGY MEDIUM'+'\n')
		f_input.write('    '*4+'ADD_LAST SYMBOLIC'+'\n')
		f_input.write('    '*4+'&EACH'+'\n')
		f_input.write('    '*5+'MD 10'+'\n')
		f_input.write('    '*4+'&END EACH'+'\n')
		f_input.write('    '*3+'&END ENERGY'+'\n')
		f_input.write('    '*2+'&END PRINT'+'\n')
		if ensemble == 'NPT':
			f_input.write('    '*2+'ENSEMBLE NPT_I'+'\n')
		elif ensemble == 'NVT':
			f_input.write('    '*2+'ENSEMBLE NVT'+'\n')
		f_input.write('    '*2+'TIMESTEP [fs] '+str(timestep)+'\n')
		f_input.write('    '*2+'STEPS '+str(steps)+'\n')
		f_input.write('    '*2+'TEMPERATURE [K]'+' '+str(T)+'\n')
		if ensemble == 'NPT':
			f_input.write('    '*2+'&BAROSTAT'+'\n')
			f_input.write('    '*3+'PRESSURE [bar] 1'+'\n')
			f_input.write('    '*3+'TIMECON 500'+'\n')
			f_input.write('    '*2+'&END BAROSTAT'+'\n')
		elif ensemble == 'NVT':
			pass
		f_input.write('    '*2+'&THERMOSTAT'+'\n')
		f_input.write('    '*3+'REGION GLOBAL'+'\n')
		f_input.write('    '*3+'TYPE NOSE'+'\n')
		f_input.write('    '*3+'&NOSE'+'\n')
		f_input.write('    '*4+'TIMECON 500'+'\n')
		f_input.write('    '*3+'&END NOSE'+'\n')
		f_input.write('    '*2+'&END THERMOSTAT'+'\n')
		f_input.write('    '+'&END MD'+'\n')
		f_input.write('&END MOTION'+'\n')


def get_dir(path):
	"""
		Description:
			get subdirectories info, the result is stored in a list named dir_list, e.g. ['900K', '950K', '1000K', ..., '1650K', 'gen_liquid_model']

		Args:
			no arg

		Returns:
			return dir_list
	"""

	dir_list = []
	file_list = os.listdir(path)
	for file in file_list:
		newpath = os.path.join(path, file)
		if os.path.isdir(newpath):
			dir_list.append(os.path.basename(newpath))

	return dir_list


def submit(sh_name, ensemble):
	"""
		Description:
			firstly a shell file is generated and then the task is submitted by yhbatch command

		Args:
			sh_name: the name of the shell file
			ensemble: the simulated condition, i.e. 'NPT' or 'NVT' 

		Returns:
			no return
	"""

	# write a shell file, and submit the task
	file_list = os.listdir()
	for file in file_list:
		if file.endswith(ensemble+'.in'):
			input_file = file
			output_file = file.split('.')[0]+'.out'
			break
	sh_file = sh_name+'-'+ensemble+'.sh'
	with open(os.path.join(os.getcwd(), sh_file), 'w') as f_sh:
		f_sh.write('#!/bin/bash'+'\n')
		f_sh.write('yhrun -N 3 -n 72 -p work cp2k.psmp -i '+input_file+' > '+output_file+'\n') # use 3 nodes and 72 cpus
	os.system('yhbatch -N 3 '+sh_file+' '+molten_salt_system)


def resubmit(slurmID):
	"""
		Description:
			resubmit failed job according to yhcontrol and delete the slurm file

		Args:
			slurmID: slurmID of the failed job

		Returns:
			no return
	"""

	for line in os.popen('yhcontrol show job '+slurmID).readlines():
		if line.find('Command') != -1:
			tmp_list = list(filter(None, line.strip().split(' ')))
			sh_name = '-'.join(tmp_list[0].split('/')[-1].split('-')[0:-1])
			ensemble = tmp_list[0].split('/')[-1].split('-')[-1].split('.')[0]
			break
	os.system('rm '+slurm_filename)
	submit(sh_name, ensemble)


def detect_job_status(seconds):
	"""
		Description:
			# detect job status

		Args:
			seconds: the frequncy to detect job status, the value should not be too small, otherwise the sacct command will encounter error

		Returns:
			no return
	"""

	# firstly, wait until slurm file has been generated, acquire slurmID
	while True:
		if glob.glob('slurm*') == []:
			time.sleep(3)
			continue
		else:
			slurm_filename = glob.glob('slurm*')[0]
			break
	slurmID = slurm_filename.split('.')[0].split('-')[1]
	# detect task status according to slurmID, i.e. yhacct -j slurmID, it'll return the job status
	while True:
		file_size = os.path.getsize(slurm_filename)
		if file_size != 0:
			with open(slurm_filename, 'r') as f_slurm:
				slurm_str = f_slurm.read()
			if ('Node fail' in slurm_str) or ('Bus error' in slurm_str) or ('Input/output error' in slurm_str):
				pass
			else:
				print('error when running!')
				exit()
		status_str = os.popen('sacct -j '+slurmID).read()
		if 'RUNNING' in status_str:
			time.sleep(seconds)			
			continue
		if 'COMPLETED' in status_str:
			break
		if ('FAILED' in status_str) or ('NODE_FAIL' in status_str):
			# resubmit this job according to yhcontrol and delete the slurm file
			time.sleep(seconds)
			resubmit(slurmID)


def get_liquid_model(molten_salt_system, ensemble, steps, Temperature_list, PIM_para_list, timestep):
	"""
		Description:
			generate liquid_model, the liquid model is a completely melted system by heating the system to a very high temperature

		Args:
			molten_salt_system: the simulated system
			ensemble: the simulated condition, i.e. 'NPT' or 'NVT'
			steps: the total simulated time, unit in fs

		Returns:
			no return
	"""

	global dir_molten_salt_system

	dir_list = get_dir(dir_molten_salt_system)
	if 'gen_liquid_model' not in dir_list:
		print("gen_liquid_model directory doesn't exists! please check!")
		exit()
	# get into gen_liquid_model directory, submit a task, waiting for task completed, acquire liquid_model
	dir_liquid_model = os.path.join(dir_molten_salt_system, 'gen_liquid_model')
	os.chdir(dir_liquid_model)
	# copy orthogonal pdb file to this directory
	os.system('cp '+dir_molten_salt_system+'/'+molten_salt_system+'.pdb .')
	length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z = get_info_from_structure_file(molten_salt_system+'.pdb')
	write_input_file(molten_salt_system, molten_salt_system+'.pdb', PIM_para_list, Temperature_list[-2], ensemble, steps, length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z, timestep)
	submit('gen_liquid_model', ensemble)
	detect_job_status(100)


def get_info_from_structure_file(pdb_file):
	"""
		Description:
			get neccesary info to generate input of cp2k from pdb file

		Args:
			pdb_file: the structure file, e.g. 'xxx.pdb'

		Returns:
			return some value
	"""

	with open(pdb_file, 'r') as f_structure:
		for line in f_structure:
			if line.startswith('CRYST1'):
				tmp_list = list(filter(None,line.strip().split(' ')))
				length_x = float(tmp_list[1])
				length_y = float(tmp_list[2])
				length_z = float(tmp_list[3])
				break
	length_min = min(length_x, length_y, length_z)
	rcut = length_min//2
	if round(length_x) % 2 == 0:
		gmax_x = round(length_x)+1
	else:
		gmax_x = round(length_x)
	if round(length_y) % 2 == 0:
		gmax_y = round(length_y)+1
	else:
		gmax_y = round(length_y)
	if round(length_z) % 2 == 0:
		gmax_z = round(length_z)+1
	else:
		gmax_z = round(length_z)

	return length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z


def intercept_last_configuration(intercepted_pdb_file, ensemble, Tstr):
	"""
		Description:
			use a trajectory file (very large file size) as input, intercept the last configuration, and generate a new pdb file (e.g. 1500K-NVT-last-configuration.pdb), which is the input configuration of next computation

		Args:
			intercepted_pdb_file: the trajectory file that is to be intercepted
			ensemble: the simulated condition, i.e. 'NPT' or 'NVT'
			Tstr: a string end with K, e.g. '1500K'

		Returns:
			no return
	"""

	all_configuration_list = []
	matched_list = []
	last_configuration_list = []
	with open(intercepted_pdb_file, 'r') as f_old:
		pointer = -1
		for line in f_old:
			all_configuration_list.append(line)
			pointer += 1
			if line.startswith('REMARK'): matched_list.append(pointer)
	last_configuration_list = all_configuration_list[matched_list[-1]:]
	del all_configuration_list
	output_pdb_file = Tstr+'-'+ensemble+'-last-configuration.pdb'
	with open(output_pdb_file, 'w') as f_new:
		for line in last_configuration_list:
			f_new.write(line)
	del last_configuration_list


def NPT_tasks_batch_submit(molten_salt_system, ensemble, steps, Temperature_list, PIM_para_list, timestep):
	"""
		Description:
			NPT tasks batch submit jobs for each calculated temperature, the results xxx.ener is used to compute specific heat capacity while xxx.cell is used to compute density and thermal expansion coefficient

		Args:
			molten_salt_system: the simulated system
			ensemble: the simulated condition, i.e. 'NPT'
			steps: the total simulated time, unit in fs

		Returns:
			no return
	"""

	global dir_molten_salt_system

	dir_list = get_dir(dir_molten_salt_system)
	if 'gen_liquid_model' in dir_list:
		dir_list.remove('gen_liquid_model')
	dir_liquid_model = os.path.join(dir_molten_salt_system, 'gen_liquid_model')
	os.chdir(dir_liquid_model)
	file_list = os.listdir(dir_liquid_model)
	for file in file_list:
		if file.endswith('pdb'):
			if file.startswith(molten_salt_system):
				pass
			else:
				intercepted_pdb_file = file
	# intercept the pdb file to acquire the last configuration, generate a last.pdb file
	intercept_last_configuration(intercepted_pdb_file, 'NPT', Temperature_list[-2])
	os.chdir(dir_molten_salt_system)
	for dirname in dir_list:
		os.chdir(os.path.join(dir_molten_salt_system, dirname))
		os.system('cp '+dir_liquid_model+'/*last-configuration.pdb '+molten_salt_system+'.pdb')
		length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z = get_info_from_structure_file(molten_salt_system+'.pdb')
		write_input_file(molten_salt_system, molten_salt_system+'.pdb', PIM_para_list, dirname, ensemble, steps, length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z, timestep)
		submit(dirname, ensemble)
		os.chdir(dir_molten_salt_system)
	print('All '+ensemble+' tasks have been submitted!')


def detect_submitted_tasks_status(ensemble):
	"""
		Description:
			detect all submitted tasks' status, only when all submitted tasks of the system are completed, the program go on

		Args:
			ensemble: the simulated condition, i.e. 'NPT' or 'NVT'

		Returns:
			no return
	"""

	global dir_molten_salt_system

	dir_list = get_dir(dir_molten_salt_system)
	if 'gen_liquid_model' in dir_list:
		dir_list.remove('gen_liquid_model')
	for dirname in dir_list:
		os.chdir(os.path.join(dir_molten_salt_system, dirname))
		detect_job_status(100)
		os.chdir(dir_molten_salt_system)
	print('All submitted tasks have COMPLETED!')
	# rename NPT slurm-xxx.out to log.ensemble-xxx
	for dirname in dir_list:
		os.chdir(os.path.join(dir_molten_salt_system, dirname))
		file_list = os.listdir()
		for file in file_list:
			if file.startswith('slurm'):
				slurmID = file.split('.')[0].split('-')[1]
				os.system('mv slurm-'+slurmID+'.out log.'+ensemble+'-'+slurmID)
				break
		os.chdir(dir_molten_salt_system)


def NVT_tasks_batch_submit(molten_salt_system, ensemble, steps, PIM_para_list, timestep):
	"""
		Description:
			NVT tasks batch submit jobs for each calculated temperature, the results trajectory xxx.pdb is used to compute RDF etc. structure info and self-diffusion coefficient

		Args:
			molten_salt_system: the simulated system
			ensemble: the simulated condition, i.e. or 'NVT'
			steps: the total simulated time, unit in fs

		Returns:
			no return
	"""

	global dir_molten_salt_system

	dir_list = get_dir(dir_molten_salt_system)
	if 'gen_liquid_model' in dir_list:
		dir_list.remove('gen_liquid_model')
	for dirname in dir_list:
		os.chdir(os.path.join(dir_molten_salt_system, dirname))
		# intercept the pdb file to acquire the last configuration, generate a last.pdb file
		intercepted_pdb_file = dirname+'-'+'NPT'+'.pdb'
		intercept_last_configuration(intercepted_pdb_file, 'NPT', dirname)
		# get the input pdb filename
		file_list = os.listdir()
		for file in file_list:
			if file.endswith('.pdb'):
				if file.strip('.pdb').endswith('NPT-last-configuration'):
					input_pdb_file = file
					break
		length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z = get_info_from_structure_file(input_pdb_file)
		write_input_file(molten_salt_system, input_pdb_file, PIM_para_list, dirname, ensemble, steps, length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z, timestep)
		submit(dirname, ensemble)
		os.chdir(dir_molten_salt_system)
	print('All '+ensemble+' tasks have been submitted!')


def detect_viscosity_tasks_status(ensemble):
	"""
		Description:
			detect viscosity tasks' status, only when all viscosity tasks of the system are completed, the program go on

		Args:
			ensemble: the simulated condition, i.e. 'NVT'

		Returns:
			no return
	"""

	global dir_molten_salt_system

	dir_list = get_dir(dir_molten_salt_system)
	if 'gen_liquid_model' in dir_list:
		dir_list.remove('gen_liquid_model')
	for dirname in dir_list:
		os.chdir(os.path.join(dir_molten_salt_system, dirname, 'viscosity'))
		detect_job_status(100)
		os.chdir(dir_molten_salt_system)
	print('All viscosity tasks have COMPLETED!')
	# rename NPT slurm-xxx.out to log.ensemble-xxx
	for dirname in dir_list:
		os.chdir(os.path.join(dir_molten_salt_system, dirname, 'viscosity'))
		file_list = os.listdir()
		for file in file_list:
			if file.startswith('slurm'):
				slurmID = file.split('.')[0].split('-')[1]
				os.system('mv slurm-'+slurmID+'.out log.viscosity-'+ensemble+'-'+slurmID)
				break
		os.chdir(dir_molten_salt_system)


def viscosity_NVT_tasks_batch_submit(molten_salt_system, ensemble, steps, PIM_para_list, timestep):
	"""
		Description:
			viscosity NVT tasks batch submit jobs for each calculated temperature, the results xxx.stress is used to compute ACF of stress and trapezoidal integral, then compute viscosity according to Green-Kubo formula

		Args:
			molten_salt_system: the simulated system
			ensemble: the simulated condition, i.e. 'NPT' or 'NVT'
			steps: the total simulated time, unit in fs

		Returns:
			no return
	"""

	global dir_molten_salt_system

	dir_list = get_dir(dir_molten_salt_system)
	if 'gen_liquid_model' in dir_list:
		dir_list.remove('gen_liquid_model')
	for dirname in dir_list:
		os.chdir(os.path.join(dir_molten_salt_system, dirname))
		# intercept the pdb file to acquire the last configuration, generate a last.pdb file
		intercepted_pdb_file = dirname+'-'+'NVT'+'.pdb'
		intercept_last_configuration(intercepted_pdb_file, 'NVT', dirname)
		# get the input pdb filename
		file_list = os.listdir()
		for file in file_list:
			if file.endswith('.pdb'):
				if file.strip('.pdb').endswith('NVT-last-configuration'):
					input_pdb_file = file
					break
		# grenerate viscosity directory
		dir_viscosity = os.path.join(dir_molten_salt_system, dirname, 'viscosity')
		if not os.path.exists(dir_viscosity):
			os.mkdir(dir_viscosity)
		os.chdir(dir_viscosity)
		os.system('cp ../'+input_pdb_file+' .')
		length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z = get_info_from_structure_file(input_pdb_file)
		write_input_file(molten_salt_system, input_pdb_file, PIM_para_list, dirname, ensemble, steps, length_x, length_y, length_z, rcut, gmax_x, gmax_y, gmax_z, timestep)
		submit('viscosity-'+dirname, ensemble)
		os.chdir(dir_molten_salt_system)
	print('All '+ensemble+' tasks have been submitted!')
