#!/bin/env python
# -*- coding: utf-8 -*-
# use numpy to accerlate computation

__author__ = 'Gechuanqi Pan'

import os
import time
import glob
import re
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import spline
import matplotlib.pyplot as plt
import MDAnalysis as MDA
from MDAnalysis.analysis.rdf import InterRDF
from bson.binary import Binary
import pickle
import pprint

global Periodic_table_of_elements, Mass_of_elements
Periodic_table_of_elements = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm','Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo')
Mass_of_elements = (1.008, 4.003, 6.941, 9.012, 10.811, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990, 24.305, 26.982, 28.086, 30.974, 32.065, 35.453, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.409, 69.723, 72.640, 74.922, 78.960, 79.904, 83.798, 85.468, 87.620, 88.906, 91.224, 92.906, 95.940, 97.907, 101.070, 102.906, 106.420, 107.868, 112.411, 114.818, 118.710, 121.760, 127.600, 126.904, 131.293, 132.905, 137.327, 138.905, 140.116, 140.908, 144.242, 145.000, 150.360, 151.964, 157.250, 158.925, 162.500, 164.930, 167.259, 168.934, 173.040, 174.967, 178.490, 180.948, 183.840, 186.207, 190.230, 192.217, 195.084, 196.967, 200.590, 204.383, 207.200, 208.980, 208.982, 209.987, 222.018, 223, 226.000, 227.000, 232.038, 231.036, 238.029, 237.000, 244.000, 243.000, 247.000, 247.000, 251.000, 252.000, 257.000, 258.000, 259.000, 262.000, 261.000, 262.000, 266.000, 264.000, 277.000, 268.000, 281.000, 272.000, 285.000, 284.000, 289.000, 288.000, 292.000, 291.000, 293.000)


#-------------------------------------------dividing line------------------------------------------------#
# the following functions are used to batch data post process

def get_dir_temperature():
	"""
		Description:
			get all temperature subdirectories named 'xxxK' and stored in a global list named dir_list

		Args:
			No args

		Returns:
			No return
	"""

	global dir_molten_salt_system
	global dir_list

	dir_list = []
	file_list = os.listdir(dir_molten_salt_system)
	for file in file_list:
		newpath = os.path.join(dir_molten_salt_system, file)
		if os.path.isdir(newpath):
			dir_list.append(os.path.basename(newpath))
	dir_list.remove('gen_liquid_model')


def dir_sort_bubble():
	"""
		Description:
			sort the dir_list so that the temperature from low to high

		Args:
			No args

		Returns:
			No return
	"""

	global dir_molten_salt_system
	global dir_list

	tmp_list = dir_list
	tmp_list = [int(x.strip('K')) for x in tmp_list]
	for i in range(len(tmp_list)):
		for j in range(i+1, len(tmp_list)):
			if tmp_list[i] > tmp_list[j]:
				tmp_list[i], tmp_list[j] = tmp_list[j], tmp_list[i]
	dir_list = [str(x)+'K' for x in tmp_list]


def read_data_from_file(filename):
	"""
		Description:
			read data from a file

		Args:
			filename: the name of the file to be read

		Returns:
			a data_list that store all needed data, 2 dimension
	"""

	data_list = []
	with open(filename, 'r') as f_read:
		for line in f_read:
			if line.find('#') >= 0: continue
			#注意文件分隔符！
			tmp_list = list(filter(None, line.strip().split(' ')))
			data_list.append(tmp_list)
	
	return data_list


def read_data_from_trajectory(filename):
	"""
		Description:
			read data from a trajectory, it is different from other files because of the pdb file format

		Args:
			filename: the name of the file to be read

		Returns:
			a data_list that store all needed data, 2 dimension
	"""

	data_list = []
	with open(filename, 'r') as f_read:
		for line in f_read:
			if (line.find('TITLE') >= 0) or (line.find('AUTHOR') >= 0) or (line.find('REMARK') >= 0) or (line.find('CRYST1') >= 0) or (line.find('END') >= 0): continue
			tmp_list = [line[0:4].strip(' '), line[4:11].strip(' '), line[11:14].strip(' '), line[30:38].strip(' '), line[38:46].strip(' '), line[46:54].strip(' '), line[54:60].strip(' '), line[60:66].strip(' '), line[66:].strip('\n').strip(' ')]
			data_list.append(tmp_list)
	
	return data_list


def changeunit(computed_property, software, T = 300, V = 10000, Nevery = 10, timestep = 1):
	"""
		Description:
			unit conversion according to computed property and used md software

		Args:
			computed_property: computed property, i.e 'density' or 'viscosity' or 'thermal conductivity' or 
							  'thermal expansion coefficient' or 'constant pressure specific heat capacity' or 
							  'self-diffusion coefficient'
			software: used md software, i.e. 'lammps' or 'cp2k'
			T: temperature, needed by compute ACF
			V: volume, needed by compute ACF
			Nevery: use input values every this many timesteps, needed by compute ACF
			timestep: timestep of md

		Returns:
			a scale number
	"""

	if computed_property == 'density':
		scale = 1/6.023e23/pow(1e-8, 3)
	elif computed_property == 'viscosity':
		dt = timestep
		if software == 'lammps':
			kB = 1.3806504e-23 # J/K Boltzman constant
			atm2Pa = 101325
			A2m = 1e-10
			fs2s = 1e-15
			convert = atm2Pa*atm2Pa*fs2s*A2m*A2m*A2m
			scale = convert/(kB*T)*V*Nevery*dt
		elif software == 'cp2k':
			kB = 1.3806504e-23 # J/K Boltzman constant
			bar2Pa = 1e5
			A2m = 1e-10
			fs2s = 1e-15
			convert = bar2Pa*bar2Pa*fs2s*A2m*A2m*A2m
			scale = convert/(kB*T)*V*Nevery*dt
		else:
			print('unkown software!')
			exit()
	elif computed_property == 'thermal conductivity':
		pass
	elif computed_property == 'thermal expansion coefficient':
		pass
	elif computed_property == 'constant pressure specific heat capacity':
		scale = 6.27509468713739e5
	elif computed_property == 'self-diffusion coefficient':
		scale = pow(1e-8, 2)/1e-15
	else:
		print('unsupported property unit conversion!')
		exit()
	
	return scale


def get_force_field_doi_from_PIM_file(molten_salt_system):
	"""
		Description:
			get force field doi from PIM_paramters_file

		Args:
			molten_salt_system: the simulated system

		Returns:
			return a string
	"""

	parameters_filename = '/WORK/nscc-gz_material_1/pangchq/parameter_files/PIM_parameters_file'
	with open(parameters_filename, 'r')  as f_para_PIM:
		# firstly, judge if system is in PIM_parameters_file(and unique), if false, report an error and exit
		flag_exist = False
		for num,line in enumerate(f_para_PIM):
			if line.strip('\n') == ('#'+molten_salt_system):
				flag_exist = True
				break
		if flag_exist == False:
			print(molten_salt_system+" is not in PIM_parameters_file, error!")
			exit()
		f_para_PIM.seek(0, 0)

		flag = False
		for num,line in enumerate(f_para_PIM):
			if line.strip('\n') == ('#'+molten_salt_system):
				flag = True
				continue
			if line.strip('\n') == ('#END '+molten_salt_system):
				break
			if flag and line.startswith('doi'):
				force_field_doi = list(filter(None, line.strip().split(' ')))[1]
				break

	return force_field_doi


def compute_density(molten_salt_system, cell_output_freq, last_steps):
	"""
		Description:
			compute density according to the last xxps of NPT ensemble, and write results into file

		Args:
			molten_salt_system: the simulated system
			cell_output_freq: the output frequency of xx.cell file
			last_steps: the last xx steps that are used to compute the equilibrium volume

		Returns:
			return a dict
	"""

	global dir_molten_salt_system
	global dir_list
	global density_list

	# determine the number of each element of this simulated system according to the molten_salt_system.cif file
	with open(os.path.join(dir_molten_salt_system, molten_salt_system+'.cif'), 'r') as f_cif:
		for line in f_cif:
			if line.startswith('_chemical_formula_sum'):
				chemical_formula_list = list(filter(None, re.split('[ \']', line.strip('\n'))))[1:]
				break
	# use re to split letter and digit, and stored in a list, the element in this list is a truple, e.g. ('Li', 288)
	atoms_list = []
	for i in range(len(chemical_formula_list)):
		element = re.findall(r"\D+", chemical_formula_list[i])[0]
		number = int(re.findall(r"\d+", chemical_formula_list[i])[0])
		atoms_list.append((element, number))

	# get total mass of the simulated system
	mass_sum = 0
	for k in range(len(atoms_list)):
		this_element, number_of_this_element = atoms_list[k][0], atoms_list[k][1]
		if this_element in Periodic_table_of_elements:
			index = Periodic_table_of_elements.index(this_element)
			mass_of_this_element = Mass_of_elements[index]
			mass_sum += mass_of_this_element*number_of_this_element

	with open(os.path.join(dir_molten_salt_system, 'property.density'), 'w') as f_property:
		f_property.write('-'*10+'density'+'-'*10+'\n')
		f_property.write('Temperature/K'+'\t'+'rho/g*cm^-3'+'\n')
		T_list, density_list, density_dict= [], [], {}
		# note unit conversion, the unit of density is g/cm^3
		scale = changeunit('density', 'cp2k')
		for dirname in dir_list:
			os.chdir(os.path.join(dir_molten_salt_system, dirname))
			T_list.append(int(dirname.strip('K')))
			filename = dirname+'-NPT-1.cell'
			data_list = read_data_from_file(filename)
			exacted_data_list = data_list[-(last_steps//cell_output_freq):]
			volume_average = np.sum(np.array([[float(x) for x in exacted_data_list[y]] for y in range(len(exacted_data_list))]), axis = 0)[11]/len(exacted_data_list)
			density = mass_sum/volume_average*scale
			density = round(density, 3)
			density_dict[dirname] = density
			density_list.append(density)
			f_property.write('\t'.join([dirname, str(density)])+'\n')
		# determine linear fitting formula of density versus T and write into results file
		coefficient_list = np.polyfit(np.array(T_list), np.array(density_list), 1)
		intercept, first_order_coefficient = format(coefficient_list[-1], '0.3f'), format(coefficient_list[-2], '0.2e')
		fitting_formula = 'rho = '+intercept+' + '+first_order_coefficient+'*T(K)'
		f_property.write('fitting of density: '+fitting_formula+'\n')
		density_dict['fiting formula'] = fitting_formula

	return density_dict


def compute_thermal_expansion_coefficient(molten_salt_system, temperature_interval):
	"""
		Description:
			compute thermal expansion coefficient, need density results, and write results into file

		Args:
			molten_salt_system: the simulated system
			temperature_interval: the interval of each temperature, e.g. 50

		Returns:
			return a dict
	"""

	global dir_molten_salt_system
	global dir_list
	global density_list

	with open(os.path.join(dir_molten_salt_system, 'property.thermal_expansion_coefficient'), 'w') as f_property:
		f_property.write('-'*10+'thermal_expansion_coefficient'+'-'*10+'\n')
		f_property.write('Temperature/K'+'\t'+'beta/K^-1'+'\n')
		T_list, thermal_expansion_coefficient_list, thermal_expansion_coefficient_dict = [], [], {}
		for i in range(1, len(dir_list)-1):
			T_list.append(int(dir_list[i].strip('K')))
			rho_current_T = density_list[i]
			rho_previous_T = density_list[i-1]
			rho_next_T = density_list[i+1]
			thermal_expansion_coefficient = ((1/rho_current_T)*(rho_previous_T-rho_next_T)/(2*temperature_interval))
			thermal_expansion_coefficient = round(thermal_expansion_coefficient, 10)
			thermal_expansion_coefficient_dict[dir_list[i]] = thermal_expansion_coefficient
			thermal_expansion_coefficient_list.append(thermal_expansion_coefficient)
			f_property.write('\t'.join([dir_list[i], str(thermal_expansion_coefficient)])+'\n')
		# determine linear fitting formula of thermal_expansion_coefficient versus T and write into results file
		coefficient_list = np.polyfit(np.array(T_list), np.array(thermal_expansion_coefficient_list), 1)
		intercept, first_order_coefficient = format(coefficient_list[-1], '0.2e'), format(coefficient_list[-2], '0.3e')
		fitting_formula = 'beta = '+intercept+' + '+first_order_coefficient+'*T(K)'
		f_property.write('fitting of thermal_expansion_coefficient: '+fitting_formula+'\n')
		thermal_expansion_coefficient_dict['fitting formula'] = fitting_formula

	return thermal_expansion_coefficient_dict
			

def compute_constant_pressure_specific_heat_capacity(molten_salt_system, ener_output_freq, last_steps):
	"""
		Description:
			compute constant pressure specific heat capacity, and write results into file

		Args:
			molten_salt_system: the simulated system
			ener_output_freq: the output frequency of xx.ener file
			last_steps: the last xx steps that are used to compute the equilibrium energy

		Returns:
			return a constane value
	"""

	global dir_molten_salt_system
	global dir_list

	# determine the number of molecules of this simulated system according to the molten_salt_system.cif file
	with open(os.path.join(dir_molten_salt_system, molten_salt_system+'.cif'), 'r') as f_cif:
		for line in f_cif:
			if line.startswith('_cell_formula_units_Z'):
				number_of_molecule = int(list(filter(None, line.strip('\n').split(' ')))[1])
				break

	with open(os.path.join(dir_molten_salt_system, 'property.constant_pressure_specific_heat_capacity'), 'w') as f_property:
		f_property.write('-'*10+'constant_pressure_specific_heat_capacity'+'-'*10+'\n')
		f_property.write('Cp/cal*K^-1*mol^-1'+'\n')
		T_list, ener_list = [], []
		# note unit conversion, the unit of constant_pressure_specific_heat_capacity is cal*K^-1*mol^-1
		scale = changeunit('constant pressure specific heat capacity', 'cp2k')
		for dirname in dir_list:
			os.chdir(os.path.join(dir_molten_salt_system, dirname))
			T_list.append(int(dirname.strip('K')))
			filename = dirname+'-NPT-1.ener'
			data_list = read_data_from_file(filename)
			exacted_data_list = data_list[-(last_steps//ener_output_freq):]
			exacted_data_list = [[float(x) for x in exacted_data_list[y]] for y in range(len(exacted_data_list))] # change string to float
			# it's worth noting that Cons Qty[a.u.] doesn't equal to Kin.[a.u.] plus Pot.[a.u.], here Kin.[a.u.] plus Pot.[a.u.] is exactly the total energy
			kin_ener_average = np.sum(np.array(exacted_data_list), axis = 0)[2]/len(exacted_data_list)
			pot_ener_average = np.sum(np.array(exacted_data_list), axis = 0)[4]/len(exacted_data_list)
			total_ener_average = kin_ener_average+pot_ener_average
			ener_list.append(total_ener_average)
		coefficient_list = np.polyfit(np.array(T_list), np.array(ener_list), 1)
		intercept, first_order_coefficient = coefficient_list[-1], coefficient_list[-2]
		constant_pressure_specific_heat_capacity = round(first_order_coefficient*scale/number_of_molecule, 3)
		f_property.write(str(constant_pressure_specific_heat_capacity)+'\n')

	return constant_pressure_specific_heat_capacity


def MSD_computer(atoms_list, atom_sum, frames_MSD, exacted_data_list, frame_interval = 1, at_per_mol = 1):
	"""
		Description:
			compute MSD, just support one element temporarily, group will be added later

		Args:
			atoms_list: contain the element and number of this element, e.g. [('Li', 288), ('Be', 144), ('F', 576)]
			atom_sum: total # of all atoms
			frames_MSD: the # of frames that are used to compute MSD
			exacted_data_list: exacted data according frames_MSD
			frame_interval: interval of frame to get reference pos
			at_per_mol: # of atoms per molecule

		Returns:
			return a dict
	"""

	trajectory_dict, MSD_dict, histogram_dict = {}, {}, {}

	# exact coordinates of each element of each frame, seperate coordinates by element
	# trajectory_dict[this_element][i] is ith frame
	# trajectory_dict[this_element][i][j] is jth atom
	# trajectory_dict[this_element][i][j][k] is kth dimension coordinate
	for i in range(len(atoms_list)):
		this_element, number_of_this_element = atoms_list[i][0], atoms_list[i][1]
		trajectory_dict[this_element] = []
		for frame in range(frames_MSD):
			tmp_list = []
			for k in range(atom_sum):
				if exacted_data_list[(atom_sum)*frame+k][0] == 'ATOM':
					if exacted_data_list[(atom_sum)*frame+k][8] == this_element:
						tmp_list.append([float(x) for x in exacted_data_list[(atom_sum)*frame+k][3:6]]) # just exact coordinates of each atom, note that data type must be transfered to float
			trajectory_dict[this_element].append(tmp_list)

	# transfer list to array to accerlate computation
	for key in trajectory_dict:
		trajectory_dict[key] = np.array(trajectory_dict[key])

	# begin to compute MSD
	for k in range(len(atoms_list)):
		this_element, number_of_this_element = atoms_list[k][0], atoms_list[k][1]
		histogram_dict[this_element] = [[] for i in range(frames_MSD)] # histogram store the msd of each frame
		starting_point_pos = []
		n_mols = len(trajectory_dict[this_element][0])//at_per_mol
		for frame, pos in enumerate(trajectory_dict[this_element]):
			# for molecule, compute position of center of mass; for atom, one atom is one molecule
			pos_mols = np.zeros((n_mols, 3))
			for i in range(n_mols):
				pos_atoms = pos[at_per_mol*i:at_per_mol*(i+1)] # get position array for all the atoms in this molecule, 2 dimension array
				if this_element in Periodic_table_of_elements:
					index = Periodic_table_of_elements.index(key)
					mass_this_element = Mass_of_elements[index]
				at_masses = np.array([mass_this_element]) # get mass array for all the atoms in the molecule, 1 dimension, 1 dimension array

				for i_at in range(at_per_mol):
					pos_mols[i, :] += at_masses[i_at]*pos_atoms[i_at, :]
				pos_mols[i, :] /= np.sum(at_masses) # compute pos of center of mass of this molecule, which is used to compute MSD of molecule, 1 dimension array

			if frame % frame_interval == 0:
				starting_point_pos.append(np.copy(pos_mols)) # get reference pos at fixed frame_interval

			for i_sp in range(len(starting_point_pos)-1, -1, -1):
				index_delta = frame - i_sp*frame_interval
				if index_delta >= frames_MSD:
					break
				# (starting_point_pos[i_sp]-pos_mols) is dx,dy,dz
				# (starting_point_pos[i_sp]-pos_mols)**2 is dx^2,dy^2,dz^2
				# np.sum((starting_point_pos[i_sp]-pos_mols)**2,axis=1) is dr^2=dx^2+dy^2+dz^2
				# np.mean compute average
				histogram_dict[this_element][index_delta].append(np.mean(np.sum((pos_mols-starting_point_pos[i_sp])**2, axis=1)))
	del trajectory_dict

	# average over the histograms of each element to adquire MSD, note that len(histogram_dict[this_element][i]) is not a constant
	for k in range(len(atoms_list)):
		this_element, number_of_this_element = atoms_list[k][0], atoms_list[k][1]
		MSD_dict[this_element] = np.zeros(frames_MSD)
		for i in range(frames_MSD):
			MSD_dict[this_element][i] = np.mean(np.array(histogram_dict[this_element][i]))
	del histogram_dict

	return MSD_dict


def self_diffusion_function(x, a, b):
	return a*np.exp(b/x)


def compute_self_diffusion_coefficient(molten_salt_system, trajectory_output_freq, last_steps, timestep = 1):
	"""
		Description:
			compute self-diffusion coefficient, and write results into file

		Args:
			molten_salt_system: the simulated system
			trajectory_output_freq: the output frequency of trajectory file
			last_steps: the last xx steps used to compute the MSD

		Returns:
			return a dict
	"""

	global dir_molten_salt_system
	global dir_list

	# determine the number of each element of this simulated system according to the molten_salt_system.cif file, then count the number of all atoms
	with open(os.path.join(dir_molten_salt_system, molten_salt_system+'.cif'), 'r') as f_cif:
		for line in f_cif:
			if line.startswith('_chemical_formula_sum'):
				chemical_formula_list = list(filter(None, re.split('[ \']', line.strip('\n'))))[1:]
				break

	# use re to split letter and digit, and stored in a dict
	atoms_list = []
	for i in range(len(chemical_formula_list)):
		element = re.findall(r"\D+", chemical_formula_list[i])[0]
		number = int(re.findall(r"\d+", chemical_formula_list[i])[0])
		atoms_list.append((element, number))

	# get the total # of all atoms
	atom_sum = 0
	for i in range(len(atoms_list)):
		if atoms_list[i][0] in Periodic_table_of_elements:
			atom_sum += atoms_list[i][1]

	with open(os.path.join(dir_molten_salt_system, 'property.self_diffusion_coefficient'), 'w') as f_property:
		f_property.write('-'*10+'self_diffusion_coefficient'+'-'*10+'\n')
		f_property.write('Temperature/K'+'\telement D/cm^2*s^-1'*len(atoms_list)+'\n')
		T_list, self_diffusion_coefficient_list, self_diffusion_coefficient_dict = [], [[] for i in range(len(atoms_list))], [{} for i in range(len(atoms_list))]
		# note unit conversion, the unit of self_diffusion_coefficient is cm^2*s^-1
		scale = changeunit('self-diffusion coefficient', 'cp2k')
		for dirname in dir_list:
			os.chdir(os.path.join(dir_molten_salt_system, dirname))
			T_list.append(int(dirname.strip('K')))
			# determine the # of frames used to compute MSD
			filename = dirname+'-NVT.pdb'
			data_list = read_data_from_trajectory(filename)
			frames_MSD = last_steps//trajectory_output_freq+1 # one frame more which is the start frame
			exacted_data_list = data_list[-atom_sum*frames_MSD:]
			# get the MSD versus frame of each element
			MSD_dict = MSD_computer(atoms_list, atom_sum, frames_MSD, exacted_data_list, 1, 1)
			# fitting MSD versus time and divided by 6 to get self_diffusion_coefficient of each temperature, just use the middle data
			time_list = [i*trajectory_output_freq*timestep for i in range(frames_MSD)] # unit fs
			f_property.write(dirname)
			for k in range(len(atoms_list)):
				this_element, number_of_this_element = atoms_list[k][0], atoms_list[k][1]
				coefficient_list = np.polyfit(np.array(time_list)[:], MSD_dict[this_element][:], 1)
				intercept, first_order_coefficient = coefficient_list[-1], coefficient_list[-2]
				self_diffusion_coefficient = first_order_coefficient/6*scale
				self_diffusion_coefficient = round(self_diffusion_coefficient, 10)
				self_diffusion_coefficient_list[k].append(self_diffusion_coefficient)
				f_property.write('\t'+this_element+' '+format(self_diffusion_coefficient, '0.3e'))
			f_property.write('\n')
			del MSD_dict
			os.chdir(dir_molten_salt_system)
		# write results into self_diffusion_coefficient_dict according to self_diffusion_coefficient_list
		for k in range(len(atoms_list)):
			this_element, number_of_this_element = atoms_list[k][0], atoms_list[k][1]
			self_diffusion_coefficient_dict[k]['element'] = this_element
			index = 0
			for dirname in dir_list:
				self_diffusion_coefficient_dict[k][dirname] = self_diffusion_coefficient_list[k][index]
				index += 1
		# determine non-linear fitting formula of self_diffusion_coefficient versus T and write into results file and self_diffusion_coefficient_dict
		x = np.array(T_list)
		for k in range(len(atoms_list)):
			this_element, number_of_this_element = atoms_list[k][0], atoms_list[k][1]
			y = np.array(self_diffusion_coefficient_list[k])
			popt, pcov = curve_fit(self_diffusion_function, x, y)
			a, b = format(popt[0], '0.3e'), format(popt[1], '0.0f')
			fitting_formula = 'D = '+a+' x exp('+b+'/T) (K)'
			f_property.write('fitting of self-diffusion coefficient of '+this_element+': '+fitting_formula+'\n')
			self_diffusion_coefficient_dict[k]['fitting formula'] = fitting_formula

	return self_diffusion_coefficient_dict


def ACF_computer_method_1(stress_data_list, Nevery, Nrepeat, Nfreq, Npair = 3):
	"""
		Description:
			compute ACF function, and the results are stored in ACF_data_list

		Args:
			stress_data_list: the list contain needed component of stress
			Nevery, Nrepeat, Nfreq: see https://lammps.sandia.gov/doc/fix_ave_correlate.html for explanation
			Npair: # of considered comonent of stress, e.g. 3 for xy, xz and yz, 5 for xy, xz, yz, xx-yy and 2zz-xx-yy

		Returns:
			return k_repeat, ACF_data_list
	"""

	# 开始计算自相关并存储在ACF_datalist中,格式模仿lammps的fix ave/correlate command,需要根据Nevery、Nrepeat、Nfreq,含义见https://lammps.sandia.gov/doc/fix_ave_correlate.html
	Index = 0
	TimeDelta = 0
	Ncount = 0
	k_repeat = int((len(stress_data_list) - 1)*Nevery/Nfreq)
	ACF_data_list= []

	# Timestep为0时单独写Cij,前2行单独写,3到every_lentgh行用for循环,格式与lammps的fix_ave_correlate命令计算结果文件完全一样
	# 第一行
	ACF_data_list.append([0, Nrepeat])
	# 第二行
	tmplist = []
	for n in range(Npair):
		tmplist.append(stress_data_list[0][n+1]*stress_data_list[0][n+1])
	ACF_data_list.append([1, 0, 1] + tmplist)
	# 第3到every_lentgh行
	for i in range(2, Nrepeat + 1):
		ACF_data_list.append([i, (i-1)*10, 0]+ [0.0]*Npair)

	# Timestep为非0时,求Cij
	# 重复k次,每次先添加Timestep Number-of-time-windows
	for k in range(k_repeat):
		ACF_data_list.append([(k+1)*Nfreq, Nrepeat])
		# 每个i对应一行,重复Nrepeat次
		for i in range(Nrepeat):
			Index = i+1
			TimeDelta = i*10
			# 每一行的相关长度Ncount
			Ncount = (Nfreq//Nevery)*(k+1)+1-i
			# 为这一行求自相关, 需要先求和, 再求平均
			tmp_list = []
			for n in range(Npair):
				tmp_sum = 0
				tmp_average_sum = 0
				# 通过是设置Nstart和Nend,使得累加不从头开始,而是从上一次的和开始
				Nstart = (Nfreq//Nevery)*k+1-i
				Nend = Ncount
				# 这里要先判断一下,以Nfreq//Nevery为长度划分,当k = 0时会遇到Nstart < 0的情况,此时需要设Nstart = 0,即k = 0时range(Nstart, Nend)这个区间的长度不足Nfreq//Nevery的长度
				if Nstart < 0: Nstart = 0
				for m in range(Nstart, Nend):
					tmp_sum += stress_data_list[m][n+1]*stress_data_list[m+i][n+1]
				tmp_sum  += ACF_data_list[(Nrepeat+1)*k+i+1][n+3]*Nstart
				tmp_average_sum = tmp_sum/Ncount
				tmp_list.append(tmp_average_sum)
			ACF_data_list.append([Index, TimeDelta, Ncount]+tmp_list)

	return k_repeat, ACF_data_list


def TRAP_computer_method_1(ACF_data_list, k_repeat, Nrepeat, Npair, scale):
	"""
		Description:
			trapezoidal integral,梯形积分

		Args:
			ACF_data_list, k_repeat: return of ACF_computer(stress_data_list, Nevery, Nrepeat, Nfreq, Npair)
			Nrepeat: see https://lammps.sandia.gov/doc/fix_ave_correlate.html for explanation
			Npair: # of considered comonent of stress, e.g. 3 for xy, xz and yz, 5 for xy, xz, yz, xx-yy and 2zz-xx-yy
			scale: unit conversion scale number

		Returns:
			return trap_data_list
	"""

	trap_data_list = []
	for k in range(k_repeat+1):
		trap_data_list.append(ACF_data_list[k*(Nrepeat+1)])
		tmp_sum = [0 for x in range(Npair)]
		for i in range(Nrepeat):
			tmp_list = [i+1]
			for n in range(Npair):
				if i == 0:
					tmp_sum[n] += 0.5*ACF_data_list[k*(Nrepeat+1)+1+i][n+3]
				elif i == Nrepeat-1:
					tmp_sum[n] += 0.5*ACF_data_list[k*(Nrepeat+1)+1+i][n+3]
				else:
					tmp_sum[n] += ACF_data_list[k*(Nrepeat+1)+1+i][n+3]
				tmp_list.append(round(tmp_sum[n]*scale, 6))
			trap_data_list.append(tmp_list)

	return trap_data_list


def viscosity_function(x, a, b):
	return a*np.exp(b/x)


def compute_viscosity_method_1(molten_salt_system, Nevery, Nrepeat, Nfreq, Npair, timstep = 1):
	"""
		Description:
			compute viscosity, firstly compute ACF of stress of xy, xz and yz by function of ACF_computer(), then do integral by TRAP_computer(), at last determine a average value

		Args:
			molten_salt_system: the simulated system
			Nevery, Nrepeat, Nfreq: see https://lammps.sandia.gov/doc/fix_ave_correlate.html for explanation
			Npair: # of considered comonent of stress, e.g. 3 for xy, xz and yz, 5 for xy, xz, yz, xx-yy and 2zz-xx-yy
			timstep: timestep of md

		Returns:
			return a dict
	"""

	global dir_molten_salt_system
	global dir_list

	with open(os.path.join(dir_molten_salt_system, 'property.viscosity'), 'w') as f_property:
		f_property.write('-'*10+'viscosity'+'-'*10+'\n')
		f_property.write('Temperature/K'+'\t'+'viscosity/mPa*s'+'\n')
		T_list, viscosity_list, viscosity_dict = [], [], {}
		for dirname in dir_list:
			os.chdir(os.path.join(os.path.join(dir_molten_salt_system, dirname), 'viscosity'))
			T_list.append(int(dirname.strip('K')))
			filename = dirname+'-NVT-1.stress'
			data_list = read_data_from_file(filename)
			# get needed stress component data, and stored in stress_data_list
			stress_data_list = []
			for tmp_list in data_list:
				tmp_list = [int(x) if x in tmp_list[0:1] else float(x) for x in tmp_list]
				step = tmp_list[0]
				stress_xy = tmp_list[3]
				stress_xz = tmp_list[4]
				stress_yz = tmp_list[7]
				stress_xx_yy = tmp_list[2]-tmp_list[6]
				stress_2zz_xx_yy = 2*tmp_list[10]-tmp_list[2]-tmp_list[6]
				# [xy, xz, yz, xx-yy, 2zz-xx-yy]
				tmp_list = [step, stress_xy, stress_xz, stress_yz, stress_xx_yy, stress_2zz_xx_yy]
				stress_data_list.append(tmp_list)
			del data_list
			# compute ACF
			k_repeat, ACF_data_list = ACF_computer_method_1(stress_data_list, Nevery, Nrepeat, Nfreq, Npair)
			# compute integral of ACF, temperature and equilibrium volume are needed and exacted from files
			T = int(dirname.strip('K'))
			filename = dirname+'-NVT-1.cell'
			data_list = read_data_from_file(filename)
			V = float(data_list[-1][-1])
			scale = changeunit('viscosity', 'cp2k', T, V, Nevery, timstep)
			trap_data_list = TRAP_computer_method_1(ACF_data_list, k_repeat, Nrepeat, Npair, scale)
			# write ACF_data_list and trap_data_list into files
			with open('stress_ACF.dat', 'w') as f_ACF:
				for  tmp_list in ACF_data_list:
					f_ACF.write(' '.join([str(x) for x in tmp_list])+'\n')
			with open('stress_trap.dat', 'w') as f_trap:
				for  tmp_list in trap_data_list:
					f_trap.write(' '.join([str(x) for x in tmp_list])+'\n')
			average_viscosity = (trap_data_list[-1][1]+trap_data_list[-1][2]+trap_data_list[-1][3])/Npair
			average_viscosity = round(average_viscosity*1000, 3)
			viscosity_dict[dirname] = average_viscosity
			viscosity_list.append(average_viscosity)
			f_property.write('\t'.join([dirname, str(average_viscosity)])+'\n')
			os.chdir(dir_molten_salt_system)
		# determine fitting formula of viscosity
		x = np.array(T_list)
		y = np.array(viscosity_list)
		popt, pcov = curve_fit(viscosity_function, x, y)
		a, b = format(popt[0], '0.3e'), format(popt[1], '0.0f')
		fitting_formula = 'eta = '+a+' x exp('+b+'/T) (K)'
		f_property.write('fitting of viscosity: '+fitting_formula+'\n')
		viscosity_dict['fitting formula'] = fitting_formula

	return viscosity_dict


def ACF_computer_method_2(stress_data_list, Nevery, Nrepeat, Nfreq, Npair = 3, frame_interval = 5):
	"""
		Description:
			compute ACF function, and the results are stored in ACF_data_list. 与方法1相比,这个方法更直观(循环层数少),但是非常占内存,frame_interval不能太小,否则会超内存!frame_interval默认值为5,因此autocorr_histogram[index_delta]样本数变少了,从而导致k_repeat减少,所以最后积分计算出来的粘度与method 1有差别,虽然差别不大,但是需慎重!

		Args:
			stress_data_list: the list contain needed component of stress
			Nevery, Nrepeat, Nfreq: see https://lammps.sandia.gov/doc/fix_ave_correlate.html for explanation
			Npair: # of considered comonent of stress, e.g. 3 for xy, xz and yz, 5 for xy, xz, yz, xx-yy and 2zz-xx-yy
			frame_interval: interval of frame to get reference stress

		Returns:
			return k_repeat, ACF_data_list
	"""

	# transfer list to array to accerlate computation
	stress_data_array = np.array(stress_data_list)
	del stress_data_list
	#k_repeat = int((len(stress_data_array) - 1)*Nevery/Nfreq)

	# begin to compute ACF
	frames_ACF = len(stress_data_array)
	autocorr_histogram = [[] for i in range(frames_ACF)]
	starting_point_stress = []
	for frame, stress in enumerate(stress_data_array):
		if frame % frame_interval == 0:
			starting_point_stress.append(np.copy(stress)) # get reference stress at fixed frame_interval

		for i_sp in range(len(starting_point_stress)-1, -1, -1):
			index_delta = frame - i_sp*frame_interval
			if index_delta >= Nrepeat:
				break
			autocorr_histogram[index_delta].append((stress*starting_point_stress[i_sp])[1:])
	del stress_data_array
	del starting_point_stress

	# write ACF_data_list
	k_repeat = len(autocorr_histogram[0])*Nevery//Nfreq + 1
	ACF_data_list = []
	for k in range(k_repeat+1):
		ACF_data_list.append([k*Nfreq, Nrepeat])
		for i in range(Nrepeat):
			Index = i+1
			TimeDelta = i*10
			# 每一行的相关长度Ncount
			Ncount = (Nfreq//Nevery)*k+1-i
			# k = 0时需要特殊处理
			if k == 0:
				if i == 0:
					ACF_data_list.append([Index, TimeDelta, 1]+autocorr_histogram[i][0].tolist())
				else:
					ACF_data_list.append([Index, TimeDelta, 0]+[0.0]*5)
			else:
				Nstart = (Nfreq//Nevery)*(k-1)+1-i
				Nend = Ncount
				if Nstart < 0: Nstart = 0
				tmp_list = np.sum(autocorr_histogram[i][Nstart:Nend], axis=0)
				tmp_list = (tmp_list+np.array(ACF_data_list[(Nrepeat+1)*(k-1)+i+1][3:])*Nstart)/Ncount
				tmp_list = [Index, TimeDelta, Ncount]+tmp_list.tolist()
				ACF_data_list.append(tmp_list)
	del autocorr_histogram

	return k_repeat, ACF_data_list


def TRAP_computer_method_2(ACF_data_list, k_repeat, Nrepeat, Npair, scale):
	"""
		Description:
			trapezoidal integral,梯形积分,矩阵运算

		Args:
			ACF_data_list, k_repeat: return of ACF_computer(stress_data_list, Nevery, Nrepeat, Nfreq, Npair)
			Nrepeat: see https://lammps.sandia.gov/doc/fix_ave_correlate.html for explanation
			Npair: # of considered comonent of stress, e.g. 3 for xy, xz and yz, 5 for xy, xz, yz, xx-yy and 2zz-xx-yy
			scale: unit conversion scale number

		Returns:
			return trap_data_list
	"""

	trap_data_list = []
	for k in range(k_repeat+1):
		trap_data_list.append(ACF_data_list[k*(Nrepeat+1)])
		for i in range(Nrepeat):
			if i == 0:
				tmp_sum_array = np.array(ACF_data_list[k*(Nrepeat+1)+1+i])*0.5
			elif i == Nrepeat-1:
				tmp_sum_array += np.array(ACF_data_list[k*(Nrepeat+1)+1+i])*0.5
			else:
				tmp_sum_array += np.array(ACF_data_list[k*(Nrepeat+1)+1+i])
			trap_data_list.append([i+1]+[round(x, 6) for x in (tmp_sum_array*scale).tolist()[3:]])

	return trap_data_list


def compute_viscosity_method_2(molten_salt_system, Nevery, Nrepeat, Nfreq, Npair, timstep = 1):
	"""
		Description:
			compute viscosity, firstly compute ACF of stress of xy, xz and yz by function of ACF_computer(), then do integral by TRAP_computer(), at last determine a average value

		Args:
			molten_salt_system: the simulated system
			Nevery, Nrepeat, Nfreq: see https://lammps.sandia.gov/doc/fix_ave_correlate.html for explanation
			Npair: # of considered comonent of stress, e.g. 3 for xy, xz and yz, 5 for xy, xz, yz, xx-yy and 2zz-xx-yy
			timstep: timestep of md

		Returns:
			No return
	"""

	global dir_molten_salt_system
	global dir_list

	with open(os.path.join(dir_molten_salt_system, 'property.results'), 'a') as f_property:
		f_property.write('-'*10+'viscosity'+'-'*10+'\n')
		f_property.write('Temperature/K'+'\t'+'viscosity/mPa*s'+'\n')
		T_list, viscosity_list = [], []
		for dirname in dir_list:
			os.chdir(os.path.join(os.path.join(dir_molten_salt_system, dirname), 'viscosity'))
			T_list.append(int(dirname.strip('K')))
			filename = dirname+'-NVT-1.stress'
			data_list = read_data_from_file(filename)
			# get needed stress component data, and stored in stress_data_list
			stress_data_list = []
			for tmp_list in data_list:
				tmp_list = [int(x) if x in tmp_list[0:1] else float(x) for x in tmp_list]
				step = tmp_list[0]
				stress_xy = tmp_list[3]
				stress_xz = tmp_list[4]
				stress_yz = tmp_list[7]
				stress_xx_yy = tmp_list[2]-tmp_list[6]
				stress_2zz_xx_yy = 2*tmp_list[10]-tmp_list[2]-tmp_list[6]
				# [xy, xz, yz, xx-yy, 2zz-xx-yy]
				tmp_list = [step, stress_xy, stress_xz, stress_yz, stress_xx_yy, stress_2zz_xx_yy]
				stress_data_list.append(tmp_list)
			del data_list
			# compute ACF
			k_repeat, ACF_data_list = ACF_computer_method_2(stress_data_list, Nevery, Nrepeat, Nfreq, Npair)
			# compute integral of ACF, temperature and equilibrium volume are needed and exacted from files
			T = int(dirname.strip('K'))
			filename = dirname+'-NVT-1.cell'
			data_list = read_data_from_file(filename)
			V = float(data_list[-1][-1])
			scale = changeunit('viscosity', 'cp2k', T, V, Nevery, timstep)
			trap_data_list = TRAP_computer_method_2(ACF_data_list, k_repeat, Nrepeat, Npair, scale)
			# write ACF_data_list and trap_data_list into files
			with open('stress_ACF.dat', 'w') as f_ACF:
				for  tmp_list in ACF_data_list:
					f_ACF.write(' '.join([str(x) for x in tmp_list])+'\n')
			with open('stress_trap.dat', 'w') as f_trap:
				for  tmp_list in trap_data_list:
					f_trap.write(' '.join([str(x) for x in tmp_list])+'\n')
			average_viscosity = (trap_data_list[-1][1]+trap_data_list[-1][2]+trap_data_list[-1][3])/Npair
			average_viscosity = round(average_viscosity, 6)
			viscosity_list.append(average_viscosity)
			os.chdir(dir_molten_salt_system)
		# determine fitting formula of viscosity
		x = np.array(T_list)
		y = np.array(viscosity_list)
		popt, pcov = curve_fit(viscosity_function, x, y)
		a, b = format(popt[0], '0.3e'), format(popt[1], '0.0f')
		fitting_formula = 'eta = '+a+' + exp('+b+'/T) (K)'
		f_property.write('fitting of viscosity: '+fitting_formula+'\n')


def deal_with_PDB(dirname, ensemble):
	traj_origin_filename = dirname+'-'+ensemble+'.pdb'
	with open(traj_origin_filename, 'r') as f:
		data_list = []
		model_count = 1
		for line in f:
			if line.startswith('REMARK'):
				data_list.append('MODEL        '+str(model_count)+'\n')
				model_count += 1
			if line.startswith('END'):
				line = 'ENDMDL\n'
			data_list.append(line)
	with open('traj.pdb', 'w') as f:
		for line in data_list:
			f.writelines(line)


def compute_microstructure(molten_salt_system, ensemble, start, stop, step, nbins=500, rdf_flag=True, cn_flag=True, adf_flag=True):
	"""
		Description:
			compute microstructure info including radial distribution fuction(rdf), coordination number(cn) and angle distribution function(adf) using MDAnalysis module. adf is not completed

		Args:
			molten_salt_system: 
			ensemble: 
			start: 
			stop: 
			step: 
			rdf_flag: 
			cn_flag: 
			adf_flag: 

		Returns:
			return dicts according to setting
	"""

	global dir_molten_salt_system
	global dir_list

	# analysis elements and get cations and inions according cif file
	os.chdir(dir_molten_salt_system)
	cations_list, anions_list = [], []
	with open(molten_salt_system+'.cif', 'r') as f_cif:
		flag = False
		atoms_list = []
		for line in f_cif:
			if line.startswith(' _atom_type_oxidation_number'):
				flag = True
				continue
			if line.startswith('loop_'):
				flag = False
				continue
			if flag:
				tmp_tuple = tuple(list(filter(None, line.strip().split(' '))))
				atoms_list.append(tmp_tuple)
	for i in range(len(atoms_list)):
		if float(atoms_list[i][1]) > 0:
			cations_list.append((re.findall(r"[A-Za-z]+", atoms_list[i][0])[0]).upper()) # transfered to uppercase because MDAnalysis can just identify capitalization
		elif float(atoms_list[i][1]) < 0:
			anions_list.append((re.findall(r"[A-Za-z]+", atoms_list[i][0])[0]).upper())
	elements_list = cations_list+anions_list
	elements_count = len(elements_list)

	# determine # of pairs(used in rdf and cn)
	pairs_count = 0
	for i in range(elements_count):
		for j in range(i, elements_count):
			pairs_count += 1

	# radial_distribution_function_dense_dict, coordination_number_dense_dict= [{} for i in range(pairs_count)], [{} for i in range(pairs_count)] # 加密的
	radial_distribution_function_dict, coordination_number_dict= [{} for i in range(pairs_count)], [{} for i in range(pairs_count)]	# 原始的
	for dirname in dir_list:
		os.chdir(os.path.join(dir_molten_salt_system, dirname))
		deal_with_PDB(dirname, ensemble) # write MODEL and ENDMDL into xxx.pdb file and stored as traj.pdb
		universe = MDA.Universe('traj.pdb')
		# determine cutoff of rdf
		tmp_list = list(filter(None, os.popen('grep CRY '+dirname+'-'+ensemble+'.pdb').readlines()[-1].strip().split(' ')))
		equilibrium_lattice_constant_a = float(tmp_list[1])
		equilibrium_lattice_constant_b = float(tmp_list[2])
		equilibrium_lattice_constant_c = float(tmp_list[3])
		equilibrium_lattice_constant_min = min(equilibrium_lattice_constant_a, equilibrium_lattice_constant_b, equilibrium_lattice_constant_c)
		cutoff = equilibrium_lattice_constant_min//2

		# begin compute rdf and cn
		columns_count = -1
		for i in range(elements_count):
			for j in range(i, elements_count):
				columns_count += 1
				atomgroup1 = universe.select_atoms('name '+elements_list[i])
				atomgroup2 = universe.select_atoms('name '+elements_list[j])
				rdf = InterRDF(atomgroup1, atomgroup2, nbins=nbins, range=(0.0, cutoff))
				rdf.run(start=start, stop=stop, step=step) # after run, rdf and count have been generated
				rdf.rdf[0] = 0 # 第一个数莫名奇妙的超级大，本该为0
				
				# This part is discarded because of the limitation of 16MB for single file in mongodb
				'''
				# 加密10倍的
				# store rdf results into a dict
				# make curve smoother
				rdf_new = np.linspace(rdf.bins.min(), rdf.bins.max(), nbins*10) # 加密十倍, 画出来的图更光滑
				rdf_smooth = spline(rdf.bins, rdf.rdf, rdf_new)
				# get rdf_rmax(first peak) and rdf_rmin(first peak valley) of rdf
				rdf_rmax_index = np.argmax(rdf_smooth)
				rdf_rmax = rdf_new[rdf_rmax_index]
				rdf_rmin_index = np.argmin(rdf_smooth[rdf_rmax_index:])
				rdf_rmin = rdf_new[rdf_rmax_index:][rdf_rmin_index]
				radial_distribution_function_dense_dict[columns_count]['pair'] = elements_list[i]+'-'+elements_list[j]
				radial_distribution_function_dense_dict[columns_count][dirname] = {'rdf_rmax': rdf_rmax, 'rdf_rmin': rdf_rmin, 'rdf_value': Binary(pickle.dumps(np.vstack((rdf_new, rdf_smooth)), protocol=2))} # numpy必须转换成二进制才能存进pymongo
				# store cn results into a dict
				rdf_count = rdf.count.copy()
				rdf_count = rdf_count/(len(atomgroup1)*rdf.__dict__['n_frames']) # average
				cn = rdf_count.cumsum() # integrate
				# make curve smoother
				cn_smooth = spline(rdf.bins, cn, rdf_new)
				# get cn_rmin according to first peak valley in rdf
				cn_rmin = cn_smooth[rdf_rmax_index:][rdf_rmin_index]
				coordination_number_dense_dict[columns_count]['pair'] = elements_list[i]+'-'+elements_list[j]
				coordination_number_dense_dict[columns_count][dirname] = {'cn_rmin': cn_rmin, 'cn_value': Binary(pickle.dumps(np.vstack((rdf_new, cn_smooth)), protocol=2))}
				'''

				# 原始的
				# rdf
				radial_distribution_function_dict[columns_count]['pair'] = elements_list[i]+'-'+elements_list[j]
				radial_distribution_function_dict[columns_count][dirname] = Binary(pickle.dumps(np.vstack((rdf.bins, rdf.rdf)), protocol=2))
				# cn
				rdf_count = rdf.count.copy()
				rdf_count = rdf_count/(len(atomgroup1)*rdf.__dict__['n_frames']) # average
				cn = rdf_count.cumsum() # integrate
				coordination_number_dict[columns_count]['pair'] = elements_list[i]+'-'+elements_list[j]
				coordination_number_dict[columns_count][dirname] = Binary(pickle.dumps(np.vstack((rdf.bins, cn)), protocol=2))
		os.system('rm traj.pdb')
		os.chdir(dir_molten_salt_system)

	return radial_distribution_function_dict, coordination_number_dict
