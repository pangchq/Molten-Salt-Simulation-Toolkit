#!/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Gechuanqi Pan'

# Structure of directories
#										dir_origin
#											|
#										 Results
#     										|
#     -----------------------------------------------------------------------
#     |					|					|				|				|
# system_1			system_2			system_3			……			system_n
#	  |                                                                     |
#   -----------------                                              -------------------------------
#   |    |          |                                              |    |   |   |   |            |
#  T_1  T_2	gen_liquid_model									   T_1  T_2 T_3 T_4 ……   gen_liquid_model

import os
import sys
import time
import glob
import argparse
import PIM_unitary_utils as utils

parser = argparse.ArgumentParser(description = 'A high-throughput MD workflow for thermophysical and microstructural properties of Polarizable Ion Model molten haldie salts.')
parser.add_argument('name', type = str, metavar = 'name', help = 'The name of the molten haldie system.')
parser.add_argument('-t', '--timestep', type = float, nargs = 1, help = "Set timestep of MD. If not set, the default value is 1 fs.")
parser.add_argument('-gmfc', '--gen_model_for_cp2k', action = 'store_false', default = True, help = "If set, do not generate model for cp2k.")
parser.add_argument('-glm', '--get_liquid_model', action = 'store_false', default = True, help = "If set, do not generate liquid model.")
parser.add_argument('-NPT', '--NPT_tasks', action = 'store_false', default = True, help = "If set, do not carry out NPT tasks.")
parser.add_argument('-NVT', '--NVT_tasks', action = 'store_false', default = True, help = "If set, do not carry out NVT tasks.")
parser.add_argument('-vNVT', '--viscosity_NVT_tasks', action = 'store_false', default = True, help = "If set, do not carry out viscosity NVT tasks")
args = parser.parse_args()

molten_salt_system = args.name
utils.molten_salt_system = molten_salt_system

total_atoms = 1000
utils.total_atoms = total_atoms
build_model_type = 'ICSD'
utils.build_model_type = build_model_type

# dir_origin is the path of workspace
dir_origin = os.getcwd()
utils.origin = dir_origin

# create Results directory
dir_Results = os.path.join(dir_origin, 'Results')
utils.dir_Results = dir_Results
if not os.path.exists(dir_Results):
	os.mkdir(dir_Results)
os.chdir(dir_Results)

# create directory for molten_salt_system
dir_molten_salt_system = os.path.join(dir_Results, molten_salt_system)
utils.dir_molten_salt_system = dir_molten_salt_system
if not os.path.exists(dir_molten_salt_system):
	os.mkdir(dir_molten_salt_system)
os.chdir(dir_molten_salt_system)

# get info from PIM_paramters_file, stored in a list para_list_ICSD
PIM_para_list = utils.get_info_from_PIM_file(molten_salt_system)

# create subdirectories by temperature for molten_salt_system
Temperature_list = utils.gen_directory(molten_salt_system, PIM_para_list)

if args.gen_model_for_cp2k:
	# generate initial configuration by pymatgen module, and transfer cif format to pdb format by Openbabel module
	utils.gen_model_for_cp2k_simple(molten_salt_system)

# determine timestep of MD
if args.timestep:
    timestep = args.timestep[0]
else:
    timestep = 1.0

if args.get_liquid_model:
	# generate liquid_model, NPT ensemble, need total steps
	utils.get_liquid_model(molten_salt_system, 'NPT', 100000, Temperature_list, PIM_para_list, timestep)

if args.NPT_tasks:
	# NPT tasks batch submit jobs for each calculated temperature, need total steps
	# the results x.ener is used to compute specific heat capacity while x.cell is used to compute density and thermal expansion coefficient
	utils.NPT_tasks_batch_submit(molten_salt_system, 'NPT', 500000, Temperature_list, PIM_para_list, timestep) # NPT ensemble lasts for 500 ps to ensure that the box was in equilibrium  
	utils.detect_submitted_tasks_status('NPT')

if args.NVT_tasks:
	# batch submit jobs for each calculated temperature, NVT ensemble, need total steps
	# the results trajectory x.pdb is used to compute RDF etc. structure info and self-diffusion coefficient
	utils.NVT_tasks_batch_submit(molten_salt_system, 'NVT', 500000, PIM_para_list, timestep) # Finally, 500 ps production run was collected for different configurations
	utils.detect_submitted_tasks_status('NVT')

if args.viscosity_NVT_tasks:
	# viscosity NVT tasks batch submit jobs for each calculated temperature, need total steps
	# the results x.stress is used to compute ACF of stress and trap, then compute viscosity according to Green-Kubo formula
	utils.viscosity_NVT_tasks_batch_submit(molten_salt_system, 'NVT', 3000000, PIM_para_list, timestep)
	utils.detect_viscosity_tasks_status('NVT')
