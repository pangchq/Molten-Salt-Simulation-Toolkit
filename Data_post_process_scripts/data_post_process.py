#!/bin/env python 
# -*- coding: utf-8 -*-

__author__ = 'Gechuanqi Pan'

import os
import sys
import time
import pymongo
import pprint
import argparse
import data_post_process_utils as utils

parser = argparse.ArgumentParser(description = 'Post process program of PIM molten salts to get thermal physical properties.')
parser.add_argument('name', type = str, metavar = 'name', help = 'The name of the molten salt system.')
parser.add_argument('-t', '--timestep', type = float, nargs = 1, help = "Set timestep of MD. If not set, the default value is 1 fs.")
parser.add_argument('-dst', '--density', action = 'store_false', default = True, help = "If set, do not compute density.")
parser.add_argument('-tec', '--thermal_expansion_coefficient', action = 'store_false', default = True, help = "If set, do not compute thermal_expansion_coefficient.")
parser.add_argument('-shc', '--constant_pressure_specific_heat_capacity', action = 'store_false', default = True, help = "If set, do not compute constant_pressure_specific_heat_capacity.")
parser.add_argument('-sdc', '--self_diffusion_coefficient', action = 'store_false', default = True, help = "If set, do not compute self_diffusion_coefficient.")
parser.add_argument('-vst', '--viscosity', action = 'store_false', default = True, help = "If set, do not compute viscosity.")
parser.add_argument('-mis', '--microstructure', action = 'store_false', default = True, help = "If set, do not compute microstructure.")
args = parser.parse_args()


molten_salt_system = args.name

# dir_origin is the path of workspace
dir_origin = os.getcwd()
utils.origin = dir_origin

dir_Results = os.path.join(dir_origin, 'Results')
utils.dir_Results = dir_Results
if not os.path.exists(dir_Results):
	print("dir_Results doesn't exits! error!")
	exit()

dir_molten_salt_system = os.path.join(dir_Results, molten_salt_system)
utils.dir_molten_salt_system = dir_molten_salt_system
if not os.path.exists(dir_molten_salt_system):
	print("dir_molten_salt_system doesn't exits! error!")
	exit()

start_time = time.time()

# mongodb connection
mongodb_addr = '12.11.70.140:20202'
db_name = 'molten_salts'
collection_name = 'molten_salts'
my_client = pymongo.MongoClient(mongodb_addr) # connect mongodb
my_db = my_client[db_name] # connect database
my_collection = my_db[collection_name] # connect collection
material_name = {'name': molten_salt_system}
if my_collection.find(material_name).count() == 0:
	my_collection.insert_one(material_name)

# get all temperature subdirectories
utils.get_dir_temperature()

# sort the subdirectories
utils.dir_sort_bubble()

# determine timestep of MD
if args.timestep:
	timestep = args.timestep[0]
else:
	timestep = 1.0

# add ID if 'msid' key is not in the mongodb, msid is short for molten salt ID
if 'msid' not in my_collection.find_one(material_name):
	ID_count = 1 # start from 1
	results = my_collection.find()
	for result in results:
		if 'msid' in result: ID_count += 1
	msid = 'ms_'+str(ID_count).zfill(4)
	my_collection.update_one(material_name, {'$set': {'msid': msid}})

# add force field doi if 'force_field_doi' key is not in the mongodb
if 'force_field_doi' not in my_collection.find_one(material_name):
	force_field_doi = utils.get_force_field_doi_from_PIM_file(molten_salt_system)
	my_collection.update_one(material_name, {'$set': {'force_field_doi': force_field_doi}})

if args.density:
	# compute density according to the last 200ps of NPT ensemble
	# def compute_density(molten_salt_system, cell_output_freq, last_steps)
	density_dict = utils.compute_density(molten_salt_system, 10, 200000)
	my_collection.update_one(material_name, {'$set': {'density': density_dict}})

if args.thermal_expansion_coefficient:
	# compute thermal expansion coefficient
	# def compute_thermal_expansion_coefficient(molten_salt_system, temperature_interval)
	thermal_expansion_coefficient_dict = utils.compute_thermal_expansion_coefficient(molten_salt_system, 50)
	my_collection.update_one(material_name, {'$set': {'thermal_expansion_coefficient': thermal_expansion_coefficient_dict}})

if args.constant_pressure_specific_heat_capacity:
	# compute constant pressure specific heat capacity
	# def compute_constant_pressure_specific_heat_capacity(molten_salt_system, ener_output_freq, last_steps)
	constant_pressure_specific_heat_capacity_value = utils.compute_constant_pressure_specific_heat_capacity(molten_salt_system, 10, 200000)
	my_collection.update_one(material_name, {'$set': {'constant_pressure_specific_heat_capacity': constant_pressure_specific_heat_capacity_value}})

if args.self_diffusion_coefficient:
	# compute self-diffusion coefficient
	# def compute_self_diffusion_coefficient(molten_salt_system, trajectory_output_freq, last_steps, timestep = 1)
	self_diffusion_coefficient_dict = utils.compute_self_diffusion_coefficient(molten_salt_system, 100, 200000, timestep)
	my_collection.update_one(material_name, {'$set': {'self_diffusion_coefficient': self_diffusion_coefficient_dict}})

if args.viscosity:
	# compute viscosity, firstly compute ACF of stress of xy, xz and yz by function of ACF_computer(), then do integral by TRAP_computer(), at last determine a average value
	# def compute_viscosity(molten_salt_system, Nevery, Nrepeat, Nfreq, Npair, timstep = 1)
	viscosity_dict = utils.compute_viscosity_method_1(molten_salt_system, 10, 1000, 10000, 3, timestep)
	my_collection.update_one(material_name, {'$set': {'viscosity': viscosity_dict}})

if args.microstructure:
	# compute microstructure info including radial distribution fuction(rdf), coordination number curve(cn) and angle distribution function(adf) using MDAnalysis module
	# compute_microstructure(molten_salt_system, ensemble, start, stop, step, nbins, rdf_flag=True, cn_flag=True, adf_flag=True):
	radial_distribution_function_dict, coordination_number_dict = utils.compute_microstructure(molten_salt_system, 'NPT', 3001, 5000, 1, 500)
	my_collection.update_one(material_name, {'$set': {'radial_distribution_function': radial_distribution_function_dict}})
	my_collection.update_one(material_name, {'$set': {'coordination_number': coordination_number_dict}})
#pprint.pprint(my_collection.find_one(material_name))

end_time = time.time()
print(end_time-start_time)