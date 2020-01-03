# Molten-Salt-Simulation-Toolkit
A high-throughput molecular dynamics workflow to compute thermophysical properties and microstructure information of molten halide salts.
# Note
Anyone would like to rerun the code should note that the HT-MD workflow is based on slurm file system of "Tianhe" No. 2 (National supercomputer center in Guangzhou, China), so you maybe need to do some modifications to adapt to other systems. Morover, the code contains some third-party softwares, i.e. pymatgen, openbabal, packmol, CP2K, MDAnalysis, scipy, numpy, matplotlib, bson etc. These softwares should be installed in advance. Also, a mongodb is required.
# Usage
1. modify the path where parameter files are stored in PIM_unitary_utils.py.
2. run the PIM_unitary.py, you can get useful info by type "python PIM_unitary.py -h".
3. run the data_post_process.py, also you can get useful info by type "data_post_process.py -h".
