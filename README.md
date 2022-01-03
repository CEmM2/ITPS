# ITPS
Micro structural optimzation of composite RVE
for the scripts to run one have to install both easyPBC and micromechanics plug ins first 
files explanation: gmsh_final reads data from geom_datafile and makes the meshing 
mech_analasys reads they mesh .it is an abaqus scripting file for running a single simulation
mech_therm analasys runs the script from python 
run_all is a batch file running gmsh final and mech_term file in a loop over the geometry data
. all files have to be in the same directory to work
