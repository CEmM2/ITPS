
import subprocess
from subprocess import run

run("abaqus cae noGUI=mech_sim.py", shell=True)
print("finished calculation")