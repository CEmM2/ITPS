FOR /L %%I IN (1,1,1600) DO (
python gmsh_final.py
python mech_therm_analasys.py)
