import numpy as np
import glob
import shutil
import os

dele_origin = 0

ii = 0
if not (os.path.exists(f'movie')):
    os.mkdir(f'movie')
for filename in glob.glob(f'equil_evolve_*/movie.png'):
    shutil.copyfile(filename, f'movie/movie_{ii}_equil.png')
    if dele_origin:
        os.remove(filename)
    ii += 1

for filename in glob.glob(f'shear_evolve_*/movie.png'):
    shutil.copyfile(filename, f'movie/movie_{ii}_shear.png')
    if dele_origin:
        os.remove(filename)
    ii += 1

for filename in glob.glob(f'stretch_evolve_*/movie.png'):
    shutil.copyfile(filename, f'movie/movie_{ii}_stretch.png')
    if dele_origin:
        os.remove(filename)
    ii += 1
