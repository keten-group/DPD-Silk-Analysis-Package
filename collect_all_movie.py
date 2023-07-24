# import cv2
import numpy as np
import glob
import shutil
import os

sys_names = ['erate0p01'] # This would allow for looping over multiple simulations.

dele_origin = 0

for ss in range(len(sys_names)): # Outer loop is unnecessary for analyzing a single simulation.
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
