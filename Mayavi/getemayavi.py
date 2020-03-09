import numpy as np
from scipy import special
from ase.io import read
import quippy
import matplotlib.pyplot as plt
from mayavi import mlab
import ase.neighborlist as nl


#function to read the xyz file
#this returns 500 xyz coordinates in a list of lists with 3 variables in each
#There are 101 sets of 500 atoms. This code reads 1 of those sets, so will need
#to set a loop

L = 7.93701
l = 6
neighbour_dist = 4.25
neighbour_dist_sqrd = neighbour_dist**2
min_count = 7
threshold_dist = 0.5
   
#creates reduced vector                    
  
def GeTe_coords(filename):    

    my_atoms = read(filename)
    
    coords = my_atoms.get_positions(wrap=False)
    
    symbols = my_atoms.get_chemical_symbols()
    
    coordinates_Ge = []
    
    coordinates_Te = []
    
    for i, s in enumerate(symbols):
        if s == 'Ge':
            coordinates_Ge.append(coords[i])
        else:
            coordinates_Te.append(coords[i])
    
    length = len(my_atoms)
      
    return coordinates_Ge, coordinates_Te

def read_xyz(filename):

    i = 0
    xyz = open(filename)
  
    while i < 2:
    
        Ge, Te = GeTe_coords(filename)
        
        x_Ge = []
        y_Ge = []
        z_Ge = []
        
        for u,v,w in Ge:
            x_Ge.append(u)
            y_Ge.append(v)
            z_Ge.append(w)
        
        x_Te = []
        y_Te = []
        z_Te = []
  
        for f,g,h in Te:
            x_Te.append(f)
            y_Te.append(g)
            z_Te.append(h)

        
        mlab.points3d(x_Ge, y_Ge, z_Ge, scale_factor=2.9,color=(0.1,0.9,1.0), mode='sphere',resolution=12,opacity=1)
        mlab.points3d(x_Te, y_Te, z_Te, scale_factor=2.9,color=(1.0,0.9,0.5), mode='sphere',resolution=12,opacity=1)
        print("Processed a configuration")
        mlab.show()
        plt.show()
        i += 1
        
    xyz.close()        
    
    return 

t= read_xyz("quenched_0.xyz")


