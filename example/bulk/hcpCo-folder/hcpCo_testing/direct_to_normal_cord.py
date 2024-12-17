import numpy as np

#Convert from direct to cartesian coord:
def dir_to_cart(dir_coord):
   a1 = np.array([0.5000000000000000, -0.8660254037844386, 0.0000000000000000], dtype=np.float64)
   a2 = np.array([0.5000000000000000, 0.8660254037844386, 0.0000000000000000], dtype=np.float64)
   a3 = np.array([0.0000000000000000, 0.0000000000000000, 1.626676471463544], dtype=np.float64)
   return dir_coord[0]*a1 + dir_coord[1]*a2 + dir_coord[2]*a3

atom1_dir = np.array([0.3333333333333333, 0.6666666666666666, 0.2500000000000000], dtype=np.float64)
atom2_dir = np.array([0.6666666666666667, 0.3333333333333334, 0.7500000000000000], dtype=np.float64)

atom1_cart = dir_to_cart(atom1_dir)
atom2_cart = dir_to_cart(atom2_dir)

print(atom1_cart)
print(atom2_cart)