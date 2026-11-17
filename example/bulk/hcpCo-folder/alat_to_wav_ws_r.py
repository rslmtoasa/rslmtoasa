import numpy as np

print("Only implemented for hcp structure with primitive unit cell\n")

#a is part of the base, and c is the height
a = float(input("Write value for a:\n"))
c = 4.0211440800000000
#c is calculated from (c/a)*a
a1 = [1/2, np.sqrt(3)/2, 0]
a2 = [1/2, -np.sqrt(3)/2, 0]
a3 = [0, 0, c/a]
#Volume for the unit cell (primitive cell) for hcp structure in Å
V = np.sqrt(3)*a**2*c/2
#Wigner-Seitz radius for each atom (assumed equal volume, 2 atoms in unit cell)
wav = (3*V/(8*np.pi))**(1/3)
# ws_r in atomic units, hence we divide by the Bohr radius a0 = 0.529177 Å
ws_r = wav/0.529177
print("")
print(f"wav = {wav} (Å) \n")
print(f"ws_r = {ws_r} (Bohr radius)")
print("")
print("The lattice vectors in terms of a is:")
print(f"a1 = {a1}")
print(f"a2 = {a2}")
print(f"a3 = {a3}")
