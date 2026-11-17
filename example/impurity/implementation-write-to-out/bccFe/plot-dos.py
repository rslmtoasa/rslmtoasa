import matplotlib.pyplot as plt
import numpy as np
import sys

def get_orb_column(orb):
   if orb == "s":
      column = 1
   elif orb == "px":
      column = 2
   elif orb == "py":
      column = 3
   elif orb == "pz":
      column = 4
   else:
      print("Orbital provided not available. Exit program")
      sys.exit()
   return column

#Finds Fermi energy
def find_fermi_energy(file_path):
   fermi_energy = None

   # Open the file and read line by line
   with open(file_path, "r") as file:
      for line in file:
         # Look for the line containing "Fermi energy:"
         if "Fermi energy:" in line:
               # Split the line to extract the numerical value
               fermi_energy = float(line.split(":")[1].strip())
               break  # Stop searching once found

   # Print the extracted Fermi energy
   if fermi_energy is not None:
      print(f"Fermi energy: {fermi_energy}")
   else:
      print("Fermi energy not found in the file.")
      sys.exit()
   return fermi_energy

fermi_energy = find_fermi_energy("report.out")
paths = []
labels = []
columns = []
filename = input("Set name of pdf:\n")
title = input("Set title:\n")
xlabel = r"$E - E_F$"
ylabel = "DOS"
list_of_orbitals = []


more = True
while more:
   print("Provide filename for data")
   path = input("Enter x in fort.x:\n")
   if path.startswith("45"):
      print("fort.45x was provided.")
      atype = input("Enter atom type (element/surface/etc...) for label:\n")
      print("Select which orbitals to include. Enter 0 when done.")
      print("Available orbitals are s, px, py ,pz, d, ... (fix this)")
      cont = True
      while cont:
         orbital = input("Enter orbital (or 0 if done):\n")
         if orbital == "0":
            cont = False
         else:
            label = atype + "_" + orbital
            paths.append("fort." + path)
            labels.append(label)
            columns.append(get_orb_column(orbital))
   else:
      label = input("Enter label:\n")
      paths.append("fort." + path)
      labels.append(label)
      columns.append(1)
   logical = True
   while logical:
      check = input("Add another plot? (y/n):\n")
      if check == "y":
         logical = False
      elif check == "n":
         more = False
         logical = False
      else:
         pass

plt.figure(figsize=(8, 6))
plt.suptitle(title)
for i, path in enumerate(paths):
   data = np.genfromtxt(path)
   en = data[:,0]
   en = (en - fermi_energy)*13.605693
   dos = data[:,columns[i]]
   plt.plot(en, dos, label=labels[i])
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()
plt.tight_layout()
plt.savefig(filename+".pdf")
   