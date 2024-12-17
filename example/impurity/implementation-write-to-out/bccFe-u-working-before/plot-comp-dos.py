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

# fermi_energy = find_fermi_energy("report.out")
fermi_energies = []
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
   path = input("Provide filename for data (it should end with .125, .251, etc.):\n")
   path_to_report = input("Input path to the report.out file (to extract fermi energy):\n")
   if path[-4:-1] == ".45":
      print("sting.45x was provided.")
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
            paths.append(path)
            fermi_energies.append(find_fermi_energy(path_to_report))
            labels.append(label)
            columns.append(get_orb_column(orbital))
   else:
      label = input("Enter label:\n")
      paths.append(path)
      fermi_energies.append(find_fermi_energy(path_to_report))
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
   en = (en - fermi_energies[i])*13.605693
   dos = data[:,columns[i]]
   plt.plot(en, dos, label=labels[i])
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()
plt.tight_layout()
plt.savefig(filename+".pdf")
   