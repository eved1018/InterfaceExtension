import subprocess

cmd = "python intercaat.py -pdb complex.pdb -fp /home/evan/Desktop/ee_research/Einstein_Project/MutationModeller/Data/1I8L/Output/Models/TYR_104/TYR_104_C/ -qc C,104 -ic A -mi 0" 
p = subprocess.run(cmd, shell = True,check = True, stdout=subprocess.PIPE, universal_newlines=True)
lines = p.stdout.rstrip().split()[5]
print(lines)