import subprocess

def runScwrl4(mutant):
    subprocess.run(f"Scwrl4 -i {mutant} -o {mutant}", shell = True, stdout = subprocess.DEVNULL) 
    return mutant 
