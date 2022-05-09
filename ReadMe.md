### Interface Extension:

![ctla4](Media/ctla4cd80.png)
------

 
### Overview:

The Interface Extension program was developed to identify pharmacologically relevant residues at the periphery of a protein’s interface from a known three-dimensional structure. First, the set of normal interface and extended interface residues are defined. The normal interface residues are then excluded from the extended interface. Lastly, the potential residues are filtered by their ability to reach the normal interface threshold upon in-silica mutation. The extended interface positions present themselves as important positions in a variety of applications, as mutations at these positions have potential to affect the protein binding without damaging the interface.

![alg](Media/alg.png)
------


### Installation and Dependencies:
* Requires python 3.6 or later. Tested on python 3.9. 
* Requires Anocanda to install Modeller easily. (https://www.anaconda.com/)
1. Install:
```sh
   pip install pericaat
```
2. Download Modeller: https://salilab.org/modeller/download_installation.html
    * for Conda enviroment:
```sh 
        conda config --add channels salilab
        conda install modeller
```
3. Download SCRWL4 (optional): http://dunbrack.fccc.edu/lab/scwrl
4. Download qhull (optional): http://www.qhull.org/

------

### Usage:
```sh
 pericaat -pdb {PDB ID} -qc {Query chain} -c {Interacting chain}
```
------

### Command Line Arguments:
* `-pdb`: RCSB PDB id, if not provided you will be prompted to select one. If it is is in the input/ folder it will be used. Otherwise it will be downloaded from the RCSB.
* `-qc`: Query chain to find extended positions on.
* `-ic`: Interacting chain.
* `-sr`: Solvent radius for extension (default 4.4).
* `-mi`: Interaction cutoff for extension (default 1).
* `-m:` Amino Acids used for extension (default TRP,ARG).
* `-r:` Output file name. 
* `-c:` Number of cores to use in parallel (default 4).
<br />

Add these flags to use certain features. 
*  `-p`: Run Interface Extention in parallel.  
* `--scwrl`: Use scrwl4 to remodel sidechain (default no sidechain remodeling) (-s also works).
* `-qh`: Use c++ qhull (default pyhull). 
* `-nomod`: (FOR TESTING ONLY) Dont use modeller and instead just change the names of the residue in the pdb. 
------

### Output:

![sr](Media/Picture1.png#gh-dark-mode-only)
![sr](Media/Picture2.png#gh-light-mode-only)

------

### Notes:
* Intercaat may not understand pdbs with insertion codes so pdb-tools fixinsert function is run to reformat the pdb. (https://github.com/haddocking/pdb-tools).
* for extended mutants argument(-m) please use upper case three letter amino acid name separated by a comma without spaces.
* if qhull is not downloaded then pyhull will be used (wrapper for qhull). For more info on pyhull see https://github.com/materialsvirtuallab/pyhull. To use qhull downloaded it (http://www.qhull.org/), update the contents of intercaat_config.ini in 'scripts/intercaatmaster' folder  and use the -qh flag in the command line.
* if Scrwl4 isnt working make sure it is properly alliased in your shell config (.bashrc or .zshrc). 

### Extension Parameters:


![sr](Media/params_dark.png#gh-dark-mode-only)
![sr](Media/params_light.png#gh-light-mode-only)
------

Written by Evan Edelstein 
<br />
Please report any questions or complaints to steven.grudman@einsteinmed.org
