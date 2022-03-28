### Interface Extension:

![output](Media/ctla4cd80.png)

------

 
### Overview:

This program uses Intercaat (https://pubmed.ncbi.nlm.nih.gov/34499117/) to identify pharmacological relevant residues at the periphery of the protein interface. This is accomplished by defining the 'wt interface', expanding that interface by modulating the minimum interaction cutoff and/or the solvent radius, mutating the expanded positions to larger mutants (ARG and TRP by default) and finally assessing whether the mutant reaches the "wt interface" threshold.

------


### Installation and Dependencies:
* Requires python 3.6 or later. Tested on python 3.9. 
* Requires Anocanda to install Modeller easily. I recommend MiniConda. https://docs.conda.io/en/latest/miniconda.html
1. Clone repo:
```sh
   git clone https://github.com/eved1018/InterfaceExtention
```
2. Download Modeller: https://salilab.org/modeller/download_installation.html
    * for Conda enviroment:
```sh 
        conda config --add channels salilab
        conda install modeller
```
3. Download python dependencies (pyhull, scipy, numpy):
```sh
    pip install -r requirements.txt 
```
4. Download SCRWL4 (optional): http://dunbrack.fccc.edu/lab/scwrl
5. Download qhull (optional): http://www.qhull.org/

------

### Usage:
1. Move to repo:
```sh
cd InterfaceExtension/
```
2. Run code:
```sh
python main.py 
```
------

### Command Line Arguments:
* `-pdb`: RCSB PDB id, if not provided you will be prompted to select one. If it is is in the input/ folder it will be used. Otherwise it will be downloaded from the RCSB.
* `-qc`: Query chain to find extended positions on.
* `-ic`: partner chain.
* `-sr`: solvent radius for extension (default 4.4).
* `-mi`: Interaction cutoff for extension (default 1).
* `-m:` Amino Acids used for extension (default TRP,ARG).
* `-r:` Output file name. 
* `-c:` Number of cores to use in parallel (default 0 ie. single threaded)
<br />

Add these flags to use Scrwl4 or qhull.
* `-s`: Use scrwl4 to remodel sidechain (default no sidechain remodeling).
* `-qh`: Use c++ qhull (default pyhull). 
------

### Output:
![output](Media/results.png)

------

### Notes:
* Intercaat may not understand pdbs with insertion codes so pdb-tools fixinsert function is run to reformat the pdb. (https://github.com/haddocking/pdb-tools).
* for extended mutants argument(-m) please use upper case three letter amino acid name separated by a comma without spaces.
* if qhull is not downloaded then pyhull will be used (wrapper for qhull). For more info on pyhull see https://github.com/materialsvirtuallab/pyhull. To use qhull downloaded it (http://www.qhull.org/),  update the content intercaat_configs.ini file in 'scripts/intercaatmaster' and use the -qh flag in the command line.

------


Written by Evan Edelstein 
<br />
Please report any questions or complaints to steve.grudman@einsteinmed.org

