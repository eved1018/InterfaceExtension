
<br />
<br />
<div align="center">
  <a href="https://github.com/eved1018/InterfaceExtention">
    <img src="Media/ctla4cd80.png" alt="Logo" width="250" height="120">
  </a>

  <h3 align="center">Interface Extention</h3>

  </p>
</div>

### Overview:

This program uses Intercaat to identify pharmacologically relevent residues at the perirephary of the strictly defined interface. This is accomplished by defining the 'wt interface', expanding that interface by modulating the minimum interaction cutoff and/or the solvent radius, mutating the expanded positions to larger mutants (ASP and TRP by defualt) and finally assesing wether the mutant reaches the "wt interface" threshold.


### Installation and Depndencies:
1. Clone repo:
```sh
   git clone https://github.com/eved1018/InterfaceExtention
```
2. Download Modeller: https://salilab.org/modeller/download_installation.html
```sh 
    conda config --add channels salilab
    conda install modeller
```
3. Download SCRWL4 (optional): http://dunbrack.fccc.edu/lab/scwrl


### Usage:

1. Run code:
```sh
python main.py -pdb `pdb file name` -qc `chain id` -ic `chain id` -sr `solvent radius expansion`
-mi `minimum interaction decrease` -r `results file` -s `use scrwl`
```

### Output:

Notes:
* Intercaat may not understand pdbs with insertion codes use pdb-tools fixinsert function to refromat the pdb. (https://github.com/haddocking/pdb-tools)
* if the pdb file is not found in the input folder a prompt will ask for permision to download it from the RCSB.


Written by Evan Edelstein
Please report any questions or complaints to steve.grudman@einsteinmed.org

<br />
<br />
<div id="Footer"></div>

[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/evan-edelstein/
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
