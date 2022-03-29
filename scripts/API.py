from concurrent.futures import ProcessPoolExecutor
import sys
import subprocess
import os
from modeller import *
from modeller.optimizers import MolecularDynamics, ConjugateGradients
from modeller.automodel import autosched
from pyhull.voronoi import VoronoiTess
import numpy as np

"""
Dependencies:
modeller
numpy 
pyhull

API arguments:
required:
    1) pdb: pdb file
    2) query_chain: chain to mutate
    3) partner_chain: interacting chain

API call will return extended interface. 
"""

def API(pdb, query_chain, partner_chain, sr=3.4, result_file="result.txt", mi=1, scrwl=False, mutants=["ARG,TRP"], qhull=False, modeller= True, cores=1):
    extended_interface = []
    pdb = fixInsert(pdb)
    intercaat_result, intercaat_result_changed, positions = intercaatRun(
        pdb, query_chain, partner_chain, sr, mi, qhull)
    results = {key: [] for key in positions}
    jobs = [(key, mutAA, pdb, positions, scrwl, qhull, sr, query_chain, partner_chain, modeller) for key in positions for mutAA in mutants if key[:3] != mutAA]    
    with ProcessPoolExecutor(max_workers=cores) as exe:
        return_vals = exe.map(parellelRun, jobs)
        for return_val in return_vals:
            rkey, is_exteneded, result = return_val
            results[rkey] = results[rkey] + result
            if is_exteneded and rkey not in extended_interface:
                extended_interface.append(rkey)  
    return extended_interface


def parellelRun(args):
    key, mutAA,pdb, positions, scrwl, qhull, sr, query_chain, partner_chain, modeller = args
    wt_interactions = int(positions[key][1])
    respos = str(positions[key][0])
    mutposition = mutAA + respos
    mutantfile = "output/mutants/" + mutposition + ".pdb"
    if modeller:
        mutant = mutateModel(pdb, respos, mutAA, query_chain, mutantfile, "input/")
    else: 
        mutant = simple_mutate(pdb, query_chain, respos, key[:3], mutAA, mutantfile)
    if scrwl:
        mutant = runScrwl4(mutant)
    mutant_interactions = mutantIntercaatRun(mutant, query_chain, partner_chain, mutposition, sr, qhull)
    return key, mutant_interactions, [f"{mutAA} {mutant_interactions}"]



def intercaatRun(in_pdb, query_chain, partner_chain, sr, mi, qhull):
    positions = {}
    intercaat_result = intercaat(
        in_pdb, query_chain, partner_chain, fp="input/", qhull = qhull)
    intercaat_result_changed = intercaat(
        in_pdb, query_chain, partner_chain, sr=sr, mi=mi, fp="input/")
    if not intercaat_result:
        print("Error: it is likley the two chains selected to not interact")
        sys.exit()

    for key in intercaat_result_changed:
        if key not in intercaat_result:
            interactions = intercaat_result_changed[key][2]
            resnum = intercaat_result_changed[key][1]
            positions[key] = [resnum, interactions]

    return intercaat_result, intercaat_result_changed, positions


def mutantIntercaatRun(mutant_pdb, query_chain, partner_chain, mut_position, sr, qhull):
    intercaat_result = intercaat(mutant_pdb, query_chain, partner_chain, qhull)
    if mut_position in intercaat_result:
        return True
    else:
        return False


def runScrwl4(mutant):
    subprocess.run(f"Scwrl4 -i {mutant} -o {mutant}", shell = True, stdout = subprocess.DEVNULL)
    return mutant


def optimize(atmsel, sched):
    #conjugate gradient
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    #md
    refine(atmsel)
    cg = ConjugateGradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


#molecular dynamics
def refine(atmsel):
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = MolecularDynamics(cap_atom_shift=0.39, md_time_step=4.0,
                           md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600,
                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                         max_iterations=its, equilibrate=equil)
            init_vel = False


#use homologs and dihedral library for dihedral angle restraints
def make_restraints(mdl1, aln):
   rsr = mdl1.restraints
   rsr.clear()
   s = Selection(mdl1)
   for typ in ('stereo', 'phi-psi_binormal'):
       rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
   for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
       rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0,
                spline_dx=0.3, spline_min_points = 5, aln=aln,
                spline_on_site=True)

#first argument
# modelname, respos, restyp, chain, = sys.argv[1:]
def mutateModel(modelname, respos, restyp, chain,filename, folder):

    modelname = folder+modelname
    # log.verbose()
    log.none()

# Set a different value for rand_seed to get a different final model
    env = Environ(rand_seed=-49837)

    env.io.hetatm = True
#soft sphere potential
    env.edat.dynamic_sphere=False
#lennard-jones potential (more accurate)
    env.edat.dynamic_lennard=True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

# Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

# Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')


# Read the original PDB file and copy its sequence to the alignment array:
    mdl1 = Model(env, file=modelname)
    ali = Alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

#set up the mutate residue selection segment
    s = Selection(mdl1.chains[chain].residues[respos])

#perform the mutate residue operation
    s.mutate(residue_type=restyp)
#get two copies of the sequence.  A modeller trick to get things set up
    ali.append_model(mdl1, align_codes=modelname)

# Generate molecular topology for mutant
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])


# Transfer all the coordinates you can from the template native structure
# to the mutant (this works even if the order of atoms in the native PDB
# file is not standard):
#here we are generating the model by reading the template coordinates
    mdl1.transfer_xyz(ali)

# Build the remaining unknown coordinates
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

#yes model2 is the same file as model1.  It's a modeller trick.
    mdl2 = Model(env, file=modelname)

#required to do a transfer_res_numb
#ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
#transfers from "model 2" to "model 1"
    mdl1.res_num_from(mdl2,ali)

#It is usually necessary to write the mutated sequence out and read it in
#before proceeding, because not all sequence related information about MODEL
#is changed by this command (e.g., internal coordinates, charges, and atom
#types and radii are not updated).

    mdl1.write(file=filename+ '.tmp')
    mdl1.read(file=filename+'.tmp')

#set up restraints before computing energy
#we do this a second time because the model has been written out and read in,
#clearing the previously set restraints
    make_restraints(mdl1, ali)

#a non-bonded pair has to have at least as many selected atoms
    mdl1.env.edat.nonbonded_sel_atoms=1

    sched = autosched.loop.make_for_model(mdl1)

#only optimize the selected residue (in first pass, just atoms in selected
#residue, in second pass, include nonbonded neighboring atoms)
#set up the mutate residue selection segment
    s = Selection(mdl1.chains[chain].residues[respos])

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(s)

    s.energy()

    s.randomize_xyz(deviation=4.0)

    mdl1.env.edat.nonbonded_sel_atoms=2
    optimize(s, sched)

#feels environment (energy computed on pairs that have at least one member
#in the selected)
    mdl1.env.edat.nonbonded_sel_atoms=1
    optimize(s, sched)

    opt_energy =s.energy()
    # print("optimal energy: ", opt_energy)
#give a proper name
    mdl1.write(file=filename)

#delete the temporary file
    os.remove(filename+ '.tmp')
    return filename

def run(fhandle, option_list):
    """
    Delete insertion codes (at specific residues).
    By default, removes ALL insertion codes on ALL residues. Also bumps
    the residue numbering of residues downstream of each insertion.
    This function is a generator.
    Parameters
    ----------
    fhandle : a line-by-line iterator of the original PDB file.
    option_list : list
        List of insertion options to act on.
        Example ["A9", "B12"]. An empty list performs the default
        action.
    Yields
    ------
    str (line-by-line)
        The modified (or not) PDB line.
    """

    option_set = set(option_list)  # empty if option_list is empty

    # Keep track of residue numbering
    # Keep track of residues read (chain, resname, resid)
    offset = 0
    prev_resi = None
    seen_ids = set()
    clean_icode = False
    records = ('ATOM', 'HETATM', 'ANISOU', 'TER')
    for line in fhandle:
        if line.startswith(records):
            res_uid = line[17:27]  # resname, chain, resid, icode
            id_res = line[21] + line[22:26].strip()  # A99, B12
            has_icode = line[26].strip()  # ignore ' ' here

            # unfortunately, this is messy but not all PDB files follow a nice
            # order of ' ', 'A', 'B', ... when it comes to insertion codes..
            if prev_resi != res_uid:  # new residue
                # Does it have an insertion code
                # OR have we seen this chain + resid combination before?
                # #2 to catch insertions WITHOUT icode ('A' ... ' ' ... 'B')
                if (has_icode or id_res in seen_ids):
                    # Do something about it
                    # if the user provided options and this residue is in them
                    # OR if the user did not provide options
                    if (not option_set) or (id_res in option_set):
                        clean_icode = True
                    else:
                        clean_icode = False
                else:
                    clean_icode = False

                prev_resi = res_uid

                if id_res in seen_ids:  # offset only if we have seen this res.
                    offset += 1

            if clean_icode:  # remove icode
                line = line[:26] + ' ' + line[27:]

            # Modify resid if necessary
            resid = int(line[22:26]) + offset
            line = line[:22] + str(resid).rjust(4) + line[26:]
            seen_ids.add(id_res)

            # Reset offset on TER
            if line.startswith('TER'):
                offset = 0
        yield line


def fixInsert(pdb_file):
    has_icode = False
    with open(f"input/{pdb_file}", "r") as pdb_fh:
        for line in pdb_fh:
            if line.startswith("ATOM"):
                if line[26] != ' ':
                    has_icode = True

    if has_icode:
        # Check Input
        pdbfh = open(f"input/{pdb_file}", 'r')
        option_list = []
        # Do the job
        new_pdb = run(pdbfh, option_list)

        with open(f"input/fixed_{pdb_file}", "w+") as _file:
            for line in new_pdb:
                _file.write(line)
        pdbfh.close()
        return f"fixed_{pdb_file}"
    else:
        return pdb_file

def simple_mutate(pdb, qc, respos, wt_res,mutAA, mutantfile):
    mutated_structure = []
    mutposition = mutAA + str(respos)
    with open(f"input/{pdb}", "r") as _pdbfh:
        for line in _pdbfh:
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                s_chain = line[21].strip()
                s_resi = line[22:26].strip()
                s_resn = line[17:20].strip()
                if s_chain == qc and s_resi == respos:
                    if s_resn == wt_res:
                        line = line[0:17]+mutAA+line[20:]
                        mutated_structure.append(line)
                    else:
                        print("Error:")
                        return None
                else:
                    mutated_structure.append(line)
            else:
                mutated_structure.append(line)
    with open(mutantfile, "w+") as _mutfile:
        for line in mutated_structure:
            _mutfile.write(line)
    return mutantfile

def intercaat(pdb: str, qc: str, ic: str, mi: int = 4,di: str = "yes",cc: str = "yes", sr: float = 1.4, vi = [], fp: str = "./", qhull = False):
    arg1 = pdb
    arg2 = qc.split(',')
    arg3 = ic.split(',')
    arg4 = int(mi)
    arg5 = di
    arg6 = cc
    arg7 = float(sr)
    arg8 = vi
    arg9 = fp
    if arg8 != []:
        arg8 = arg8.split(',')
    try:
        pdb  = parse(arg1, (arg2+arg3 + arg8), arg9)
    except FileNotFoundError:
        sys.exit("filename not found")
    coordinates = []
    match = []

    # Parse PDB file into list of PDB lines containing atoms and hetatms excluding hydrogen and water.
    for line in pdb:
        coordinates.append([line[8], line[9], line[10]])

    # Creates 3D voronoi diagram and returns indices of neighboring atoms

    
    contacts = voroPyhull(coordinates)
    # Creates a list (pdbAtomClass) that contains classes for all atoms analayzed
    pdbAtomClass = appendAtomClasses(pdb)

    # Creates list of all atomic interactions
    count1 = 0
    count2 = 0
    for buddies in contacts:
        while count2 < len(buddies):
            # XYZ = atom, X, Y, Z
            XYZ1 = [pdb[count1][2][0:2], float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])]
            XYZ2 = [pdb[buddies[count2]][2][0:2], \
            float(pdb[buddies[count2]][8]), float(pdb[buddies[count2]][9]), float(pdb[buddies[count2]][10])]
            # Returns the distance between two atoms and min distance reqeuired for solvent molecule
            Ad, Vd = inter(XYZ1,XYZ2,arg7)
            # Finds class of atom1 and atom2
            class1 = pdbAtomClass[count1]
            class2 = pdbAtomClass[buddies[count2]]
            # atom = residue, residue #, chain, atom
            atom1 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[count1][4], pdb[count1][6], pdb[count1][5], pdb[count1][2])
            atom2 = '{0:<3} {1:>5} {2} {3:<4}'.format(pdb[buddies[count2]][4], pdb[buddies[count2]][6],\
                     pdb[buddies[count2]][5], pdb[buddies[count2]][2])
            Line = '{0} | {1} | {2:<4} |    {3}   {4}'.format(atom1, atom2, str(round(Ad,2)), str(class1), str(class2))
            # Only appends if classes are compatible or the user inputs that class compatibility does not matter
            # if class compatibility is unknown, the atomic interaction will be shown
            if compatible(class1,class2) == True or arg6 == 'no':
                # line 1: appends if distance between atoms is within bound and accounts for occupancy < 1
                # line 2: appends only for specific chain
                if Ad < Vd and any((atom1 + atom2) in sub for sub in match) == False \
                and pdb[count1][5] == arg2[0]:
                    # appends all residues from specified chain
                    if len(arg2) == 1:
                        # appends only against queried neighbor chains
                        if pdb[buddies[count2]][5] in arg3:
                            match.append(Line)
                    #appends only specified residues from specified chain
                    elif pdb[count1][6] in arg2:
                        # appends only against queried neighbor chains
                        if pdb[buddies[count2]][5] in arg3:
                            match.append(Line)
            count2 += 1
        count2 = 0
        count1 += 1
    # Filters the list generated above (match) based on arg2 and arg4
    # Also displays interactions matrix based on arg5
    out1, out2, out3 = filterMatch(match, pdb, arg2, arg4, arg5)
    return to_dict(out2, out3)

def to_dict(newInteractionRes, newInteraction):
    return {f"{i}{l[0]}": [i,l[0],l[1]] for i, l in zip(newInteractionRes, newInteraction)}

def voroPyhull(points):
	points = [[float(i) for i in j] for j in points] 	
	v = VoronoiTess(np.array(points))
	contacts = v.ridges
	contacts = [i for i in contacts]
	contacts = np.array(contacts)
	count1 = 0
	count2 = 0
	hold = []
	neighbors = []
	while count1 < len(points):
		while count2 < len(contacts):
			if contacts[count2,0] == count1:
				hold.append(contacts[count2,1])
			if contacts[count2,1] == count1:
				hold.append(contacts[count2,0])
			count2 += 1
		hold = sorted(hold)
		neighbors.append(hold)
		hold = []
		count2 = 0
		count1 += 1
	return neighbors

def aClass(coordinates, atom, noPop):
	"""
	Assigns each atom a class. Class definitions were copied from CSU
	Always assume nitrogen is class III. it can sometimes be class I but it is unclear when
	DOI: 10.1093/bioinformatics/15.4.327
	arg1:   (list of list of floats) list of coordinates of each atom in a residue. Sometimes the coordinates of the atom
	in the next residue. For example when assigning class to amino acids, you need to consider the N in the next residue.
	arg2:   (list of strings) list of atoms types for each atom in a residue
	arg3:	(boolean) Decides whether to return the atom in the next residue or not
	return: (dict) returns a dictionary with the atom index and the class it belongs to for
				   a single residue
	return: (string) '?' if unsure of class type.
	NUMBER  ATOM CLASS
	I       Hydrophilic
	II      Acceptor
	III     Donor
	IV      Hydrophobic
	V       Aromatic
	VI      Neutral
	VII     Neutral-donor
	VIII    Neutral-acceptor
	"""
	count1 = 0
	count2 = 0
	distMatrix = np.zeros((len(coordinates[:,1]), len(coordinates[:,1])))

	# create distance matrix (distance of residue atom to every other residue atom)
	# only keep distances less than 1.6 angstroms: covalent bond.
	for i in coordinates:
		while count1 < len(coordinates[:,1]):
			hold = dist(i[0], i[1], i[2], coordinates[count1,0], coordinates[count1,1], coordinates[count1,2])
			if hold < 1.7:
				distMatrix[count1, count2] = hold
			count1 += 1
		count1 = 0
		count2 += 1

	# Determines if residue is planar
	pla = False
	if np.count_nonzero(np.tril(distMatrix)) >= len(coordinates):
		pla = planar(coordinates)

	# Assigns atom class to oxygen atoms and as a default assings all carbons to class 4
	count1 = 0
	count2 = 0
	atomClasses = {}
	while count1 < len(atom):
		if atom[count1][0:1] == 'O':
			atomClasses[count1] = '1'
			if np.count_nonzero(distMatrix[count1,:]) == 1:
				while count2 < len(distMatrix):
					if distMatrix[count1,count2] != 0 and distMatrix[count1,count2] < 1.3:
						atomClasses[count1] = '2'
					count2 += 1
			if np.count_nonzero(distMatrix[count1,:]) > 1:
				atomClasses[count1] = '2'
			# The program incorrectly assigns atom O3' to class 1. This atom is specifically found in DNA and RNA
			# It assigns O3' to class 1 because the query residue does not contain the phosphate from the next residue
			# I could not think of a simple way to do this.
			# Can be fixed by upgrading function planar to be able to look at more than 1 ring at once on a single residue
			if atom[count1] == 'O3\'':
				atomClasses[count1] = '2'
		if atom[count1][0:1] == 'C':
			atomClasses[count1] = '4'
		count1 += 1
		count2 = 0

	# Assigns atom class to sulfur, potassium, fluorine atoms
	# Assigns atom class to chlorine, bromine, and iodine atoms
	# Assigns atom class to nitrogen.
	count1 = 0
	while count1 < len(atom):
		if atom[count1][0:1] == 'S' or atom[count1][0:1] == 'P' or atom[count1][0:1] == 'F':
			atomClasses[count1] = '6'
		if atom[count1][0:2] == 'Cl' or atom[count1][0:2] == 'Br' or atom[count1][0:1] == 'I' or\
		   atom[count1][0:2] == 'CL' or atom[count1][0:2] == 'BR':
			atomClasses[count1] = '4'
		if atom[count1][0:1] == 'N':
			atomClasses[count1] = '3'
		count1 += 1

	# Assings '?' to all atoms not givin a class
	count1 = 0
	if len(atomClasses) < len(atom):
		while count1 < len(atom):
			try:
				atomClasses[count1]
			except:
				atomClasses[count1] = '?'
			count1 += 1

	# Assigns atom class to carbon
	count1 = 0
	count2 = 0
	holdClass= []
	while count1 < len(atom):
		if atom[count1][0:1] == 'C' and (atom[count1][0:2] != 'CL' or atom[count1][0:2] != 'Cl'):
			while count2 < len(distMatrix):
				if distMatrix[count1,count2] != 0:
					holdClass.append(atomClasses[count2])
				count2 += 1
			if holdClass.count('3') == 1:
				atomClasses[count1] = '7'
			if holdClass.count('2') == 1:
				atomClasses[count1] = '8'
			if holdClass.count('1') == 1:
				atomClasses[count1] = '6'
			if holdClass.count('2') + holdClass.count('3') > 1:
				atomClasses[count1] = '6'
		holdClass = []
		count1 += 1
		count2 = 0

	# Assings atom class to planar carbon atoms
	count1 = 0
	if pla != False:
		while count1 < len(pla):
			if atom[pla[count1]][0:1] == 'C':
				atomClasses[pla[count1]] = '5'
			count1 += 1

	# Decides whether to return the atom in the next residue or not
	if noPop == False:
		atomClasses.pop(len(atomClasses) - 1)

	return atomClasses


def appendAtomClasses(pdb):
	"""
	Calls aClass to append all the atoms classes of each residue to a list (pdbAtomClass)
	arg1:   (list of lists of strings) parsed pdb file being analyzed
	return: (list) The atom classes of all atoms in a pdb file. Should be same length as pdb file
	"""
	pdb.append(['','','','','','','','',1,1,1])
	count1 = 1
	count2 = 0
	dictHold = {}
	pdbAtomClass = []
	classCoordinates = [[float(pdb[0][8]), float(pdb[0][9]), float(pdb[0][10])]]
	noPop = False
	atm = [pdb[0][2]]
	while count1 < len(pdb):
		# Second half of if statent is needed because the function planar can only handle a residue with one ring
		# Since DNA/RNA has at least two rings, the residue must be seperated into pieces, the backbone, and the nitrogenous base
		if pdb[count1 - 1][6] == pdb[count1][6] and \
		(('\'' in pdb[count1 - 1][2] and '\'' in pdb[count1][2]) or ('\'' not in pdb[count1 - 1][2] and '\'' not in pdb[count1][2]) \
		or ('\'' not in pdb[count1 - 1][2] and '\'' in pdb[count1][2])):
			classCoordinates.append([float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])])
			atm.append(pdb[count1][2])
		else:
			# For DNA, does not include the phosphate in the next residue with the nitrogenous base.
			if pdb[count1][2] != 'P':
				classCoordinates.append([float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])])
				atm.append(pdb[count1][2])
			else:
				noPop = True
			dictHold = aClass(np.array(classCoordinates),atm, noPop)
			while count2 < len(dictHold):
				try:
					pdbAtomClass.append(int(dictHold[count2]))
				except:
					pdbAtomClass.append(dictHold[count2])
				count2 += 1
			dictHold = {}
			classCoordinates = [[float(pdb[count1][8]), float(pdb[count1][9]), float(pdb[count1][10])]]
			atm = [pdb[count1][2]]
			noPop = False
		count1 += 1
		count2 = 0
	pdb.pop(len(pdb)-1)
	return pdbAtomClass

def parse(pdb_filename, include, dir):
	'''
	Parse PDB file into list of artibutes with no spaces
	arg1:   (list of strings) line seperated pdb file
	arg2:   (list of strings) list of chains to include for voronoi calculation
	arg3:	(string) path to directory that hold PDB file
	return: (list of list of strings) PDB lines containing atoms and hetatms. 
									  Exludes hydrogen and water atoms. Only includes atoms with highest occupancy
	'''
	file = dir + pdb_filename
	fh = open(file)
	atomTemp = []
	atom = []
	listTemp = []
	count = 0

	for line in fh:
		if 'MODEL' in line[0:6] and '2' in line[12:16]:
			break
		if ('ATOM' in line[0:6] or 'HETATM' in line[0:6]) and line[21:22] in include \
		and 'H' not in line[76:78] and 'HOH' not in line[17:20] and 'NA' not in line[17:20]\
		and (line[16:17] == ' ' or line[16:17] == 'A'):
			atomTemp.append(line)
	fh.close()

	for l in atomTemp:
		listTemp = [l[0:6]] + [l[6:11]] + [l[12:16]] + [l[16:17]] \
		+ [l[17:20]] + [l[21:22]] + [l[22:26]] + [l[26:27]] \
		+ [l[30:38]] + [l[38:46]] + [l[46:54]] + [l[54:60]] + [l[60:66]]
		while count < len(listTemp):
			listTemp[count] = listTemp[count].strip()
			count += 1
		atom.append(listTemp)
		count = 0
	return atom

def filterMatch(match, pdb, arg2, arg4, arg5):
	"""
	Filters the list generated above (match) based on arg2 and arg4
	Also displays interactions matrix based on arg5
	arg1:   (list of list of strings) All atomic interactions
	arg2:   (list of list of strings) parsed pdb file being analyzed
	arg3:   (list of string(s)) Query chains and optional residues
	arg4:   (int) minimum number of atomic interactions
	arg5:   (string) option to display interaction matrix. 'yes' or 'no'
	return: (list of list of strings) filterd match
	"""
	# parse match (the list generated in the last for loop)
	l = []
	hold = []
	res = []
	for line in match:
		l = line.split(' ')
		for element in l:
			if len(element) > 0:
				hold.append(element)
		res.append(hold)
		hold = []

	# need to add empty list to the end of res for next loop to work
	res.append(['','','','','','','','','','','','','',''])
	res = np.array(res)

	# Sorts arg2. Need to turn into int to sort. Also remove chain from arg2
	count4 = 0
	# if arg2 only contains chain identifier
	if len(arg2) == 1:
		for AA in pdb:
			if AA[5] == arg2[0] and int(AA[6]) not in arg2:
				arg2.append(int(AA[6]))
		arg2.remove(arg2[0])
	# if arg2 contains the chain identifier and specific residues
	else:
		arg2Hold = []
		for AA in arg2[1:]:
			arg2Hold.append(int(AA))
		arg2 = arg2Hold
	arg2 = sorted(arg2)
	# turns values in arg2 back to strings
	while count4 < len(arg2):
		arg2[count4] = str(arg2[count4])
		count4 += 1

	# obtain matrix of # of interactions/residue
	# first column: residue number, second column: # of interactions
	# count1 loops through while loop, count2 is an index for arg2, count3 counts # of interactions
	count1 = 1
	count2 = 0
	count3 = 1
	interactions = []
	while count1 < len(res):
		if res[count1 -1][1] != res[count1][1]:
			interactions.append([int(res[count1 - 1][1]), count3])
			count2 += 1
			count3 = 1
		else:
			count3 +=1
		count1 += 1
	interactions = np.array(interactions)

	# Creates new list of matches that have min # of interactions
	count1 = 0
	count2 = 0
	newMatch = []
	while count1 < len(res):
		while count2 < len(interactions):
			if res[count1][1] == str(int(interactions[count2][0])):
				if interactions[count2][1] >= arg4:
					newMatch.append(match[count1])
			count2 += 1
		count2 = 0
		count1 += 1

	# Prints a matrix of interactions/residue
	newInteractions = []
	count1 = 0
	if arg5 == 'yes':
		for i in interactions:
			if i[1] >= arg4:
				newInteractions.append(list(i))
			count1 += 1
		newMatch.append('--------------------------------')
		newInteractionRes = []
		count1 = 0
		while count1 < (len(newMatch)-1):
			if newMatch[count1][5:9] != newMatch[count1+1][5:9]:
				newInteractionRes.append(newMatch[count1][0:3])
			count1 += 1
		count1 = 0
		# print('Res #   Interactions')
		# while count1 < len(newInteractions):
			# print('{0} {1:<5}  {2:<4}'.format(newInteractionRes[count1], newInteractions[count1][0], newInteractions[count1][1]))
			# count1 += 1
		newMatch.pop(len(newMatch)-1)

	return newMatch, newInteractionRes, newInteractions

def planar(z):
	'''
	Determines if a residue contains a ring and if that ring is planar
	If a residue has more than 1 ring, function assumes both rings are planar
	arg1:   (np matrix of 3 floats) xyz coordinates of every atoms in a single residue
	return: (list of ints) list of the indices of atoms that are planar in a residue
	return: False if the residue has a ring but it is not planar (ex. proline)
	'''
	from math import ceil

	count1 = 0
	count2 = 0
	distMatrix = np.zeros((len(z[:,1]), len(z[:,1])))
	holdDistMatrix1 = np.zeros((len(z[:,1]), len(z[:,1])))
	holdDistMatrix2 = np.zeros((len(z[:,1]), len(z[:,1])))

	# create distance matrix (distance of residue atom to every other residue atom)
	# only keep distances less than 1.6 angstroms: covalent bond.
	for i in z:
		while count1 < len(z[:,1]):
			hold = dist(i[0], i[1], i[2], z[count1,0], z[count1,1], z[count1,2])
			if hold < 1.6:
				distMatrix[count1, count2] = hold
			count1 += 1
		count1 = 0
		count2 += 1

	holdDistMatrix1 = distMatrix

	moreThanOneRing = False
	# If a residue has more than 1 ring the program assumes both are planar
	if np.count_nonzero(np.tril(distMatrix)) > len(z):
		moreThanOneRing = True

	# If the below statement is True, the residue does not contain any rings
	if np.count_nonzero(np.tril(distMatrix)) < len(z):
		return False

	# iterativley remove atoms unitl all remaining atoms have at least two neighbots
	# determines which atoms are in a ring
	count1 = 0
	count2 = 0
	while np.array_equal(holdDistMatrix1, holdDistMatrix2) == False:
		holdDistMatrix2 = np.copy(holdDistMatrix1)
		while count1 < len(holdDistMatrix1):
			if np.count_nonzero(holdDistMatrix1[count1,:]) == 1:
				while count2 < len(holdDistMatrix1):
					if holdDistMatrix1[count1,count2] != 0:
						holdDistMatrix1[count1,count2] = 0
						holdDistMatrix1[count2,count1] = 0
					count2 += 1
			count1 += 1
			count2 = 0
		count1 = 0

	# the indices of the matrix represent each atom in a residue
	# only need lower left half of matrix since top right half is redundant
	bondMatrix = np.tril(holdDistMatrix1)

	# Obtain indices of atoms in the ring
	count1 = 0
	count2 = 0
	indices = []
	while count1 < len(bondMatrix):
		while count2 < len(bondMatrix):
			if bondMatrix[count1][count2] > 0:
				indices.append([count1,count2])
			count2 += 1
		count2 = 0
		count1 += 1

	# Creates lists of 4 atoms to check dihedral angle of
	count1 = 0
	hold1 = 0
	hold2 = 0
	dihedralAtoms = []
	while count1 < ceil(len(indices)/2):
		for index in indices:
			if index == indices[count1]:
				continue
			elif indices[count1][0] in index or indices[count1][1] in index:
				if index[0] == indices[count1][0]:
					hold1 = index[1]
				if index[1] == indices[count1][0]:
					hold1 = index[0]
				if index[0] == indices[count1][1]:
					hold2 = index[1]
				if index[1] == indices[count1][1]:
					hold2 = index[0]
		dihedralAtoms.append([hold1,indices[count1][0],indices[count1][1],hold2])
		count1 += 1

	# Checks the dihedral angle and appends angle to dihedralAngle
	count1 = 0
	dihedralAngle = []
	while count1 < len(dihedralAtoms):
		coordinatesHold = np.array([[z[dihedralAtoms[count1][0]][0], z[dihedralAtoms[count1][0]][1], z[dihedralAtoms[count1][0]][2]], \
									[z[dihedralAtoms[count1][1]][0], z[dihedralAtoms[count1][1]][1], z[dihedralAtoms[count1][1]][2]], \
									[z[dihedralAtoms[count1][2]][0], z[dihedralAtoms[count1][2]][1], z[dihedralAtoms[count1][2]][2]], \
									[z[dihedralAtoms[count1][3]][0], z[dihedralAtoms[count1][3]][1], z[dihedralAtoms[count1][3]][2]]])
		dihedralAngle.append(dihe(coordinatesHold))
		count1 += 1

	# if the average of the dihedral angles is less that 3 returns list of planar atoms
	planarAtoms = []
	if sum(dihedralAngle)/len(dihedralAngle) < 3 or moreThanOneRing == True:
		for atomIndex in indices:
			planarAtoms.append(atomIndex[0])
			planarAtoms.append(atomIndex[1])
		return list(set(planarAtoms))
	else:
		return False

def inter(atom1, atom2, solv):
	'''
	Takes two atoms are returns the distance between them and the minimum
	distance required to fit a solvent molecule in between the two atoms.
	Arbitrary default van der waal radius = 1.8
	Radii taken from doi: 10.1021/j100785a001
	arg1:   (list of string and 3 floats) atom1 [atom, x coordinate, y coordinate, z coordinate]
	arg2:   (list of string and 3 flaots) atom2 [atom, x coordinate, y coordinate, z coordinate]
	arg3:   (float) (optional) solvent van der waal radius
	arg1 or arg2 example: ['C' 20.32 19.45 17.53]
	return: (float) distance between two atoms
	return: (float) minimum distance required to fit a sovlent molecule between the two atoms
	'''
	rwS = 1.8
	rwO = 1.52
	rwN = 1.55
	rwC = 1.8
	rwI = 2.01
	rwB = 1.84
	rwF = 1.4
	rwCl = 1.73
	rwSolvent = float(solv)
	atom1rw = 1.8
	atom2rw = 1.8

	if atom1[0][0:1] == 'C':
		atom1rw = rwC
	if atom1[0][0:1] == 'N':
		atom1rw = rwN
	if atom1[0][0:1] == 'S':
		atom1rw = rwS
	if atom1[0][0:1] == 'O':
		atom1rw = rwO
	if atom1[0][0:1] == 'B':
		atom1rw = rwB
	if atom1[0][0:1] == 'F':
		atom1rw = rwF
	if atom1[0][0:1] == 'I':
		atom1rw = rwI
	if atom1[0][0:2] == 'CL' or atom2[0][0:2] == 'Cl':
		atom1rw = rwCl

	if atom2[0][0:1] == 'C':
		atom2rw = rwC
	if atom2[0][0:1] == 'N':
		atom2rw = rwN
	if atom2[0][0:1] == 'S':
		atom2rw = rwS
	if atom2[0][0:1] == 'O':
		atom2rw = rwO
	if atom2[0][0:1] == 'B':
		atom2rw = rwB
	if atom2[0][0:1] == 'F':
		atom2rw = rwF
	if atom2[0][0:1] == 'I':
		atom2rw = rwI
	if atom2[0][0:2] == 'CL' or atom2[0][0:2] == 'Cl':
		atom2rw = rwCl

	atom1atom2Dist = dist(atom1[1],atom1[2],atom1[3],atom2[1],atom2[2],atom2[3])
	vanDerWaalDist = atom1rw + atom2rw + 2*rwSolvent
	return(atom1atom2Dist, vanDerWaalDist)

def dihe(p):
	'''
	Obtains the dihedral angle from the coordinates of 4 atoms
	arg1:   (np matrix of floats) contains the X, Y, and Z coordinates of 4 atoms
	return: The dihedral angle between 4 points
	'''
	# obtain 3 vectors from 4 points
	vectors = p[:-1] - p[1:]
	vectors[0] *= -1
	# calculate normal to two vectors
	normalVec = np.zeros((2,3))
	normalVec[0,:] = np.cross(vectors[0],vectors[1])
	normalVec[1,:] = np.cross(vectors[2],vectors[1])
	# normalize vectors
	normalVec /= np.sqrt(np.einsum('...i,...i', normalVec, normalVec)).reshape(-1,1)
	# return dihedrall angle
	return np.degrees(np.arccos( normalVec[0].dot(normalVec[1])))

def dist(x1, y1, z1, x2, y2, z2):
	'''
	Calculates distance between two points in 3D
	arg1:   (float) x coordinate of atom 1
	arg2:   (float) y coordinate of atom 1
	arg3:   (float) z coordinate of atom 1
	arg4:   (float) x coordinate of atom 2
	arg5:   (float) y coordinate of atom 2
	arg6:   (float) z coordinate of atom 2
	return: (float) Distance between two atoms
	'''
	D = ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)**0.5
	return D


def compatible(c1,c2):
	'''
	Determines if two classes are compatible according to CSU
	Table 1 from DOI: 10.1093/bioinformatics/15.4.327
	arg1:   (int) class of atom 1
	arg2:   (int) class of atom 2
	return: (boolean) True/False
	'''
	if c1 == '?' or c2 == '?':
		return True

	compatibleMatrix = np.array([[1,1,1,0,1,1,1,1],\
								 [1,0,1,0,1,1,1,0],\
								 [1,1,0,0,1,1,0,1],\
								 [0,0,0,1,1,1,1,1],\
								 [1,1,1,1,1,1,1,1],\
								 [1,1,1,1,1,1,1,1],\
								 [1,1,0,1,1,1,0,1],\
								 [1,0,1,1,1,1,1,0]])

	if compatibleMatrix[c1-1][c2-1] == 0:
		return False
	if compatibleMatrix[c1-1][c2-1] == 1:
		return True