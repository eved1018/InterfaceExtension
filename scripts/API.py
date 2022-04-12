from scripts.pdb_fixinsert import fixInsert
from scripts.CLI import download_pdb
import os
import sys
from scripts.mutate_model import mutateModel
from scripts.Scwrl import runScwrl4
from scripts.simple_mutate import simple_mutate
from scripts.intercaatInterface import intercaatRun, mutantIntercaatRun
from concurrent.futures import ProcessPoolExecutor
import shutil




def interfaceExtention(pdb,query_chain,partner_chain, sr = 3.4, mi =1 , mutants = "ARG,TRP",result_file = False,scwrl = True, qhull= False, nomod = False,parallel = True, cores = 4):
    files = [i for i in os.listdir("input/") if i.endswith(".pdb")]
    mutants = mutants.split(",")
    mutants = [i.upper() for i in mutants]
    for i in mutants:
        if len(i) != 3:
            print("Error please use three letter amino acid code separated by a comma")
            sys.exit()
    if mi > 4:
        print("Error please choose `mi` less than 4")
        sys.exit()
    pdb, pdb_file = pdbManager(pdb, files)
    pdb_file = fixInsert(pdb_file)
    print(files)
    print(pdb_file)
    os.makedirs("output/", exist_ok=True)
    os.makedirs(f"output/{pdb}/mutants/", exist_ok=True)
    extended_interface = []
    intercaat_result, intercaat_result_changed, positions = intercaatRun(pdb, query_chain, partner_chain, sr, mi, qhull)
    results = {key: [] for key in positions}
    jobs = [(key, mutAA, pdb, positions, scwrl, qhull, sr, query_chain, partner_chain, nomod) for key in positions for mutAA in mutants if key[:3] != mutAA]
    with ProcessPoolExecutor(max_workers=cores) as exe:
        return_vals = exe.map(parelleljob, jobs)
        for return_val in return_vals:
            rkey, is_exteneded, result = return_val
            results[rkey] = results[rkey] + result
            if is_exteneded and rkey not in extended_interface:
                extended_interface.append(rkey)
    outputWriter(result_file, pdb_file, query_chain, partner_chain, intercaat_result, intercaat_result_changed, extended_interface, results, positions, sr, mi)
    return extended_interface

def pdbManager(pdb, files):
    pdb_file = pdb + ".pdb"
    if pdb_file not in files:
        download_pdb(pdb, "input/")
    return pdb, pdb_file


def parelleljob(args):
    key, mutAA,pdb, positions, scrwl, qhull, sr, query_chain, partner_chain, nomod = args
    wt_interactions = int(positions[key][1])
    respos = str(positions[key][0])
    mutposition = mutAA + respos
    mutantfile = f"output/{pdb}/mutants/" + mutposition + ".pdb"
    if nomod:
        mutant = simple_mutate(pdb, query_chain, respos, key[:3], mutAA, mutantfile)
    else:
        mutant = mutateModel(pdb, respos, mutAA, query_chain, mutantfile, "input/")

    if scrwl:
        mutant = runScwrl4(mutant)
    mutant_interactions = mutantIntercaatRun(mutant, query_chain, partner_chain, mutposition, sr, qhull)
    return key, mutant_interactions, [f"{mutAA} {mutant_interactions}"]

def outputWriter(result_file, pdb, query_chain, partner_chain, intercaat_result, intercaat_result_changed, extended_interface, results, positions, sr, mi):
    with open(f"output/{pdb}/{result_file}", "a+") as outfile:
        outfile.write("-------------------\n")
        outfile.write(f"Protein: {pdb} qc: {query_chain} ic {partner_chain}")
        outfile.write("\nInterface:\n")
        outfile.write("Res\t#Interactions\n")
        for i in intercaat_result:
            outfile.write(f"{i}\t{intercaat_result[i][2]}\n")
        outfile.write(
            f"\nExtended Interface (solvent radius {sr} | minimum interactions {mi}):\n")
        outfile.write("Res\t#Interactions\n")
        for i in intercaat_result_changed:
            if i in positions:
                outfile.write(f"*{i}\t{intercaat_result_changed[i][2]}\n")
            else:
                outfile.write(f" {i}\t{intercaat_result_changed[i][2]}\n")
        outfile.write("\nExtened Interface Position(s): ")
        outfile.write(" ".join(extended_interface))
        outfile.write("\nMutation results: (True if made contact after mutation)\n")
        for i in results:
            j = " ".join(results[i])
            outfile.write(f"{i} {j}\n")
    shutil.make_archive(f"output/{pdb}/mutants_{pdb}" , 'zip', f"output/{pdb}/mutants")
    shutil.rmtree(f"output/{pdb}/mutants")
    return


