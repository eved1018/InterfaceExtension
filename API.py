import os
import sys
from scripts.mutate_model import mutateModel
from scripts.Scwrl import runScwrl4
from scripts.simple_mutate import simple_mutate
from scripts.intercaatInterface import mutantIntercaatRun
from concurrent.futures import ProcessPoolExecutor
import shutil
from scripts.intercaatmaster.intercaatWrapper import intercaat

def interfaceExtention(pdb,pdb_file, query_chain,partner_chain, sr = 3.4, mi =1 , mutants = "ARG,TRP",result_file = False,scwrl = True, qhull= False, nomod = False,parallel = True, cores = 4, intercaat_result= None):
    mutants = mutants.split(",")
    mutants = [i.upper() for i in mutants]
    extended_interface = []
    intercaat_result_changed, positions = intercaatRun(pdb_file, query_chain, partner_chain, sr, mi, qhull, intercaat_result)
    print(positions)
    results = {key: [] for key in positions}
    jobs = [(key, mutAA, pdb, pdb_file, positions, scwrl, qhull, sr, query_chain, partner_chain, nomod) for key in positions for mutAA in mutants if key[:3] != mutAA]
    with ProcessPoolExecutor(max_workers=cores) as exe:
        return_vals = exe.map(parelleljob, jobs)
        for return_val in return_vals:
            rkey, is_exteneded, result = return_val
            results[rkey] = results[rkey] + result
            if is_exteneded and rkey not in extended_interface:
                extended_interface.append(rkey)
    outputWriter(result_file,pdb, pdb_file, query_chain, partner_chain, intercaat_result, intercaat_result_changed, extended_interface, results, positions, sr, mi)
    return extended_interface


def parelleljob(args):
    key, mutAA,pdb,pdb_file, positions, scrwl, qhull, sr, query_chain, partner_chain, nomod = args
    respos = str(positions[key][0])
    mutposition = mutAA + respos
    mutantfile = f"output/{pdb}/mutants/" + mutposition + ".pdb"
    # if nomod:
    #     mutant = simple_mutate(pdb_file, query_chain, respos, key[:3], mutAA, mutantfile)
    # else:
    mutant = mutateModel(pdb_file, respos, mutAA, query_chain, mutantfile, "input/")

    if scrwl:
        mutant = runScwrl4(mutant)
    mutant_interactions = mutantIntercaatRun(mutant, query_chain, partner_chain, mutposition, sr, qhull)
    return key, mutant_interactions, [f"{mutAA} {mutant_interactions}"]

def outputWriter(result_file, pdb,pdb_file ,query_chain, partner_chain, intercaat_result, intercaat_result_changed, extended_interface, results, positions, sr, mi):
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
        outfile.write("\nExtended Interface Position(s): ")
        outfile.write(" ".join(extended_interface))
        outfile.write("\nMutation results: (True if made contact after mutation)\n")
        for i in results:
            j = " ".join(results[i])
            outfile.write(f"{i} {j}\n")
    return

def intercaatRun(in_pdb, query_chain, partner_chain, sr, mi, qhull, intercaat_result):
    positions = {}
    intercaat_result_changed = intercaat(
        in_pdb, query_chain, partner_chain, sr=sr, mi=mi, fp="input/")
    if not intercaat_result:
        print("Error: it is likley the two chains selected do not interact")
        sys.exit()

    for key in intercaat_result_changed:
        if key not in intercaat_result:
            interactions = intercaat_result_changed[key][2]
            resnum = intercaat_result_changed[key][1]
            positions[key] = [resnum, interactions]

    return intercaat_result_changed, positions

