from scripts.mutate_model import mutateModel
from scripts.CLI import CLI
from scripts.Scwrl import runScwrl4
from scripts.simple_mutate import simple_mutate
from scripts.intercaatInterface import intercaatRun, mutantIntercaatRun
from concurrent.futures import ProcessPoolExecutor
import shutil


def main(pdb, query_chain, partner_chain, sr, result_file, mi, scrwl, mutants, qhull, nomod, cores):
    extended_interface = []
    intercaat_result, intercaat_result_changed, positions = intercaatRun(
        pdb, query_chain, partner_chain, sr, mi, qhull)
    results = {key: [] for key in positions}
    jobs = [(key, mutAA, pdb, positions, scrwl, qhull, sr, query_chain, partner_chain, nomod) for key in positions for mutAA in mutants if key[:3] != mutAA]    
    print(f"{len(jobs)} jobs")
    with ProcessPoolExecutor(max_workers=cores) as exe:
        return_vals = exe.map(parellelRun, jobs)
        for return_val in return_vals:
            rkey, is_exteneded, result = return_val
            results[rkey] = results[rkey] + result
            if is_exteneded and rkey not in extended_interface:
                extended_interface.append(rkey)  
    outputWriter(result_file, pdb, query_chain, partner_chain, intercaat_result,
                 intercaat_result_changed, extended_interface, results, positions)
    print("extended interface positions: ", extended_interface)
    return extended_interface


def parellelRun(args):
    key, mutAA,pdb, positions, scrwl, qhull, sr, query_chain, partner_chain, nomod = args
    wt_interactions = int(positions[key][1])
    respos = str(positions[key][0])
    mutposition = mutAA + respos
    mutantfile = "output/mutants/" + mutposition + ".pdb"
    if nomod:
        mutant = simple_mutate(pdb, query_chain, respos, key[:3], mutAA, mutantfile)
    else: 
        mutant = mutateModel(pdb, respos, mutAA, query_chain, mutantfile, "input/")

    if scrwl:
        mutant = runScwrl4(mutant)
    mutant_interactions = mutantIntercaatRun(mutant, query_chain, partner_chain, mutposition, sr, qhull)
    return key, mutant_interactions, [f"{mutAA} {mutant_interactions}"]

def singlethreadRun(pdb, query_chain, partner_chain, sr, result_file, mi, scrwl, qhull, nomod, mutants):
    extended_interface = []
    print(pdb, query_chain, partner_chain)
    intercaat_result, intercaat_result_changed, positions = intercaatRun(
        pdb, query_chain, partner_chain, sr, mi, qhull)
    results = {key: [] for key in positions}
    for mutAA in mutants:
        for key in positions:
            wt_res = key[:3]
            if wt_res!= mutAA:
                wt_interactions = int(positions[key][1])
                respos = str(positions[key][0])
                mutposition = mutAA + respos
                mutantfile = "output/mutants/" + mutposition + ".pdb"
                if nomod:
                    mutant = simple_mutate(pdb, query_chain, respos, wt_res,mutAA, mutantfile)
                else: 
                    mutant = mutateModel(pdb, respos, mutAA,query_chain, mutantfile, "input/")
                if scrwl:
                    mutant = runScwrl4(mutant)
                mutant_interactions = mutantIntercaatRun(
                    mutant, query_chain, partner_chain, mutposition, sr, qhull)
                results[key] = results[key] + [f"{mutAA} {mutant_interactions}"]
                if mutant_interactions and key not in extended_interface:
                    extended_interface.append(key)
    outputWriter(result_file, pdb, query_chain, partner_chain, intercaat_result,
                 intercaat_result_changed, extended_interface, results, positions)
    print("extended interface positions: ", extended_interface)
    return extended_interface


def outputWriter(result_file, pdb, query_chain, partner_chain, intercaat_result, intercaat_result_changed, extended_interface, results, positions):
    with open(f"output/{result_file}", "a+") as outfile:
        outfile.write("-------------------\n")
        outfile.write(f"Protein: {pdb} qc: {query_chain} ic {partner_chain}")
        outfile.write("\nInterface:\n")
        outfile.write("Res\t#Interactions\n")
        for i in intercaat_result:
            outfile.write(f"{i}\t{intercaat_result[i][2]}\n")
        outfile.write(
            f"\Extended Interface (solvent radius {sr} | minimum interactions {mi}):\n")
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
    shutil.make_archive(f"output/mutants_{pdb}" , 'zip', "output/mutants")
    shutil.rmtree("output/mutants")
    return


if __name__ == '__main__':
    pdb, query_chain, partner_chain, sr, result_file, mi, scrwl, mutants, qhull,nomod, cores, parallel = CLI()
    if parallel:
        print(f"using {cores} cores:")
        extended_interface = main(pdb, query_chain, partner_chain, sr, result_file, mi, scrwl, mutants, qhull,nomod,  cores)
    else:
        extended_interface = singlethreadRun(pdb, query_chain, partner_chain, sr, result_file, mi, scrwl, qhull,nomod, mutants)
