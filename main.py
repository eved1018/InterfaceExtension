from scripts.mutate_model import mutateModel
from scripts.CLI import CLI
from scripts.Scrwl import runScrwl4
from scripts.intercaatInterface import intercaatRun, mutantIntercaatRun


def main(pdb, query_chain, partner_chain, sr, result_file, mi, scrwl, mutants):
    extened_interface = []
    print(pdb, query_chain, partner_chain)
    intercaat_result, intercaat_result_changed, positions = intercaatRun(
        pdb, query_chain, partner_chain, sr, mi)
    results = {key: [] for key in positions}
    for mutAA in mutants:
        for key in positions:
            if key[:3] != mutAA:
                wt_interactions = int(positions[key][1])
                respos = str(positions[key][0])
                mutposition = mutAA + respos
                mutantfile = "output/mutants/" + mutposition + ".pdb"
                mutant = mutateModel(pdb, respos, mutAA,
                                     query_chain, mutantfile, "input/")
                if scrwl:
                    mutant = runScrwl4(mutant)
                mutant_interactions = mutantIntercaatRun(
                    mutant, query_chain, partner_chain, mutposition, sr)
                print("\n\n", key, wt_interactions,
                      mutant_interactions, "\n\n")
                results[key] = results[key] + \
                    [f"{mutAA} {wt_interactions} {mutant_interactions}"]
                if mutant_interactions and key not in extened_interface:
                    extened_interface.append(key)
    outputWriter(result_file, pdb, query_chain, partner_chain, intercaat_result,
                 intercaat_result_changed, extened_interface, results, positions)
    return


def outputWriter(result_file, pdb, query_chain, partner_chain, intercaat_result, intercaat_result_changed, extened_interface, results, positions):
    with open(f"output/{result_file}", "a+") as outfile:
        outfile.write("-------------------\n")
        outfile.write(f"Protein: {pdb}\tqc: {query_chain}\tic {partner_chain}")
        outfile.write(f"\nW.T. interface:\n")
        outfile.write("Res #   Interactions\n")
        for i in intercaat_result:
            outfile.write(f"{i}\t{intercaat_result[i][2]}\n")
        outfile.write(
            f"\ninterface (solvent radius {sr} | minimum interactions {mi}\t):\n")
        outfile.write("Res #   Interactions\n")
        for i in intercaat_result_changed:
            if i in positions:
                outfile.write(f"*{i}\t{intercaat_result_changed[i][2]}\n")
            else:
                outfile.write(f" {i}\t{intercaat_result_changed[i][2]}\n")
        outfile.write("\npotential extened interface positions: ")
        outfile.write(" ".join(extened_interface))
        outfile.write("\n all results:\n")
        for i in results:
            j = " ".join(results[i])
            outfile.write(f"{i} {j}\n")
    return


if __name__ == '__main__':
    pdb, query_chain, partner_chain, sr, result_file, mi, scrwl, mutants = CLI()
    main(pdb, query_chain, partner_chain, sr, result_file, mi, scrwl, mutants)
