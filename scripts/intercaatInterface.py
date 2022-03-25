from .intercaatmaster.intercaatWrapper import intercaat
import sys


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