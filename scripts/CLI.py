import urllib.request
import os
import argparse
import sys

def CLI():
    # pdb = "6waq.pdb"
    # partner_chain = "A"
    # query_chain = "B"
    files = [i for i in os.listdir("input/") if i.endswith(".pdb")]
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb',required =True, help= f"please choose an input file from {files} or add to input folder")
    parser.add_argument('-qc'  , '--query_chain', help= f"please choose a query chain", default = None)
    parser.add_argument('-ic' ,'--interacting_chain', help= f"please choose the interacting chain", default = None)
    parser.add_argument('-sr' ,'--solvent_radius', help= f"please choose the solvent radius extention (defualt 2.4)", default = 4.4)
    parser.add_argument('-r' ,'--result_file', help= f"result file", default = "result.txt")
    parser.add_argument('-mi' ,'--min_ints', help= f"please choose the minimum interaction (defualt 4)", default = 1)
    parser.add_argument('-s' ,'--scrwl', help= f"use Scrwl4 to remodel sidechains", default = True)
    parser.add_argument('-m' ,'--mutants', help= f"mutants to change to: ex TRP,ARG", default = "TRP,ARG")
    args = parser.parse_args()
    pdb = args.pdb
    pdb_file = pdb + ".pdb"
    partner_chain = args.interacting_chain
    query_chain = args.query_chain
    sr = args.solvent_radius
    result_file = args.result_file
    mi = args.min_ints
    scrwl = args.scrwl
    mutants = args.mutants
    mutants = mutants.split(",")
    if pdb_file not in files:
        download = input(f"Download {pdb} from the RCSB? [y,n]: ")
        if download == "y":
            download_pdb(pdb, "input/")
        else:
            print(f"{pdb} not downloaded please add {pdb} to input folder")
            sys.exit()
    if query_chain is None or partner_chain is None:
        query_chain, partner_chain = getChains(pdb)
    os.makedirs("output/mutants/", exist_ok = True)
    return pdb_file, query_chain, partner_chain, sr, result_file, mi, scrwl, mutants

def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None

def getChains(pdb):
    chains = set()
    with open(f"input/{pdb}.pdb", "r") as pdb_fh:
        for line in pdb_fh:
            if line.startswith("ATOM"):
                chains.add(line[21])
    qc = input(f"Please select query [{chains}]: ")
    if qc not in chains:
        print("chain not found")
        sys.exit()

    chains.remove(qc)
    ic = input(f"Please select interacting chain [{chains}]: ")
    if ic not in chains:
        print("chain not found")
        sys.exit()
    return qc, ic
