from API import interfaceExtention
import pandas as pd
import matplotlib.pyplot as plt
from scripts.intercaatmaster.intercaatWrapper import intercaat
import os
from scripts.pdb_fixinsert import fixInsert
from scripts.CLI import download_pdb
import shutil

def main():
    data = []
    proteins = {
            "6waq": ["A","B"],
            "1i8l": ["B","D"],
            "1a3r":["H","P"],
            "1a5h":["A","C"],
            "1a81":["A","B"],
            "1ab9":["B","A"],
            "1axd":["A","C"],
            "1aya":["A","P"],
            "1b8h":["B","D"]
            # "1joj":["B","Q"],
            # "1jp5":["A","C"],
            # "3mhp":["A","C"],
            # "2l7u":["A","B"],
            # "3dxc":["A","B"],
            # "3rtx":["A","C"],
            # "1i3z":["A","B"],
            # "2vif":["A","P"],
            # "148l":["E","S"],
            # "1b0r":["A","C"],
            # "1bqh":["A","C"],
            # "1f4v":["B","E"]
            }
    # radii = [2.4,3.4,4.4,5.5]
    radii = [5.4]
    mi_thresholds = [3,2,1]
    for pdb in proteins:
        pdb, pdb_file = preprocess(pdb)
        try:
            wt_interface = intercaatRun(pdb_file, proteins[pdb][0], proteins[pdb][1])
        except:
            print("error")
            continue
        for sr in radii:
            for mi in mi_thresholds:
                extended_interface = interfaceExtention(pdb,pdb_file,proteins[pdb][0], proteins[pdb][1], sr, mi,result_file = f"{pdb}_{sr}_{mi}.txt", intercaat_result = wt_interface)
                print(pdb, pdb_file, sr, mi, extended_interface)
                extended_interface = " ".join(extended_interface)
                data.append([pdb,str(sr),str(mi),extended_interface])
        shutil.make_archive(f"output/{pdb}/mutants_{pdb}" , 'zip', f"output/{pdb}/mutants")
        shutil.rmtree(f"output/{pdb}/mutants")
    with open(f"benchmark_data.csv","a+") as bmark_fh:
        for line in data:
            bmark_fh.write(",".join(line) + "\n")
    print(data)
    return

def intercaatRun(in_pdb, query_chain, partner_chain, qhull =False):
    intercaat_result = intercaat(in_pdb, query_chain, partner_chain, fp="input/", qhull = qhull)
    return intercaat_result

def preprocess(pdb):
    files = [i for i in os.listdir("input/") if i.endswith(".pdb")]
    pdb, pdb_file = pdbManager(pdb, files)
    pdb_file = fixInsert(pdb_file)
    print(pdb, pdb_file)
    os.makedirs("output/", exist_ok=True)
    os.makedirs(f"output/{pdb}/mutants/", exist_ok=True)
    return pdb, pdb_file

def pdbManager(pdb, files):
    pdb_file = pdb + ".pdb"
    if pdb_file not in files:
        download_pdb(pdb, "input/")
    return pdb, pdb_file


# main()

def cleanup():
    data = []
    proteins = ["6waq","1i8l","1i85","1a3r","1a5h","1a81","1ab9","1axd","1aya","1b8h","1bq9","1joj","1jp5"]

    for folder in os.listdir("output/"):
        if folder in proteins:
            for filename in os.listdir(f"output/{folder}"):
                if filename.startswith(folder):
                    sr = filename.split("_")[1]
                    mi = filename.split("_")[2]
                    mi = mi.split(".")[0]
                    with open(f"output/{folder}/{filename}", "r") as fh:
                        for line in fh.readlines():
                            if line.startswith("Extened Interface Position(s): ") or line.startswith("Extended Interface Position(s):"):
                                line = line.rstrip()
                                positions = line.split(":")[1]
                                data.append([folder,sr, mi, positions])
    print(data)
    df = pd.DataFrame(data, columns = ["pdb","sr","mi","extended_interface"])
    print(df)
    df.to_csv("benchmark_data.csv")
    return

# cleanup()


def plot():
    df: pd.DataFrame = pd.read_csv("benchmark_data.csv")  # type: ignore
    df["extention_length"] = [len(i.split()) if isinstance(i, str) else 0 for i in df["extended_interface"]]
    print(df)
    for pdb, sdf in df.groupby("pdb"):
        fig, axes = plt.subplots(1,4, sharey = True)
        for c, (sr, ssdf) in enumerate(sdf.groupby("sr")):
            ssdf.plot(kind = "scatter", x = "mi", y = "extention_length",title = f"{pdb} | sr: {sr}",  ax = axes[c])

        # for c, (mi, ssdf) in enumerate(sdf.groupby("mi")):
            # ssdf.plot(kind = "scatter", x = "sr", y = "extention_length",title = f"{pdb} | mi: {mi}",  ax = axes[1,c])
        plt.tight_layout()
        plt.savefig(f"output/benchmark/plots/{pdb}_benchmark.png")
    return

# plot()


def findMinParams():
    data = {}
    hist = {
            (2.4,1):0,
            (2.4,2):0,
            (2.4,3):0,
            (3.4,1):0,
            (3.4,2):0,
            (3.4,3):0,
            (4.4,1):0,
            (4.4,2):0,
            (4.4,3):0,
            (5.4,1):0,
            (5.4,2):0,
            (5.4,3):0}

    df: pd.DataFrame = pd.read_csv("benchmark_data.csv")  # type: ignore
    df["extention_length"] = [len(i.split()) if isinstance(i, str) else 0 for i in df["extended_interface"]]
    for pos, sdf in df.groupby("pdb"):
        ssdf = sdf[sdf["extention_length"] == sdf["extention_length"].max()]
        mi_sdf = ssdf[ssdf["mi"] == ssdf["mi"].max()]
        sr_sdf = mi_sdf[ssdf["sr"] == mi_sdf["sr"].min()]
        value  =  (float(sr_sdf["sr"]), float(sr_sdf["mi"]))
        data[pos] = value
        hist[value]  = hist[value] + 1  # type: ignore

    x = [str(i) for i in hist]
    y = hist.values()
    plt.bar(x, y)
    plt.ylabel("frequency")
    plt.xlabel("(sr, mi)")
    plt.show()

    x = list(data.values())
    print(x)
    sr = [i[0] for i in x]
    mi = [i[1] for i in x]
    print(sr, mi)
    h = plt.hist2d(mi,sr, cmap = "Greys")
    plt.ylabel("SR")
    plt.xlabel("MI")
    plt.colorbar(h[3])
    plt.xticks(np.arange(min(mi), max(mi)+1, 1.0))
    plt.yticks(np.arange(min(sr), max(sr)+1, 1.0))
    plt.tight_layout()
    plt.show()
    return
findMinParams()
