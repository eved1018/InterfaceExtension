from intercaatWrapper import intercaat, writeFile




# pdbs = ["","",""]

# for pdb in pdbs:
#     filename = ""
#     newMatch, newInteractionRes, newInteraction = intercaat(pdb, ligand, receptor, mi = mi  , sr = sr)
#     writeFile(filename, newMatch, newInteractionRes, newInteraction)
newMatch, newInteractionRes, newInteraction = intercaat("1cph.pdb", "A", "B")
writeFile("testout.txt",newMatch, newInteractionRes, newInteraction  )
