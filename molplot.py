from openeye import oechem
mol = oechem.OEGraphMol()
oechem.OESmilesToMol(mol, "c1cc(N)cc(S(=O)(=O)O)c1 3-aminobenzenesulfonic acid")
oedepict.OEPrepareDepiction(mol)
oedepict.OERenderMolecule("DepictMolSimple.png", mol)