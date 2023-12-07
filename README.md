# PeptideMaker
Based on the given sequence, model it through PeptideBuilder, add hydrogen atoms and end-capping groups (ACE and NME), and then perform a short energy minimization in vacumm through Openmm<br>
usage: python PeptideMaker.py -s {sequence} -o {output_path} -remove (remove process files or not) <br>
-s: the sequence of peptide you want to build, str <br>
-o: specify a path for output (defult path is the current path), str <br>
-remove: remove the files generated during the modelling process (defult is True), bool
