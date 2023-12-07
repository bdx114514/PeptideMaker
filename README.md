# PeptideMaker
Based on the given sequence, model it through PeptideBuilder, add hydrogen atoms and end-capping groups (ACE and NME), and then perform a short energy minimization in vacumm through Openmm<br>
usage: python PeptideMaker.py -s &lt;str&gt; -o &lt;str&gt; -remove &lt;bool&gt; <br>
-s: the sequence of peptide you want to build <br>
-o: specify a path for output (if the path do not exsit, it will be created) <br>
-remove: remove the files generated during the modelling process (defult is True)
