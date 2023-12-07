# PeptideMaker
Based on the given sequence, model it through PeptideBuilder, add hydrogen atoms and end-capping groups (ACE and NME), and then perform a short energy minimization in vacumm through Openmm<br>
usage: python PeptideMaker.py -s &lt;str&gt; -o &lt;str&gt; &#91;-helix&#93; &#91;-retain&#93; <br>
-s: the sequence of peptide you want to build <br>
-o: specify a path for output (if the path do not exsit, it will be created) <br>
-helix: optional flag, add this if you want the peptide be modelled into helix (otherwise it will be loop) <br>
-retain: optional flag, add this if you want to retain all the files generated during the modelling process (otherwise they will be deleted)
