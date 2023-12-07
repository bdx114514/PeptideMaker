from Capped import add_ACE
import os
import Bio.PDB
import PeptideBuilder
from PeptideBuilder import Geometry
from pdbfixer import PDBFixer
from simtk import openmm, unit
from openmm import app
from openmm.app import PDBFile

def build_pep(seq,out_dir,remove_process):
    #1.build peptide backbone
    os.mkdir(out_dir)
    seq = "A" + seq + "A"
    reslist = list(seq)
    for i in range(len(reslist)):
        geo = Geometry.geometry(reslist[i])
        geo.phi = -60
        geo.psi_im1 = -40
        if i == 0:
            structure = PeptideBuilder.initialize_res(geo)
        else:
            PeptideBuilder.add_residue(structure, geo)
    # PeptideBuilder.add_terminal_OXT(structure)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    save_path = f"{out_dir}/step1.pdb"
    out.save(save_path)
    
    #2.add end-capping groups and H
    add_ACE(save_path,f"{out_dir}/step2.pdb")
    fixer = PDBFixer(filename = f"{out_dir}/step2.pdb")
    fixer.addMissingHydrogens(7.0)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(f"{out_dir}/step3.pdb", 'w'))

    #3.short EM
    forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')
    pdb = app.PDBFile(f"{out_dir}/step3.pdb") 
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff)

    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 2 * unit.femtoseconds)
    simulation = app.Simulation(pdb.topology, system, integrator)

    simulation.context.setPositions(pdb.positions)

    print("Minimizing energy...")
    simulation.minimizeEnergy()

    minimized_state = simulation.context.getState(getEnergy=True, getPositions=True)
    minimized_energy = minimized_state.getPotentialEnergy()
    minimized_positions = minimized_state.getPositions()

    with open(f"{out_dir}/ACE-{seq}-NME.pdb", 'w') as pdb_file:
        app.PDBFile.writeFile(simulation.topology, minimized_positions, pdb_file)

    print(f"Minimized energy: {minimized_energy}")
    print(f"Minimized structure saved to {out_dir}/ACE-{seq}-NME.pdb")
    if remove_process == True:
        os.system(f"rm {out_dir}/step*")

  
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--seq', type=str, required=True)
    parser.add_argument('-o', '--out_dir', type=str, default='./')    
    parser.add_argument('-remove', '--remove_process', action='store_true', default=True)
    args = parser.parse_args()
    build_pep(args.seq,args.out_dir,args.remove_process)
    os.system

