import sys
import PeptideBuilder
import Bio
import os
import math

def _fsgn(a):
 if a>0: return 1
 if a<0: return -1
 return 0

def VECADD(c,a,b):
 c[0]=a[0]+b[0]
 c[1]=a[1]+b[1]
 c[2]=a[2]+b[2]

def VECSUB(c,a,b):
 c[0]=a[0]-b[0]
 c[1]=a[1]-b[1]
 c[2]=a[2]-b[2]

def VECLEN(a):     
 return (a[0]*a[0]+a[1]*a[1]+a[2]*a[2])**0.5

def VECOUT(c,a,b):
  c[0]=a[1]*b[2]-a[2]*b[1]
  c[1]=a[2]*b[0]-a[0]*b[2]
  c[2]=a[0]*b[1]-a[1]*b[0]
  
def VECELONG(c,a,v):
  c[0]=a[0]*(v)
  c[1]=a[1]*(v)
  c[2]=a[2]*(v)

def calxyz(c1, c2, c3,  b,  a,  d):
  ci=[0.,0.,0.]
  r=[0.,0.,0.]
  h=[0.,0.,0.]
  r21=[0.,0.,0.]
  r23=[0.,0.,0.]
  rb=[0.,0.,0.]
  rb21=[0.,0.,0.]
  rc=[0.,0.,0.]
  
  PI=3.141592654
  ca=a;
  cb=b;
  cd=d;
  ca=ca/180.*PI;
  cd=cd/180.*PI;
             
  VECSUB(r21, c1, c2);
  VECSUB(r23, c3, c2);
           
  VECOUT(r, r21, r23);

  VECOUT(h, r21, r);

  rv=VECLEN(r);
  hv=VECLEN(h);
  r21v=VECLEN(r21);

  VECELONG(r,r,1./rv); 
  VECELONG(h,h,1./hv);
  VECELONG(r21,r21,1./r21v);


  rbv=cb*math.sin(ca);

  VECELONG(rb21, r21, math.cos(PI-ca)*cb); 
           
  VECELONG(rc, h, math.cos(PI-math.fabs(cd)) );
  VECELONG(rb, r, math.sin(math.fabs(cd))*_fsgn(cd) );
  VECADD (rb, rb, rc);
  VECELONG (rb, rb, rbv);

  VECADD( ci, rb, rb21);
  VECADD( ci, ci, c1);

  return ci        

# 1.im build pdb
def peptide_builder(peps_list, save_path='./', no_add_ALA=False):
  for p in peps_list:
    if not no_add_ALA:
      _p = 'A' + p + 'A'
    else:
      _p = p
    structure = PeptideBuilder.initialize_res(_p[0])
    for aa in _p[1:]:
      PeptideBuilder.add_residue(structure, aa)
    # PeptideBuilder.add_terminal_OXT(structure)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save(os.path.join(save_path, "./%s.pdb"%p))
    out.save(os.path.join(save_path, "./%s2.pdb"%p))

# 2. read and write pdb
def _pdb_atom_string(atom):
  return _pdb_string(atom[0],atom[1],atom[2],atom[3],atom[4],atom[5],atom[6],atom[7],atom[-1])

def _pdb_string(aid, anm, rnm, cnm, rid, x, y, z, snm):
  if len(anm)<3: 
    anm1=anm+' '
  else:
    anm1=anm
  if len(rnm)<4: 
    rnm1=rnm+' '
  else:
    rnm1=rnm
  if aid>99999: aid=99999
  return "ATOM  %5d%5s %4s%1s"%(aid, anm1,rnm1,cnm)+\
        '%4d    %8.3f%8.3f%8.3f'%(rid,x,y,z)+\
        '%18s%3s'%(' ',snm)

class PDBdata:
  def __init__(self, file=None):
    self.atoms = []
    if file is not None:
      self.readPDB(file)

  def readPDB(self, file):
    with open(file, 'r') as f:
      lines = [t.strip().split() for t in f.readlines()]
    for line in lines:
      if line[0] == 'ATOM' or line[0] == 'HETATM':
        self.atoms.append(line)
    for idx in range(len(self.atoms)):
      self.atoms[idx][1] = int(self.atoms[idx][1])
      self.atoms[idx][5] = int(self.atoms[idx][5])
      self.atoms[idx][6] = float(self.atoms[idx][6])
      self.atoms[idx][7] = float(self.atoms[idx][7])
      self.atoms[idx][8] = float(self.atoms[idx][8])
      self.atoms[idx] = self.atoms[idx][1:]
  
  def find_atom_xyz_by_name(self, atomName):
    for idx in range(len(self.atoms)):
      z = []
      if self.atoms[idx][1] == atomName:
          z.append(self.atoms[idx][5])
          z.append(self.atoms[idx][6])
          z.append(self.atoms[idx][7])
          break
    return z
  def update(self, atoms):
    self.atoms = atoms

  def to_pdb(self, file):
    with open(file, 'w') as f:
      for a in self.atoms:
        f.write(_pdb_atom_string(a) + '\n')
      f.write('END\n')
  



# 3. add FMOC
FMOC_builder_PACE = [
  ['C12', 'CF1', 'OF1', 'C', 1.52, 109., 180., 'C'],
  ['C11', 'C12', 'CF1', 'OF1', 1.39, 109., 180., 'C'],
  ['C10', 'C11', 'C12', 'CF1', 1.39, 120., -60., 'C'],
  ['C9', 'C10', 'C11', 'C12', 1.39, 120., 180., 'C'],
  ['C8', 'C9', 'C10', 'C11', 1.39, 120., 0., 'C'],
  ['C7', 'C8', 'C9', 'C10', 1.39, 120., 0., 'C'],
  ['C6', 'C7', 'C8', 'C9', 1.39, 120., 0., 'C'],
  ['C5', 'C6', 'C7', 'C8', 1.39, 130., 180., 'C'],
  ['C4', 'C5', 'C6', 'C7', 1.39, 130., 0., 'C'],
  ['C3', 'C4', 'C5', 'C6', 1.39, 120., 180., 'C'],
  ['C2', 'C3', 'C4', 'C5', 1.39, 120., 0., 'C'],
  ['C1', 'C2', 'C3', 'C4', 1.39, 120., 0., 'C'],
  ['C13', 'C1', 'C2', 'C3', 1.39, 120., 0., 'C']
  ]

def find_atom(atoms, name, aa_index=1):
  for i, atom in enumerate(atoms):
    if atom[1] == name and atom[4] == aa_index:
      return i

def add_ACE(file, save_path):
  oriPDB = PDBdata(file)
  # del ALA1-CB
  oriAtoms = oriPDB.atoms[1:4] + oriPDB.atoms[5:-3]
  # change ALA1-CA to OM, ALA1-N to C14
  oriAtoms[-1][2] = 'NME'
  oriAtoms[-2][2] = 'NME'
  oriAtoms[-1][1] = 'CH3'
  oriAtoms[0][1] = 'CH3'
  # update index
  for i in range(len(oriAtoms)):
    oriAtoms[i][0] = i
  # update aa type
  for i, atom in enumerate(oriAtoms):
    if atom[4] == 1:
      oriAtoms[i][2] = 'ACE'
    else:
      break
  oriPDB.update(oriAtoms)
  oriPDB.to_pdb(save_path)
  return oriPDB.atoms

def add_Fmoc(file, save_path, builder):
  oriPDB = PDBdata(file)
  # del ALA1-CB
  oriAtoms = oriPDB.atoms[0:4] + oriPDB.atoms[5:-3]
  # change ALA1-CA to OM, ALA1-N to C14
  oriAtoms[-1][2] = 'NME'
  oriAtoms[-1][1] = 'CH3'
  oriAtoms[-2][2] = 'NME'
  # update index
  '''
  for i in range(len(oriAtoms)):
    oriAtoms[i][0] = i
  '''

  # change ALA1-CA to OF1, ALA1-N to CF1
  oriAtoms[0][1] = 'CF1'
  oriAtoms[0][-1] = 'C'
  oriAtoms[1][1] = 'OF1'
  oriAtoms[1][-1] = 'O'
  oriAtoms[3][1] = 'OF2'

  oriPDB.update(oriAtoms)
  for build in builder:
    anm, anm1, anm2, anm3, distance, angel, dihedral, snm = build
    c1 = oriPDB.find_atom_xyz_by_name(anm1)
    c2 = oriPDB.find_atom_xyz_by_name(anm2)
    c3 = oriPDB.find_atom_xyz_by_name(anm3)
   # print(c1,c2,c3)
    ci = calxyz(c1, c2, c3, distance, angel, dihedral)
    x = ci [0]
    y = ci [1]
    z = ci [2]
    atom = [0, anm, 'FMO', 'A', 1, x, y ,z, snm]
    oriAtoms.insert(0, atom)
    oriPDB.update(oriAtoms)
    # update aa type
  for i, atom in enumerate(oriAtoms):
    if atom[4] == 1:
      oriAtoms[i][2] = 'FMO'
    else:
      break
  for i in range(len(oriAtoms)):
    oriAtoms[i][0] = i
  oriPDB.update(oriAtoms)
  oriPDB.to_pdb(save_path)
  return oriPDB.atoms



if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('-p', '--save_path', type=str, default='./')
  parser.add_argument('-s', '--sequence', type=str, required=True)
  parser.add_argument('--no_add_ALA', action='store_true')
  args = parser.parse_args()
  
  os.makedirs(args.save_path, exist_ok=True)
  # build peptide PDB
  peps_list = args.sequence.split('|')

  peptide_builder(peps_list, args.save_path, args.no_add_ALA)
  for p in peps_list:
    add_ACE(os.path.join(args.save_path, '%s.pdb'%p),os.path.join(args.save_path, "ACE_%s.pdb"%p))
    add_Fmoc(os.path.join(args.save_path, '%s2.pdb'%p),os.path.join(args.save_path, "Fmoc_%s.pdb"%p), FMOC_builder_PACE)
