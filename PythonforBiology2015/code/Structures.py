from Alignments import sequenceAlign, BLOSUM62
from Maths import calcTorsionAngle, getRotationMatrix
from Modelling import getStructuresFromFile, Structure, Chain, Residue, Atom

from math import pi, sqrt, degrees
from matplotlib import pyplot
from numpy import array, dot, zeros, ones, cross, sqrt, linalg, exp, identity
  
try:
  # Python 3
  from urllib.request import urlopen
except ImportError:
  # Python 2
  from urllib2 import urlopen

try:
  from Bio import PDB
  BIOPYTHON_INSTALLED = True
  
except ImportError:
  print('\n* * * * * BioPython not installed * * * * * ')
  BIOPYTHON_INSTALLED = False

try:
  import pymol
  PYMOL_INSTALLED = True
  
except ImportError:
  print('\n* * * * * PyMol not installed * * * * * ')
  PYMOL_INSTALLED = False

PDB_URL = 'http://www.rcsb.org/pdb/cgi/export.cgi/' \
          '%s.pdb?format=PDB&compression=None'

ATOMIC_NUMS = {'H':1, 'C':12, 'N':14, 'O':16, 'P':31, 'S':32}

THREE_LETTER_TO_ONE = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E',
                       'PHE':'F','GLY':'G','HIS':'H','ILE':'I',
                       'LYS':'K','LEU':'L','MET':'M','ASN':'N',
                       'PRO':'P','GLN':'Q','ARG':'R','SER':'S',
                       'THR':'T','VAL':'V','TRP':'W','TYR':'Y',
                       'G':'G','C':'C','A':'A','T':'T','U':'U'}



def downloadPDB(pdbId, fileName=None):

  if not fileName:
    fileName = '%s.pdb' % pdbId

  response = urlopen(PDB_URL % pdbId)
  data = response.read().decode('utf-8')

  fileObj = open(fileName, 'w')
  fileObj.write(data)
  fileObj.close()

  return fileName
  
  
def writeStructureToFile(structure, fileName):
  
  fileObj = open(fileName, 'w') # Should really check if the file exists and warn
  
  fileObj.write('HEADER    %s\n' % structure.pdbId)
  fileObj.write('TITLE     %s\n' % structure.name)
  pdbFormat = '%-6.6s%5.1d %4.4s %3.3s %s%4.1d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n'
  
  prefix = 'ATOM' 
  
  i = 0
  for chain in structure.chains:
    for residue in chain.residues:
      for atom in residue.atoms:
        i += 1
        x, y, z = atom.coords
        data = (prefix, i, '%-3.3s' % atom.name, residue.code,
                chain.code, residue.seqId, x, y, z, 0.0, 0.0,   # B-factor not used in our model
                atom.element)
        line = pdbFormat % data
        fileObj.write(line)

  fileObj.close()
  

def getCenterOfMass(structure):

  centerOfMass = zeros(3, float)
  totalMass = 0.0

  for chain in structure.chains:
    for residue in chain.residues:
      for atom in residue.atoms:
        mass = ATOMIC_NUMS.get(atom.element, 12.0)
        centerOfMass += mass * atom.coords
        totalMass += mass
  
  centerOfMass /= totalMass

  return centerOfMass


def rotateStructure(structure, axis=(1,0,0), angle=0):

  rMatrix = array(getRotationMatrix(axis, angle))
  
  for chain in structure.chains:
    for residue in chain.residues:
      for atom in residue.atoms:
        newCoords = dot(rMatrix, atom.coords)      
        atom.coords = newCoords
        
        
def getAtomCoords(structure):

  coords = []
  atoms = []
  
  for chain in structure.chains:
    for residue in chain.residues:
      for atom in residue.atoms:
        coords.append(atom.coords)
        atoms.append(atom)

  return atoms, array(coords) 
  
  
def rotateStructureNumPy(structure, axis=array([1,0,0]), angle=0):

  rMatrix = array(getRotationMatrix(axis, angle))

  atoms, coords = getAtomCoords(structure)
  
  coords = dot(coords, rMatrix.T)
  
  for index, atom in enumerate(atoms):
    atom.coords = list(coords[index])


def affineTransformStructure(structure, transform=identity(3),
                             translate=(0,0,0)):

  atoms, coords = getAtomCoords(structure)

  coords = dot(coords, transform.T)
  
  coords = coords + translate

  for index, atom in enumerate(atoms):
    atom.coords = coords[index]


def findCloseAtoms(structure, xyz, limit=5.0):

  closeAtoms = []
  xyz = array(xyz)
  atoms, coords = getAtomCoords(structure)
  limit2 = limit * limit
  
  deltas = coords - xyz
  squares = deltas * deltas
  sumSquares = squares.sum(axis=1)  
  
  boolArray = sumSquares < limit2
  indices = boolArray.nonzero()[0]
  
  closeAtoms = [atoms[i] for i in indices]

  return closeAtoms


def getPhiPsi(residue, inDegrees=True):

  phi = None
  psi = None

  chain = residue.chain
  residues = chain.residues

  atomN  = residue.getAtom('N')
  atomCa = residue.getAtom('CA')
  atomC  = residue.getAtom('C')
  
  coordsN  = atomN.coords
  coordsCa = atomCa.coords
  coordsC  = atomC.coords

  index = residues.index(residue)
  
  if index > 0:

    residuePrev = residues[index-1]
    
    atomC0 = residuePrev.getAtom('C')
    coordsC0  = atomC0.coords
    
    phi = calcTorsionAngle(coordsC0, coordsN, coordsCa, coordsC)
    
    if inDegrees:
      phi = degrees(phi)
    
  if index < ( len(residues)-1 ):
    residueNext = residues[index+1]
    
    atomN2 = residueNext.getAtom('N')
    coordsN2  = atomN2.coords
  
    psi = calcTorsionAngle(coordsN, coordsCa, coordsC, coordsN2)
    
    if inDegrees:
      psi = degrees(psi)
  
  return phi, psi


def filterSubStructure(structure, chainCodes=None,
                       residueIds=None, atomNames=None):

  name = structure.name + '_filter'
  
  if chainCodes:
    name += ' ' + ','.join(chainCodes)
    chainCodes = set(chainCodes)
  
  if residueIds:
    name += ' ' + ','.join([str(x) for x in residueIds])
    residueIds = set(residueIds)

  if atomNames:
    name += ' ' + ','.join(atomNames)
    atomNames = set(atomNames)

  conf = structure.conformation
  pdbId = structure.pdbId

  filterStruc = Structure(name=name, conformation=conf, pdbId=pdbId)

  for chain in structure.chains:
    if chainCodes and (chain.code not in chainCodes):
      continue
  
    includeResidues = []
    
    for residue in chain.residues:
      if residueIds and (residue.seqId not in residueIds):
        continue

      includeAtoms = []

      for atom in residue.atoms:
        if atomNames and (atom.name not in atomNames):
          continue

        includeAtoms.append(atom)

      if includeAtoms:
        includeResidues.append( (residue, includeAtoms) )

    if includeResidues:
      filterChain = Chain(filterStruc, chain.code, chain.molType)
      
      for residue, atoms in includeResidues:
        filterResidue = Residue(filterChain, residue.seqId,
                                residue.code)
    
        for atom in atoms:
          coords = array(atom.coords)
          Atom(filterResidue, atom.name, coords, atom.element)
  
  return filterStruc


def copyStructure(structure):

  return filterSubStructure(structure, None, None, None)


def centerCoords(coords, weights):

  wCoords = coords.transpose() * weights

  xyzTotals = wCoords.sum(axis=1)

  center = xyzTotals/sum(weights)

  coords -= center
  
  return coords, center  


def alignCoords(coordsA, coordsB, weights=None):

  n = len(coordsA)
  if weights is None:
    weights = ones(n)

  rMat = dot(coordsB.transpose()*weights, coordsA)

  rMat1, scales, rMat2 = linalg.svd(rMat)

  sign = linalg.det(rMat1) * linalg.det(rMat2)

  if sign < 0:
    rMat1[:,2] *= -1
  
  rotation = dot(rMat1, rMat2)
    
  coordsB = dot(coordsB, rotation)
  
  return rotation, coordsB


def calcRmsds(refCoords, allCoords, weights):

  rmsds = []
  totalWeight = sum(weights)
  totalSquares = zeros(refCoords.shape)
  
  for coords in allCoords:
    delta = coords-refCoords
    squares = delta * delta
    totalSquares += squares
    sumSquares = weights*squares.sum(axis=1)
    rmsds.append( sqrt(sum(sumSquares)/totalWeight) )
  
  nStruct = len(allCoords)
  atomRmsds = sqrt(totalSquares.sum(axis=1)/nStruct)

  return rmsds, atomRmsds
  
  
def superimposeCoords(allCoords, weights, threshold=5.0):

  nStruct = len(allCoords)

  refCoords = allCoords[0]
  meanCoords = zeros(refCoords.shape)
  rotations = []

  for index, coords in enumerate(allCoords):
    if index == 0:
      rotation = identity(3)
    else:  
      rotation, coords = alignCoords(refCoords, coords, weights)
      allCoords[index] = coords # Update to aligned
      
    rotations.append(rotation)
    meanCoords += coords
    
  meanCoords /= nStruct
  
  rmsds, atomRmsds = calcRmsds(meanCoords, allCoords, weights)  
  
  bestRmsd = min(rmsds)
  bestIndex = rmsds.index(bestRmsd)
  bestCoords = allCoords[bestIndex]

  weightScale = atomRmsds/threshold
  weights *= exp(-weightScale*weightScale)
  
  meanCoords = bestCoords.copy()
  
  for index, coords in enumerate(allCoords):
    if index != bestIndex:
      rotation, coords = alignCoords(bestCoords, coords, weights)
      rotations[index] = rotation
      allCoords[index] = coords # Update to aligned
      meanCoords += coords
      
  meanCoords /= nStruct     
  weights = ones(len(weights))
  rmsds, atomRmsds = calcRmsds(meanCoords, allCoords, weights)
  
  return allCoords, rmsds, atomRmsds, rotations


def superimposeStructures(structures):

  weights = None
  allCoords = []

  for structure in structures:
    atoms, coords = getAtomCoords(structure)
    
    if weights is None:
      weights = [ATOMIC_NUMS[atom.element] for atom in atoms]
      weights = array(weights)

    coords, center = centerCoords(coords, weights)
    allCoords.append(coords)
  
  results = superimposeCoords(allCoords, weights)
  allCoords, rmsds, atomRmsds, rotations = results  
  
  for i, structure in enumerate(structures):
    atoms, oldCoords = getAtomCoords(structure)
    
    for j, atom in enumerate(atoms):
      atom.coords = allCoords[i][j]

  return rmsds, atomRmsds



def getChainSequence(chain):

  letters = []
    
  for residue in chain.residues:
    code =  residue.code
    letter = THREE_LETTER_TO_ONE.get(code, 'X')
    letters.append(letter)
    
  seq = ''.join(letters)

  return seq


def seqStructureBackboneAlign(chainA, chainB, 
                              atomNames=set(['CA','C','N']),
                              simMatrix=BLOSUM62):

  structureA = chainA.structure
  structureB = chainB.structure
  residuesA = chainA.residues
  residuesB = chainB.residues

  seqA = getChainSequence(chainA)
  seqB = getChainSequence(chainB)

  score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix)

  pairedPosA = []
  pairedPosB = []  
  posA = 0
  posB = 0

  for i in range(len(alignA)):
    # No dashes in both at same location
  
    if alignA[i] == '-':
      posB += 1
   
    elif alignB[i] == '-':
      posA += 1
    
    else:
      pairedPosA.append(posA)
      pairedPosB.append(posB)
      
      posA += 1
      posB += 1

  filterIdsA = [residuesA[p].seqId for p in pairedPosA]
  filterIdsB = [residuesB[p].seqId for p in pairedPosB]

  backboneStrucA = filterSubStructure(structureA, [chainA.code], 
                                      filterIdsA, atomNames)
  backboneStrucB = filterSubStructure(structureB, [chainB.code], 
                                      filterIdsB, atomNames)
  
  atomsA, coordsA = getAtomCoords(backboneStrucA)
  atomsB, coordsB = getAtomCoords(backboneStrucB)
  weights = ones(len(atomsA))
  
  coordsA, centerA = centerCoords(coordsA, weights)
  coordsB, centerB = centerCoords(coordsB, weights)
    
  coords = [coordsA, coordsB]
  c, rmsds, atomRmsds, rotations = superimposeCoords(coords, weights)

  affineTransformStructure(structureA, rotations[0], -centerA)  
  affineTransformStructure(structureB, rotations[1], -centerB)  
    
  return rmsds, atomRmsds

  
def findCloseAtomsBioPy(structure, xyz, conf=0, limit=5.0):

  if not isinstance(structure, PDB.Structure.Structure):
    raise Exception('Structure must be Bio.PDB.Structure class')

  closeAtoms = []
  xyz = array(xyz)
  limit2 = limit * limit
  
  coords = []
  atoms = []
  confModel = structure[conf]
 
  for chain in confModel:
    for residue in chain:
      for atom in residue:
        coords.append(atom.coord)
        atoms.append(atom)

  deltas = coords - xyz
  squares = deltas * deltas
  sumSquares = squares.sum(axis=1)  
  
  boolArray = sumSquares < limit2
  indices = boolArray.nonzero()[0]
  
  closeAtoms = [atoms[i] for i in indices]
   
  return closeAtoms


if __name__ == '__main__':

  print('\nGet a PDB file')

  fileName = 'examples/1A12.pdb'

  print('\nMake structure object')
  strucObjs = getStructuresFromFile(fileName)

  print('\nGet centre of mass')
  struc = getStructuresFromFile(fileName)[0]
  print(getCenterOfMass(struc))

  print('\nStructure tranformations')
  rotateStructure(struc, (0,1,0), pi/2)

  rotateStructureNumPy(struc, (0,1,0), pi/2)

  mirrorTransform = array([[-1,0,0], [0,-1,0], [0,0,-1]])

  translate = array([10.0,10.0,0.0])

  affineTransformStructure(struc, mirrorTransform, translate)

  writeStructureToFile(struc, 'testTransform.pdb')

  print('\nCalculate a distance')
  residues = struc.getChain('A'). residues
  atomA = residues[0].getAtom('CA')
  atomB = residues[1].getAtom('CA')

  deltas = atomA.coords - atomB.coords
  sumSquares = dot(deltas, deltas)
  distance = sqrt( sumSquares )

  print(distance)


  print('\nFind close atoms')

  atoms = findCloseAtoms(struc, (18.89, 0.36, 1.24), 5.0)
  for atom in atoms:
    print(atom.residue.code, atom.name)


  print('\nPlot Phi and Psi backbone angles')

  phiList = []
  psiList = []
  for chain in struc.chains:
    for residue in chain.residues[1:-1]:
      phi, psi = getPhiPsi(residue)
      phiList.append(phi)
      psiList.append(psi)

  pyplot.scatter(phiList, psiList)
  pyplot.axis([-180,180,-180,180])
  pyplot.show()


  print('\nGet backbone sub-structure')

  chainCodes = set(['A'])
  residueIds = None                # No residue filter: all of them
  atomNames  = set(['N','CA','C']) # Heavy backbone atoms (not H)

  chain_A_backbone = filterSubStructure(struc, chainCodes,
                                        residueIds, atomNames)


  print('\nTest same structure RMSDs')

  strucA = getStructuresFromFile(fileName)[0]
  strucB = getStructuresFromFile(fileName)[0]

  rotateStructureNumPy(strucA, (1,0,0), pi/3)

  rmsds, atomRmsds  = superimposeStructures([strucA, strucB])

  print(rmsds)  # Hopefully all close to zero


  print('\nRMSDs for bundle of conformational models')

  fileName = downloadPDB('1UST')
  strucObjs = getStructuresFromFile(fileName)
  rmsds, atomRmsds = superimposeStructures(strucObjs)

  print(rmsds)


  print('\nHomologue alignment')
  struc1 = getStructuresFromFile(downloadPDB('1UST'))[0]
  struc2 = getStructuresFromFile(downloadPDB('1HST'))[0]

  chain1 = struc1.getChain('A')
  chain2 = struc2.getChain('A')

  rmsds, atomRmsds = seqStructureBackboneAlign(chain1, chain2)

  print(rmsds)


  if BIOPYTHON_INSTALLED:
    print('\nBioPython PDB handling')
 
    fileName = 'examples/1UST.pdb'
    parser = PDB.PDBParser()
    struc  = parser.get_structure('Name', fileName)

    conformation = struc[0]

    for chain in conformation:
      for residue in chain:
        atomNames = [a.name for a in residue]
        print(chain.id, residue.id[1], residue.resname, atomNames)
        caAtom = residue['CA']
        print(caAtom.name, caAtom.coord, caAtom.bfactor)
 
    outFileName = 'test.pdb'
    writer = PDB.PDBIO()
    writer.set_structure(struc)
    writer.save(outFileName)

    print('\nBioPython find close atoms')

    atoms = findCloseAtomsBioPy(struc, (18.89, 0.36, 1.24), 5.0)
    for atom in atoms:
      print(atom.name, atom.coord)


  if PYMOL_INSTALLED:
    print('\nPyMol examples')
 
    #Maybe need something like the following:
    #import os
    #os.environ['PYMOL_PATH'] = os.path.dirname(pymol.__file__)

    pymol.finish_launching()

    fileName = downloadPDB('1AFO')
    strucName = 'Glycophorin'
    pymol.cmd.load(fileName, strucName)

    pymol.cmd.select('bb','name n+c+o+ca')
    pymol.cmd.select('tmd', 'resi 71-100 and name n+c+o+ca ')

    pymol.cmd.color('gray', strucName)
    pymol.cmd.color('red', 'bb')
    pymol.cmd.show('cartoon', 'bb')
    pymol.cmd.color('blue', 'tmd')
    pymol.cmd.hide('lines','hydro')

    outFileName = strucName + '.pdb'
    pymol.cmd.save(outFileName, strucName, 0, 'pdb')

    pymol.cmd.png(strucName+'.png')

    pymol.cmd.quit()



