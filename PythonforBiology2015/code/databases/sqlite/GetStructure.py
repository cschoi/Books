from Modelling import Structure, Chain, Residue, Atom

def formatStatement(text, placeHolder):

  if placeHolder == '%s':
    return text

  numInserts = text.count('%s')

  return text % (numInserts*(placeHolder,))

times = []
def getStructureFromDb(connection, pdbId, conformation=1, placeHolder='%s'):

  cursor = connection.cursor()

  try:
    # get structure record from database

    t0 = time.time()
    stmt = "SELECT structure.name, chain.molType, chain.code, residue.seqId, residue.code, atom.name, atom.x, atom.y, atom.z, atom.element FROM structure, chain, residue, atom WHERE structure.pdbId=%s AND structure.conformation=%s AND structure.id=chain.structureId AND chain.id=residue.chainId AND residue.id=atom.residueId ORDER BY chain.id, residue.id, atom.id"
    stmt = formatStatement(stmt, placeHolder)

    cursor.execute(stmt, (pdbId, conformation))
    result = cursor.fetchall()
    if not result:
      raise Exception('structure with (pdbId=%s, conformation=%s) not known' % (pdbId, conformation))
    t1 = time.time()
    times.append(t1-t0)

    structure = chain = residue = atom = None
    for (structureName, chainMolType, chainCode, residueSeqId, residueCode, atomName, atomX, atomY, atomZ, atomElement) in result:

      # create structure object if needed
      if not structure:
        structure = Structure(structureName, conformation, pdbId)

      # create chain object if needed
      if not chain or chain.code != chainCode:
        chain = Chain(structure, chainCode, chainMolType)

      # create residue object if needed
      if not residue or residue.chain != chain or residue.seqId != residueSeqId:
        residue = Residue(chain, residueSeqId, residueCode)

      # create atom object
      coords = (atomX, atomY, atomZ)
      Atom(residue, atomName, coords, atomElement)

  finally:
    cursor.close()

  return structure

if __name__ == '__main__':

  import sys
  import time

  if len(sys.argv) not in (3, 4):
    print('need to specify database, PDB id, [conformation=1]')
    sys.exit(1)

  database = sys.argv[1]
  pdbId = sys.argv[2]

  if len(sys.argv) == 3:
    conformation = 1
  else:
    conformation = int(sys.argv[3])

  import sqlite3
  connection = sqlite3.connect(database)
  placeHolder = '?'

  t0 = time.time()

  try:
    for conformation in range(1,21):
      structure = getStructureFromDb(connection, pdbId, conformation, placeHolder)
  finally:
    connection.close()

  t1 = time.time()
  print('time taken = %.3f' % (t1-t0))
  print('sql time = %.3f' % sum(times))

