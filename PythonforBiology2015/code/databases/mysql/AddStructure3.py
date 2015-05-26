def addStructureToDb3(connection, structure):

  cursor = connection.cursor()

  pdbId = structure.pdbId
  conformation = structure.conformation

  # check if structure already in database

  stmt = "select * from structure where pdbId=%s and conformation=%s"

  cursor.execute(stmt, (pdbId, conformation))
  if cursor.fetchone():
    cursor.close()
    raise Exception('structure with (pdbId=%s, conformation=%s) already known' % (pdbId, conformation))
 
  # if remaining commands are successful, then commit
  # but if any exception thrown then rollback

  try:
    # insert structure into database

    stmt = "insert into structure (name, pdbId, conformation) values (%s, %s, %s);"

    cursor.execute(stmt, (structure.name, pdbId, conformation))
    structureId = cursor.lastrowid

    # insert chains into database
    values = []

    for chain in structure.chains:
      molType = chain.molType
      code = chain.code
      values.extend([structureId, molType, code])
    n = len(values) / 3
    ss = '(%s, %s, %s)'
    ss = n * [ss]
    ss = ','.join(ss)
    stmt = "insert into chain (structureId, molType, code) values " + ss
    cursor.execute(stmt, values)

    stmt = "select id, molType, code from chain where structureId=%s"
    cursor.execute(stmt, (structureId,))
    result = cursor.fetchall()
    chainDict = {}
    for chainId, molType, code in result:
      chainDict[(molType, code)] = chainId

    values = []
    for chain in structure.chains:
      molType = chain.molType
      code = chain.code
      chainId = chainDict[(molType, code)]

      # insert residues into database

      for residue in chain.residues:
        seqId = residue.seqId
        values.extend([chainId, seqId, residue.code])
    n = len(values) / 3
    ss = '(%s, %s, %s)'
    ss = n * [ss]
    ss = ','.join(ss)
    stmt = "insert into residue (chainId, seqId, code) values " + ss
    cursor.execute(stmt, values)

    stmt = "select chain.id, residue.id, residue.seqId from chain, residue where chain.structureId=%s and chain.id=residue.chainId"
    cursor.execute(stmt, (structureId,))
    result = cursor.fetchall()
    residueDict = {}
    for chainId, residueId, seqId in result:
      residueDict[(chainId, seqId)] = residueId

    values = []
    for chain in structure.chains:
      molType = chain.molType
      code = chain.code
      chainId = chainDict[(molType, code)]
      
      for residue in chain.residues:
        seqId = residue.seqId
        residueId = residueDict[(chainId, seqId)]
        # insert atoms into database

        for atom in residue.atoms:
          # insert atom into database

          (x, y, z) = atom.coords

          values.extend([residueId, atom.name, x, y, z, atom.element])

    n = len(values) / 6
    ss = '(%s, %s, %s, %s, %s, %s)'
    ss = n * [ss]
    ss = ','.join(ss)
    stmt = "insert into atom (residueId, name, x, y, z, element) values " + ss
    cursor.execute(stmt, values)

    cursor.close()
    connection.commit()

  except Exception as e:
    cursor.close()
    try:
      connection.rollback()
    except:
      pass

    raise e # re-raise original exception

if __name__ == '__main__':

  import sys
  import time

  if len(sys.argv) != 3:
    print('need to specify database and PDB file')
    sys.exit(1)

  database = sys.argv[1]
  pdbFile = sys.argv[2]

  from Modelling import getStructuresFromFile
  structures = getStructuresFromFile(pdbFile)

  import MySQLdb
  import getpass
  user = getpass.getuser()
  pwd = getpass.getpass()
  connection = MySQLdb.connect(db=database, user=user, passwd=pwd)

  t0 = time.time()

  try:
    for structure in structures:
      addStructureToDb3(connection, structure)
  finally:
    connection.close()

  t1 = time.time()
  print('time taken = %.1f' % (t1-t0))
