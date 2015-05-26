def formatStatement(text, placeHolder):

  if placeHolder == '%s':
    return text

  numInserts = text.count('%s')

  return text % (numInserts*(placeHolder,))

def addStructureToDb(connection, structure, placeHolder='%s'):

  cursor = connection.cursor()

  pdbId = structure.pdbId
  conformation = structure.conformation

  # check if structure already in database

  stmt = "select * from structure where pdbId=%s and conformation=%s"
  stmt = formatStatement(stmt, placeHolder)

  cursor.execute(stmt, (pdbId, conformation))
  if cursor.fetchone():
    cursor.close()
    raise Exception('structure with (pdbId=%s, conformation=%s) already known' % (pdbId, conformation))
 
  # if remaining commands are successful, then commit
  # but if any exception thrown then rollback

  try:
    # insert structure into database

    stmt = "insert into structure (name, pdbId, conformation) values (%s, %s, %s);"
    stmt = formatStatement(stmt, placeHolder)

    cursor.execute(stmt, (structure.name, pdbId, conformation))
    structureId = cursor.lastrowid

    # insert chains into database

    for chain in structure.chains:
      molType = chain.molType
      code = chain.code

      # insert chain into database

      stmt = "insert into chain (structureId, molType, code) values (%s, %s, %s)"
      stmt = formatStatement(stmt, placeHolder)

      cursor.execute(stmt, (structureId, molType, code))
      chainId = cursor.lastrowid

      # insert residues into database

      for residue in chain.residues:
        seqId = residue.seqId

        # insert residue into database

        stmt = "insert into residue (chainId, seqId, code) values (%s, %s, %s)"
        stmt = formatStatement(stmt, placeHolder)

        cursor.execute(stmt, (chainId, seqId, residue.code))
        residueId = cursor.lastrowid

        # insert atoms into database

        for atom in residue.atoms:
          # insert atom into database

          (x, y, z) = atom.coords

          stmt = "insert into atom (residueId, name, x, y, z, element) values (%s, %s, %s, %s, %s, %s)"
          stmt = formatStatement(stmt, placeHolder)

          cursor.execute(stmt, (residueId, atom.name, x, y, z, atom.element))

    cursor.close()
    connection.commit()

  except Exception as e:
    cursor.close()
    try:
      connection.rollback()
    except:
      pass

    raise e # re-raise original exception

def addStructureToDb2(connection, structure):

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

    for chain in structure.chains:
      molType = chain.molType
      code = chain.code

      # insert chain into database

      stmt = "insert into chain (structureId, molType, code) values (%s, %s, %s)"

      cursor.execute(stmt, (structureId, molType, code))
      chainId = cursor.lastrowid

      # insert residues into database

      for residue in chain.residues:
        seqId = residue.seqId

        # insert residue into database

        stmt = "insert into residue (chainId, seqId, code) values (%s, %s, %s)"

        cursor.execute(stmt, (chainId, seqId, residue.code))
        residueId = cursor.lastrowid

        # insert atoms into database

        values = []
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

def getStructureFromDb(connection, pdbId, conformation=1, placeHolder='%s'):

  from Modelling import Structure, Chain, Residue, Atom

  cursor = connection.cursor()

  try:
    # get structure record from database

    stmt = "SELECT structure.name, chain.molType, chain.code, residue.seqId, residue.code, atom.name, atom.x, atom.y, atom.z, atom.element FROM structure, chain, residue, atom WHERE structure.pdbId=%s AND structure.conformation=%s AND structure.id=chain.structureId AND chain.id=residue.chainId AND residue.id=atom.residueId ORDER BY chain.id, residue.id, atom.id"
    stmt = formatStatement(stmt, placeHolder)

    cursor.execute(stmt, (pdbId, conformation))
    result = cursor.fetchall()
    if not result:
      raise Exception('structure with (pdbId=%s, conformation=%s) not known' % (pdbId, conformation))

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

def countChainMolType(connection, molType, placeHolder='%s'):

  cursor = connection.cursor()

  try:
    # get matching chain records from database

    stmt = "select * from chain where molType=%s"
    stmt = formatStatement(stmt, placeHolder)

    cursor.execute(stmt, (molType,))
    result = cursor.fetchall()

    count = len(result)

  finally:
    cursor.close()

  return count
