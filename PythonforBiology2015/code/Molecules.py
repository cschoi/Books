class Molecule:

  def __init__(self, name):
    if not name:
      raise Exception('name must be set to something')
    self.name = name

  def getName(self):
    return self.name

  def getCapitalisedName(self):
    name = self.getName()
    return name.capitalize()

class Protein(Molecule):

  def __init__(self, name, sequence):
    Molecule.__init__(self, name)
    self.aminoAcids = []

    for code in sequence:
      aminoAcid = AminoAcid(code)
      self.aminoAcids.append(aminoAcid)

  def getAminoAcids(self):
    return self.aminoAcids

  def getSequence(self):
    return [aminoAcid.code for aminoAcid in self.aminoAcids]

  def getMass(self):
    mass = 18.02  # N-terminus H and C-terminus OH
    aminoAcids = self.getAminoAcids()
    for aminoAcid in aminoAcids:
      mass += aminoAcid.getMass()
    return mass

class AminoAcid:

  massDict = { "A": 71.07, "R":156.18, "N":114.08, "D":115.08,
               "C":103.10, "Q":128.13, "E":129.11, "G": 57.05,
               "H":137.14, "I":113.15, "L":113.15, "K":128.17,
               "M":131.19, "F":147.17, "P": 97.11, "S": 87.07,
               "T":101.10, "W":186.20, "Y":163.17, "V": 99.13 }

  acceptableCodes = set(massDict.keys())

  def __init__(self, code):
    if code not in self.acceptableCodes:
      text = 'code = "%s", must be in list %s'
      raise Exception(text % (code, sorted(self.acceptableCodes)))
    self.code = code

  def getMass(self):
    return self.massDict[self.code]

if __name__ == '__main__':

  molecule = Molecule('moleculeName')
  print('molecule attributes')
  print('molecule name =', molecule.name)

  print('molecule function calls')
  print('molecule name =', molecule.getName())
  print('molecule capitalisedName =', molecule.getCapitalisedName())

  protein = Protein('proteinName', 'ACGCATH')
  print('protein attributes')
  print('protein name =', protein.name)
  print('protein aminoAcids =', protein.aminoAcids)

  print('protein function calls')
  print('protein name =', protein.getName())
  print('protein aminoAcids =', protein.getAminoAcids())
  print('protein sequence =', protein.getSequence())
  print('protein mass = ', protein.getMass())

