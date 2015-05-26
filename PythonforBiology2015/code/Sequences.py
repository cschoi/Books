
STANDARD_GENETIC_CODE = { 
          'UUU':'Phe', 'UUC':'Phe', 'UCU':'Ser', 'UCC':'Ser',
          'UAU':'Tyr', 'UAC':'Tyr', 'UGU':'Cys', 'UGC':'Cys',
          'UUA':'Leu', 'UCA':'Ser', 'UAA':None,  'UGA':None,
          'UUG':'Leu', 'UCG':'Ser', 'UAG':None,  'UGG':'Trp',
          'CUU':'Leu', 'CUC':'Leu', 'CCU':'Pro', 'CCC':'Pro',
          'CAU':'His', 'CAC':'His', 'CGU':'Arg', 'CGC':'Arg',
          'CUA':'Leu', 'CUG':'Leu', 'CCA':'Pro', 'CCG':'Pro',
          'CAA':'Gln', 'CAG':'Gln', 'CGA':'Arg', 'CGG':'Arg',
          'AUU':'Ile', 'AUC':'Ile', 'ACU':'Thr', 'ACC':'Thr',
          'AAU':'Asn', 'AAC':'Asn', 'AGU':'Ser', 'AGC':'Ser',
          'AUA':'Ile', 'ACA':'Thr', 'AAA':'Lys', 'AGA':'Arg',
          'AUG':'Met', 'ACG':'Thr', 'AAG':'Lys', 'AGG':'Arg',
          'GUU':'Val', 'GUC':'Val', 'GCU':'Ala', 'GCC':'Ala',
          'GAU':'Asp', 'GAC':'Asp', 'GGU':'Gly', 'GGC':'Gly',
          'GUA':'Val', 'GUG':'Val', 'GCA':'Ala', 'GCG':'Ala', 
          'GAA':'Glu', 'GAG':'Glu', 'GGA':'Gly', 'GGG':'Gly'}

GES_SCALE = {'F':-3.7,'M':-3.4,'I':-3.1,'L':-2.8,'V':-2.6,
             'C':-2.0,'W':-1.9,'A':-1.6,'T':-1.2,'G':-1.0,
             'S':-0.6,'P': 0.2,'Y': 0.7,'H': 3.0,'Q': 4.1,
             'N': 4.8,'E': 8.2,'K': 8.8,'D': 9.2,'R':12.3}

def proteinTranslation(seq, geneticCode):
  """ This function translates a nucleic acid sequence into a
      protein sequence, until the end or until it comes across
      a stop codon """
   
  seq = seq.replace('T','U') # Make sure we have RNA sequence
  proteinSeq = []

  i = 0
  while i+2 < len(seq):
    codon = seq[i:i+3]
    aminoAcid = geneticCode[codon]

    if aminoAcid is None: # Found stop codon
      break

    proteinSeq.append(aminoAcid)
    i += 3

  return proteinSeq
 
 
def estimateMolMass(seq, molType='protein'):
  """Calculate the molecular weight of a biological sequence assuming
     normal isotopic ratios and protonation/modification states
  """

  residueMasses = {"DNA":
                   {"G":329.21,"C":289.18,"A":323.21,"T":304.19,},
                   "RNA":
                   {"G":345.21,"C":305.18,"A":329.21,"U":302.16,},
                   "protein":
                   {"A": 71.07,"R":156.18,"N":114.08, "D":115.08,
                    "C":103.10,"Q":128.13, "E":129.11,"G": 57.05,
                    "H":137.14, "I":113.15,"L":113.15,"K":128.17,
                    "M":131.19,"F":147.17,"P": 97.11, "S": 87.07,
                    "T":101.10,"W":186.20, "Y":163.17,"V": 99.13}}
  
  massDict = residueMasses[molType]

  # Begin with mass of extra end atoms H + OH
  molMass = 18.02

  for letter in seq:
    molMass += massDict.get(letter, 0.0)

  return molMass


def matchDnaProfile(seq, profile):
  """ Find the best matching position and score when comparing a DNA      
      sequence with a DNA sequence profile""" 

  bestScore = 0
  bestPosition = None # Just to start with

  width = len(profile['A'])

  for i in range(len(seq)-width):
    score = 0

    for j in range(width):
      letter = seq[i+j]
      score += profile[letter][j]
    
    if score > bestScore:
      bestScore = score
      bestPosition = i

  return bestScore, bestPosition


def calcGcContent(seq, winSize=10):
  
  gcValues = []
  
  for i in range(len(seq)-winSize):

    subSeq = seq[i:i+winSize]

    numGc = subSeq.count('G') + subSeq.count('C')

    value = numGc/float(winSize)

    gcValues.append(value)

  return gcValues


def hydrophobicitySearch(seq, scale, winSize=15):
  """Scan a protein sequence for hydrophobic regions using the GES
     hydrophobicity scale.
  """

  score = None
  scoreList  = []

  for i in range(len(seq)- winSize):
    j = i + winSize
    
    if score is None:
      score = 0
      for k in range(i,j):
        score += scale[seq[k]]

    else:
      score += scale[seq[j-1]]
      score -= scale[seq[i-1]]

    scoreList.append(score)

  return scoreList


def relativeEntropySearch(seq, winSize, isProtein=False):
  """Scan a sequence for repetitiveness by calculating relative
     information entropy.
  """

  lenSeq = len(seq)
  scores = [0.0] * lenSeq
  
  extraSeq = seq[:winSize]
  seq += extraSeq

  if isProtein:
    resCodes = 'ACDEFGHIKLMNPQRSTVWY'
  else:
    resCodes = 'GCAT'

  for i in range(lenSeq):
    subSeq = seq[i:i+winSize]
    scores[i] = calcRelativeEntropy(subSeq, resCodes)

  return scores


def estimateCharge(sequence, pH):
  """Using pKa values estimate the charge of a sequence of
     amino acids at a given pH"""

  pKaDict = {'+': 8.0,'-': 3.1,'K':10.0,'R':12.0,
             'H': 6.5,'E': 4.4,'D': 4.4,'Y':10.0,'C': 8.5}

  isAcid = {'+':False,'-':True,'K':False,'R':False,
            'H':False,'E':True,'D':True,'Y':True,'C':True}

  total = 0.0
  
  for aminoAcid in sequence:
    pKa = pKaDict.get(aminoAcid)

    if pKa is not None:
      r = 10.0 ** (pH-pKa)
      dissociated = r/(r+1.0)

      if isAcid[aminoAcid]:
        charge = -1.0 * dissociated 
      else:
        charge = 1.0 - dissociated 

      total += charge

  return total


def estimateIsoelectric(sequence):
  """Estimate the charge neutral pH of a protein sequence.
     This is just a guess as pKa values will vary according to 
     protein sequence, conformation and conditions.
  """

  sequence = '+' + sequence + '-' # assumes seq is a string
  bestValue = 0.0
  minCharge = estimateCharge(sequence, bestValue)
  increment = 7.0

  while abs(minCharge) > 0.001:
    pHtest = bestValue + increment
    charge = estimateCharge(sequence, pHtest)

    if abs(charge) < abs(minCharge):
      minCharge = charge
      bestValue = pHtest


    else:
      increment = abs(increment)/2.0
      if minCharge < 0.0:
        increment *= -1

  return bestValue
  
  
def calcRelativeEntropy(seq, resCodes):
  """Calculate a relative entropy value for the residues in a     
     sequence compared to a null hypothesis where each residue
     appears in a randomly and unbiased.
  """
  
  from math import log

  N = float(len(seq))

  base = 1.0/len(resCodes)
  
  prop = {}
  for r in resCodes:
    prop[r] = 0

  for r in seq:
    prop[r] += 1

  for r in resCodes:
    prop[r] /= N

  DKL = 0
  for r in resCodes:
    if prop[r] != 0.0:
      d = prop[r]* log(prop[r]/base, 2.0)
      DKL += d

  return DKL


def relativeEntropySearch(seq, winSize, isProtein=False):
  """Scan a sequence for repetitiveness by calculating relative
     information entropy.
  """

  lenSeq = len(seq)
  scores = [0.0] * lenSeq
  
  extraSeq = seq[:winSize]
  seq += extraSeq

  if isProtein:
    resCodes = 'ACDEFGHIKLMNPQRSTVWY'
  else:
    resCodes = 'GCAT'

  for i in range(lenSeq):
    subSeq = seq[i:i+winSize]
    scores[i] = calcRelativeEntropy(subSeq, resCodes)

  return scores


def estimateCharge(sequence, pH):
  """Using pKa values estimate the charge of a sequence of
     amino acids at a given pH"""

  pKaDict = {'+': 8.0,'-': 3.1,'K':10.0,'R':12.0,
             'H': 6.5,'E': 4.4,'D': 4.4,'Y':10.0,'C': 8.5}

  isAcid = {'+':False,'-':True,'K':False,'R':False,
            'H':False,'E':True,'D':True,'Y':True,'C':True}

  total = 0.0
  
  for aminoAcid in sequence:
    pKa = pKaDict.get(aminoAcid)

    if pKa is not None:
      r = 10.0 ** (pH-pKa)
      dissociated = r/(r+1.0)

      if isAcid[aminoAcid]:
        charge = -1.0 * dissociated 
      else:
        charge = 1.0 - dissociated 

      total += charge

  return total


def estimateIsoelectric(sequence):
  """Estimate the charge neutral pH of a protein sequence.
     This is just a guess as pKa values will vary according to 
     protein sequence, conformation and conditions.
  """

  sequence = '+' + sequence + '-' # assumes seq is a string
  bestValue = 0.0
  minCharge = estimateCharge(sequence, bestValue)
  increment = 7.0

  while abs(minCharge) > 0.001:
    pHtest = bestValue + increment
    charge = estimateCharge(sequence, pHtest)

    if abs(charge) < abs(minCharge):
      minCharge = charge
      bestValue = pHtest


    else:
      increment = abs(increment)/2.0
      if minCharge < 0.0:
        increment *= -1

  return bestValue

if __name__ == '__main__':

  dnaSeq = 'ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTG'

  protein3LetterSeq = proteinTranslation(dnaSeq, STANDARD_GENETIC_CODE)

  rnaSeq = dnaSeq.replace('T','U')

  proteinSeq = 'IRTNGTHMQPLLKLMKFQKFLLELFTLQKRKPEKGYNLPIISLNQ'
  proteinMass = estimateMolMass(proteinSeq)

  dnaMass = estimateMolMass(dnaSeq, molType='DNA')

  seq = 'AGCTCGCTCGCTGCGTATAAAATCGCATCGCGCGCAGC'
  position1 = seq.find('TATAAA')
  position2 = seq.find('GAGGAG')

  profile = {
    'A':[ 61, 16,352,  3,354,268,360,222,155, 56, 83, 82, 82, 68, 77],
    'C':[145, 46,  0, 10,  0,  0,  3,  2, 44,135,147,127,118,107,101],
    'G':[152, 18,  2,  2,  5,  0, 10, 44,157,150,128,128,128,139,140],
    'T':[ 31,309, 35,374, 30,121,  6,121, 33, 48, 31, 52, 61, 75, 71]}
 
  score, position = matchDnaProfile(dnaSeq, profile)
  print(score, position, dnaSeq[position:position+15])

  from matplotlib import pyplot

  gcResults = calcGcContent(dnaSeq)

  pyplot.plot(gcResults)
  pyplot.show()


  from matplotlib import pyplot

  scores = hydrophobicitySearch(proteinSeq, GES_SCALE)

  pyplot.plot(scores)
  pyplot.show()


  from matplotlib import pyplot

  dnaScores = relativeEntropySearch(dnaSeq, 6)

  proteinScores = relativeEntropySearch(proteinSeq, 10, isProtein=True)

  pyplot.plot(dnaScores)
  pyplot.plot(proteinScores)
  pyplot.show()

  try:
    from Bio import SeqIO
 
  except ImportError:
    print('\n* * * BioPython not installed * * * * * * * * \n')
    print('\n* * * Remaining examples will not work  * * * \n')
    
    import sys
    sys.exit(0)
    

  fileObj = open("examples/demoSequences.fasta", "rU")

  for protein in SeqIO.parse(fileObj, 'fasta'):
    print(protein.id)
    print(protein.seq)
    print(estimateIsoelectric(protein.seq))

  fileObj.close()

  from Bio.SeqRecord import SeqRecord
  from Bio.Seq import Seq
  from Bio.Alphabet import IUPAC

  fileObj = open("output.fasta", "w")


  seqObj = Seq(proteinSeq, IUPAC.protein)
  proteinObj = SeqRecord(seqObj, id="TEST")

  SeqIO.write([proteinObj,], fileObj,  'fasta')

  fileObj.close()


  from Bio import Entrez

  Entrez.email = 'mickey@disney.com'
  socketObj = Entrez.efetch(db="protein", rettype="fasta",
                            id="71066805")

  dnaObj = SeqIO.read(socketObj, "fasta")
  socketObj.close()

  print(dnaObj.description)
  print(dnaObj.seq)

  from Bio import ExPASy

  socketObj = ExPASy.get_sprot_raw('HBB_HUMAN')
  proteinObj = SeqIO.read(socketObj, "swiss")
  socketObj.close()

  print(proteinObj.description)
  print(proteinObj.seq)

