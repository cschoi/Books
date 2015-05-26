from math import log, exp

from MultipleAlign import profile, profileAlign
from Sequences import STANDARD_GENETIC_CODE as SGC
from Alignments import DNA_1, BLOSUM62, sequenceAlign, calcSeqSimilarity

AA_CATEGORIES = [('G','G'),('A','A'),('I','I'),('V','V'),
                 ('L','L'),('M','M'),('F','F'),('Y','Y'),
                 ('W','W'),('H','H'),('C','C'),('P','P'),
                 ('K','K'),('R','R'),('D','D'),('E','E'),
                 ('Q','Q'),('N','N'),('S','S'),('T','T'),
                 ('-','-'),
                 ('acidic',     'DE'),
                 ('hydroxyl',   'ST'),
                 ('aliphatic',  'VIL'),
                 ('basic',      'KHR'),
                 ('tiny',       'GAS'),
                 ('aromatic',   'FYWH'),
                 ('charged',    'KHRDE'),
                 ('small',      'AGSVTDNPC'),
                 ('polar',      'KHRDEQNSTC'),
                 ('hydrophobic','IVLFYWHAGMC'),
                 ('turnlike',   'GASHKRDEQNSTC'),
                 ('undef',      'AGVILMFYWHCPKRDEQNST-')]


def getConservation(align, simMatrix):

  conservation = []
  prof = profile(align)
  
  for compDict in prof:
    
    items = list(compDict.items())  # do not need list() in Python 2

    items.sort( key=lambda x: x[1] )
        
    score = 0.0
    
    for resA, compA in items:
      for resB, compB in items:
        score += compA * compB * simMatrix[resA][resB]
 
    bestLetter = items[-1][0]
    maxScore = simMatrix[bestLetter][bestLetter]
   
    score /= maxScore
    conservation.append(score)
  
  return conservation


def makeSimilarityString(align, simMatrix, thresholds):

  simString = ''
  conservation = getConservation(align, simMatrix)
  t1, t2, t3 = thresholds
  
  for score in conservation:
        
    if score >= t1:
      symbol = '*'
    elif score > t2:
      symbol = ':'
    elif score > t3: 
      symbol = '.'
    else:
      symbol = ' '
  
    simString += symbol
  
  return simString


def getAlignProperties(align, categories):

  properties = []
  
  prof = profile(align)

  for fracDict in prof:
  
    letters = fracDict.keys()

    for name, group in categories:

      for letter in letters:
        
        if letter not in group:
          break                     # quit inner loop

      else:                         # all letters are in group
        properties.append(name)     
        break 				# quit outer loop

  return properties   


def calcSubstitutionMatrix(alignments, alphabet, maxVal, smooth=5):

  matrix = {}
  counts = {}
  
  for letterA in alphabet:
    subDict = {}
    
    for letterB in alphabet:
      subDict[letterB] = 0
  
    matrix[letterA] = subDict
    counts[letterA] = 0
  
  totalRes = 0.0
  totalSub = 0.0

  for align in alignments:
 
    numPos = len(align[0])

    for i in range(numPos):
 
      letters = []
      
      for seq in align:

        letter = seq[i]
        if letter == '-':
          continue
    
        letters.append(letter)

      for letterA in letters:
        counts[letterA] += 1
      
        for letterB in letters:          
          matrix[letterA][letterB] += 1

      numLetters = len(letters)
      totalRes += numLetters    
      totalSub += numLetters * numLetters

  averageComp = {}    
  for letter in alphabet:
    averageComp[letter] = counts[letter]/totalRes      

  maxScore = None
  for resA in alphabet:
    for resB in alphabet:

      expected = averageComp[resA] * averageComp[resB]
      
      if not expected:
        continue

      observed = matrix[resA][resB]
      weight = 1.0 / (1.0+(observed/smooth))

      observed /= totalSub
      observed = weight*expected + (1-weight)*observed

      logOdds = log(observed/expected)
                  
      if (maxScore is None) or (logOdds>maxScore):
        maxScore = logOdds
      
      matrix[resA][resB] = logOdds

  maxScore = abs(maxScore)

  for resA in alphabet:
    for resB in alphabet:
      matrix[resA][resB] = int(maxVal*matrix[resA][resB]/maxScore)

  return matrix


def getDistanceMatrix(seqs, simMatrix):

  n = len(seqs)
  matrix = [[0.0] * n for x in range(n)]
  maxScores = [calcSeqSimilarity(x, x, simMatrix) for x in seqs]

  for i in range(n-1):
    seqA = seqs[i]
  
    for j in range(i+1,n):
      seqB = seqs[j]
      
      score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix)
      maxScore = max(maxScores[i],maxScores[j])
      dist = maxScore - score
      
      matrix[i][j] = dist
      matrix[j][i] = dist

  return matrix


def getJoinPair(distMatrix):

  n = len(distMatrix)

  minQ = None
  joinPair = None

  for i in range(n-1):
    sumRow = sum(distMatrix[i])
  
    for j in range(i+1, n):
      sumCol = sum(distMatrix[j])
    
      dist = distMatrix[i][j]
      q = (n-2)*dist - sumRow - sumCol
      
      if (minQ is None) or (q < minQ):
        minQ = q
        joinPair = [i,j]

  return joinPair


def getDistToJunction(distMatrix, i, j):
  
  n = len(distMatrix)
  row = distMatrix[i]
  column = distMatrix[j]
 
  dist = distMatrix[i][j] + (sum(row)-sum(column))/(n-2)
  dist *= 0.5

  return dist


def neighbourJoinTree(distMatrix):
 
  joinOrder = []
  n = len(distMatrix)
  tree = list(range(n))  # do not need list() in Python 2
  
  while n > 2:

    x, y = getJoinPair(distMatrix)

    node = (tree[x], tree[y])
    joinOrder.append(node)
    tree.append(node)

    del tree[y]
    del tree[x]

    distX = getDistToJunction(distMatrix, x, y)
    distY = getDistToJunction(distMatrix, y, x)
  
    distMatrix.append([0] * (n+1))

    for i in range(n):
      if i not in (x,y):

        dist = (distMatrix[x][i]-distX) + (distMatrix[y][i]-distY)
        dist *= 0.5
  
        distMatrix[i].append(dist)
        distMatrix[n][i] = dist

    del distMatrix[y]
    del distMatrix[x]
  
    for row in distMatrix:
      del row[y]
      del row[x]

    n -= 1

  tree = tuple(tree)
  joinOrder.append(tree)
  
  return tree, joinOrder


def treeProfileMultipleAlign(seqs, simMatrix):

  multipleAlign = {}
  for i, seq in enumerate(seqs):
    multipleAlign[i] = [seq,] 
  
  distMatrix = getDistanceMatrix(seqs, simMatrix)
  tree, treeOrder = neighbourJoinTree(distMatrix)
  joinStack = treeOrder[:]

  while joinStack:
    keyA, keyB  = joinStack.pop(0)
 
    subAlignA = multipleAlign[keyA]
    subAlignB = multipleAlign[keyB]

    profA = profile(subAlignA)
    profB = profile(subAlignB)

    s, profAlignA, profAlignB = profileAlign(profA, profB, simMatrix)

    gaps = []
    for i, fractions in enumerate(profAlignA):
      if fractions is None:
        gaps.append(i)
  
    for j, seq in enumerate(subAlignA):
      for gap in gaps:
        seq = seq[:gap] + '-' + seq[gap:]
        
      subAlignA[j] = seq

    gaps = []
    for i, fractions in enumerate(profAlignB):
      if fractions is None:
        gaps.append(i)
    
    for j, seq in enumerate(subAlignB):
      for gap in gaps:
        seq = seq[:gap] + '-' + seq[gap:]
        
      subAlignB[j] = seq
    
    newKey = (keyA, keyB)  
    multipleAlign[newKey] = subAlignA + subAlignB

    del multipleAlign[keyA]
    del multipleAlign[keyB]

  return multipleAlign[newKey], tree, treeOrder


def ancestorResidue(codesA, codesB):

  common = codesA.intersection(codesB)

  if common:
    return common, False

  else:
   union = codesA.union(codesB)
 
   if codesA and codesB:
     return union, True
   else:
     return union, False


def calcSubstitutionRates(align, treeOrder):

  n = len(align[0])
  numNodes = float(len(treeOrder))
  absRates = [0.0] * n
  relRates = [0.0] * n
  
  treeSeqs = {}
  for i, seq in enumerate(align):
    sets = []

    for letter in seq:
      if letter == '-':
        sets.append(set())
      else:
        sets.append(set([letter]))

    treeSeqs[i] = sets
  
  while treeOrder:

    a, b = treeOrder.pop(0)
    seqA = treeSeqs[a]
    seqB = treeSeqs[b]
    seqC = []
   
    for i in range(n):

      residueSet, swapped = ancestorResidue(seqA[i], seqB[i])
      seqC.append(residueSet)
      
      if swapped:
        absRates[i] += 1.0
      
    treeSeqs[(a, b)] = seqC
    del treeSeqs[a]
    del treeSeqs[b]
  
  meanRate = sum(absRates)/(numNodes*n)
     
  for i in range(n):
    rate = absRates[i] / numNodes
    absRates[i] = rate
    relRates[i] = (rate-meanRate)/meanRate
      
  return absRates, relRates


def calcActivePassive(align, treeOrder, geneticCode):

  n = len(align[0])
  numNodes = float(len(treeOrder))
  active = 0
  passive = 0
  
  treeSeqs = {}
  for i, seq in enumerate(align):
    sets = []

    for letter in seq:
      if letter == '-':
        sets.append(set())
      else:
        sets.append(set([letter]))

    treeSeqs[i] = sets
  
  while treeOrder:

    a, b = treeOrder.pop(0)
    seqA = treeSeqs[a]
    seqB = treeSeqs[b]
    seqC = []
    
    for i in range(n):

      residues, swapped = ancestorResidue(seqA[i], seqB[i])
      seqC.append(residues)
      
      if swapped:
        codonStart = (i//3)*3 # // is integer division
        
        subSeqA = seqA[codonStart:codonStart+3]
        subSeqB = seqB[codonStart:codonStart+3]

        aminoAcidsA = set()
        for x in subSeqA[0]:
          for y in subSeqA[1]:
            for z in subSeqA[2]:
              codon = x+y+z
              codon.replace('T','U')
              aminoAcidsA.add( geneticCode.get(codon) )

        aminoAcidsB = set()
        for x in subSeqB[0]:
          for y in subSeqB[1]:
            for z in subSeqB[2]:
              codon = x+y+z
              codon.replace('T','U')
              aminoAcidsB.add( geneticCode.get(codon) )

        if aminoAcidsA.intersection(aminoAcidsB):
          passive += 1
        else:
          active += 1  
          
    treeSeqs[(a, b)] = seqC
    del treeSeqs[a]
    del treeSeqs[b]
      
  return active, passive


if __name__ == '__main__':

  print("\nConservation scores along alignment")

  align1 = ['AAGCCGCACACAGACCCTGAG',
            'AAGCTGCACGCAGACCCTGAG',
            'AGGCTGCACGCAGACCCTGAG',
            'AAGCTGCACGTGGACCCTGAG',
            'AGGCTGCACGTGGACCCTGAG',
            'AGGCTGCACGTGGACCCTGAG',
            'AAGCTGCATGTGGACCCTGAG']
 
  conserv1 = getConservation(align1, DNA_1)

  print([round(v, 2) for v in conserv1])

  align2 = ['QPVHPFSRPAPVVIILIILCVMAGVIGTILLISYGIRLLIK-------------',
            'QLVHRFTVPAPVVIILIILCVMAGIIGTILLISYTIRRLIK-------------',
            'QLAHHFSEPE---ITLIIFGVMAGVIGTILLISYGIRRLIKKSPSDVKPLPSPD',
            'QLVHEFSELV---IALIIFGVMAGVIGTILFISYGSRRLIKKSESDVQPLPPPD',
            'MLEHEFSAPV---AILIILGVMAGIIGIILLISYSIGQIIKKRSVDIQPPEDED',
            'PIQHDFPALV---MILIILGVMAGIIGTILLISYCISRMTKKSSVDIQSPEGGD',
            'QLVHIFSEPV---IIGIIYAVMLGIIITILSIAFCIGQLTKKSSLPAQVASPED',
            '-LAHDFSQPV---ITVIILGVMAGIIGIILLLAYVSRRLRKRP-----PADVP-',
            'SYHQDFSHAE---ITGIIFAVMAGLLLIIFLIAYLIRRMIKKPLPVPKPQDSPD']

  conserv2 = getConservation(align2, BLOSUM62)

  print([round(v, 2) for v in conserv2])


  print("\nMultiple alignment with similarity symbols")

  symbols = makeSimilarityString(align2, BLOSUM62, (1.0, 0.5, 0.3))

  for seq in align2:
    print(seq)

  print(symbols)

  print("\nConsensus amino acid categories")
  
  catNames = getAlignProperties(align2, AA_CATEGORIES)

  for i, category in enumerate(catNames):
    print(i, category)

  
  print("\nCalculate substitution matrix from alignment")

  aminoAcids = BLOSUM62.keys()
  print(calcSubstitutionMatrix([align2,], aminoAcids, 10))


  seqs = ['QPVHPFSRPAPVVIILIILCVMAGVIGTILLISYGIRLLIK',
          'QLVHRFTVPAPVVIILIILCVMAGIIGTILLISYTIRRLIK',
          'QLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKKSPSDVKPLPSPD',
          'QLVHEFSELVIALIIFGVMAGVIGTILFISYGSRRLIKKSESDVQPLPPPD',
          'MLEHEFSAPVAILIILGVMAGIIGIILLISYSIGQIIKKRSVDIQPPEDED',
          'PIQHDFPALVMILIILGVMAGIIGTILLISYCISRMTKKSSVDIQSPEGGD',
          'QLVHIFSEPVIIGIIYAVMLGIIITILSIAFCIGQLTKKSSLPAQVASPED',
          'LAHDFSQPVITVIILGVMAGIIGIILLLAYVSRRLRKRPPADVP',
          'SYHQDFSHAEITGIIFAVMAGLLLIIFLIAYLIRRMIKKPLPVPKPQDSPD']


  distMatrix = getDistanceMatrix(seqs, BLOSUM62)
  distMatrix = getDistanceMatrix(align1, DNA_1)

  for row in distMatrix:
    print(['%.2f' % x for x in row])

  
  print("\nGet closest sequence pair from a distance matrix")

  print(getJoinPair(distMatrix))

  
  print("\nEstimate pylogenetic tree by neighbour joining")

  distMatrix = getDistanceMatrix(seqs, BLOSUM62)
  tree, treeJoinOrder = neighbourJoinTree(distMatrix)
 

  print(tree) # Result : (((7, (0, 1)), (4, 5)), ((2, 3), (6, 8)))


  print("\nMultiple alignment using pylogenetic tree")
 
  align, tree, order = treeProfileMultipleAlign(seqs, BLOSUM62)

  for seq in align:
    print(seq)

  print("\nCalculate residue substitution rates from pylogenetic tree")

  align, tree, joinOrder = treeProfileMultipleAlign(seqs, BLOSUM62)
  absRates, relRates = calcSubstitutionRates(align, joinOrder)

  for i, absRate in enumerate(absRates):
    print('%2d %6.2f %6.2f' % (i, absRate, relRates[i]))


  print("\nCalculate passive and active DNA substitutions with tree")

  seqs = ['AAAGTGGATGAAGTTGGTGCTGAGGCCCTGGGCAGGCTG',
          'AAAGTGGATGATGTTGGTGCTGAGGCCCTGGGCAGGCTG',
          'AAAGTGGATGAAGTTGGTGCTGAGGCCCTGGGCAGGCTG',
          'AAAGTGGATGAAGTTGGTGCTGAAGCCCTGGGCAGGCTG',
          'AAAGTGGATGAAGTTGGTGCTGAGGCCCTGGGCAGGCTG',
          'CATGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTG',
          'AAAGTGGACGAAGTTGGTGCTGAGGCCCTGGGCAGGCTG',
          'CATGTGGATGAAATTAGTGGTGAGGTCCTGGGCAGGCTG',
          'AACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTG']


  align, tree, joinOrder = treeProfileMultipleAlign(seqs, DNA_1)
  active, passive = calcActivePassive (align, joinOrder, SGC)

  print('Active: %d Passive: %d' % (active, passive))


