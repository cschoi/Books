
def consensus(alignment, threshold=0.25):

  n = len(alignment[0])
  nSeq = float(len(alignment))
  consensus = ''

  for i in range(n):
    counts = {}

    for seq in alignment:
      letter = seq[i]
      if letter == '-':
        continue
    
      counts[letter] = counts.get(letter, 0) + 1
    
    fractions = []
    for letter in counts:
      frac = counts[letter]/nSeq
      fractions.append([frac, letter])
      
    fractions.sort()
    bestFraction, bestLetter = fractions[-1]
    
    if bestFraction <  threshold:
      consensus += 'X'
    
    else:
      consensus += bestLetter

  return consensus


def profile(alignment):

  n = len(alignment[0])
  nSeq = float(len(alignment))
  prof = []
  
  for i in range(n):
    counts = {}
    
    for seq in alignment:
      letter = seq[i]
      if letter == '-':
        continue
    
      counts[letter] = counts.get(letter, 0) + 1
    
    for letter in counts:
      counts[letter] /= nSeq
    
    prof.append(counts)

  return prof

  
def profileAlign(profileA, profileB, simMatrix, insert=8, extend=4):
 
  numI = len(profileA) + 1
  numJ = len(profileB) + 1
  
  scoreMatrix = [[0] * numJ for x in range(numI)]
  routeMatrix = [[0] * numJ for x in range(numI)]
  
  for i in range(1, numI):
    routeMatrix[i][0] = 1
  
  for j in range(1, numJ):
    routeMatrix[0][j] = 2
  
  for i in range(1, numI):
    for j in range(1, numJ):
    
      penalty1 = insert
      penalty2 = insert
      
      if routeMatrix[i-1][j] == 1:
        penalty1 = extend
        
      elif routeMatrix[i][j-1] == 2:
        penalty2 = extend
        
      fractionsA = profileA[i-1]
      fractionsB = profileB[j-1]

      similarity = 0.0
      totalWeight = 0.0
      for residueA in fractionsA:
        for residueB in fractionsB:
          weight = fractionsA[residueA] * fractionsB[residueB]
          totalWeight += weight
          similarity += weight * simMatrix[residueA][residueB]
      
      penalty1 *= totalWeight
      penalty2 *= totalWeight
      
      paths = [scoreMatrix[i-1][j-1] + similarity, # Route 0 
               scoreMatrix[i-1][j] - penalty1, #Route 1
               scoreMatrix[i][j-1] - penalty2]  # Route 2
               
      best = max(paths)
      route = paths.index(best)           

      scoreMatrix[i][j] = best
      routeMatrix[i][j] = route
      
  profileOutA = []
  profileOutB = []
  
  i = numI-1
  j = numJ-1
  score = scoreMatrix[i][j]
    
  while i > 0 or j > 0:
    route = routeMatrix[i][j]    
      
    if route == 0: # Diagonal
      profileOutA.append(profileA[i-1])
      profileOutB.append(profileB[j-1])
      i -= 1
      j -= 1
      
    elif route == 1: # Gap in profile B
      profileOutA.append(profileA[i-1])
      profileOutB.append(None)
      i -= 1      
  
    elif route == 2: # Gap in profile A
      profileOutA.append(None)
      profileOutB.append(profileB[j-1])
      j -= 1      
  
  profileOutA.reverse()
  profileOutB.reverse()

  return score, profileOutA, profileOutB 

def simpleProfileMultipleAlign(seqs, simMatrix):

  n = len(seqs)

  score, alignA, alignB = sequenceAlign(seqs[0], seqs[1], simMatrix)
  multipleAlign = [alignA, alignB]

  for i in range(2,n):
    profA = profile(multipleAlign)
    toAdd = [seqs[i],]
    profB = profile(toAdd)

    score, alignA, alignB = profileAlign(profA, profB, simMatrix)

    gaps = []
    for j, fractions in enumerate(alignA):
      if fractions is None:
        gaps.append(j)
      
    for j, seq in enumerate(multipleAlign):
      for gap in gaps:
        seq = seq[:gap] + '-' + seq[gap:]
        
      multipleAlign[j] = seq

    gaps = []
    for j, fractions in enumerate(alignB):
      if fractions is None:
        gaps.append(j)
      
    for j, seq in enumerate(toAdd):
      for gap in gaps:
        seq = seq[:gap] + '-' + seq[gap:]
        
      toAdd[j] = seq
      
    multipleAlign.extend(toAdd)
       
  return multipleAlign
  

def consensusMultipleAlign(seqs, threshold, simMatrix):

  n = len(seqs)
  multipleAlign = []
  
  i = 0
  for j in range(i+1,n):
    seqB = seqs[j]

    if not multipleAlign:
      seqA = seqs[i]
      score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix)
      multipleAlign.append(alignA)
      multipleAlign.append(alignB)
    
    else:
      seqA = consensus(multipleAlign, threshold)
      score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix)
      
      gaps = []
      for k, letter in enumerate(alignA):
        if letter == '-':
          gaps.append(k)
      
      for k, seq in enumerate(multipleAlign):
        for gap in gaps:
          seq = seq[:gap] + '-' + seq[gap:]
        
        multipleAlign[k] = seq
      
      multipleAlign.append(alignB)
      
  for k, seq in enumerate(multipleAlign):
    print(k, seq)


if __name__ == '__main__':

  print('\nConsensus sequence')
  
  alignment = ['SRPAPVVIILIILCVMAGVIGTILLISYGIRLLIK',
               'TVPAPVVIILIILCVMAGIIGTILLISYTIRRLIK',
               'HHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIK',
               'HEFSELVIALIIFGVMAGVIGTILFISYGSRRLIK']

  print(consensus(alignment))
  
  print('\nSequence profile dictonary')

  alignment = ['SRPAPVVIILIILCVMAGVIGTILLISYGIRLLIK',
               'TVPAPVVIILIILCVMAGIIGTILLISYTIRRLIK',
               'HHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIK',
               'HEFSELVIALIIFGVMAGVIGTILFISYGSRRLIK']

  print(profile(alignment))

  # First sub-dict: {'H': 0.5, 'S': 0.25, 'T': 0.25},

  alignA = ['SRPAPVV--LII',
            'TVPAPVVIILII']
  alignB = ['HHFSEPEITLIIF',
            'H-FSELVIALIIF']

  from Alignments import BLOSUM62, sequenceAlign

  score, profA, profB = profileAlign( profile(alignA), profile(alignB), BLOSUM62)

  print('\nProfile alignment with None as gaps')
  
  print(score)
  print(profA)
  print(profB)

  seqs = ['SRPAPVVLIILCVMAGVIGTILLISYGIRLLIK',
          'TVPAPVVIILIILCVMAGIIGTILLLIISYTIRRLIK',
          'HHFSEPEITLIIFGVMAGVIGTILLLIISYGIRLIK',
          'HFSELVIALIIFGVMAGVIGTILFISYGSRLIK']

  align = simpleProfileMultipleAlign(seqs, BLOSUM62)
  for k, seq in enumerate(align):
    print(k, seq)

  print('\nConsensus paired alignment')

  consensusMultipleAlign(seqs, 0.25, BLOSUM62)


  try:
    from Bio import SeqIO
 
  except ImportError:
    print('\n* * * BioPython not installed * * * * * * * * \n')
    print('\n* * * Remaining examples will not work  * * * \n')
    
    import sys
    sys.exit(0)
    
  import os
  
  from Bio.Seq import Seq
  from Bio.SeqRecord import SeqRecord
  from Bio.Alphabet import IUPAC
  from Bio import AlignIO
  from subprocess import call

  fastaFileName = "test2.fasta"
  alignFileName = "test2.aln"


  records = []
  for i, seq in enumerate(seqs):
    seqObj = Seq(seq, IUPAC.protein)
    name = 'test%d' % i
    recordObj = SeqRecord(seqObj, id=name, description='demo only')
    records.append(recordObj)

  outFileObj = open(fastaFileName, "w")
  SeqIO.write(records, outFileObj, "fasta")
  outFileObj.close()

  cmdArgs = ['clustalw',
             '-INFILE=' + fastaFileName,
             '-OUTFILE=' + alignFileName]
  call(cmdArgs)

  fileObj = open(alignFileName)
  alignment = AlignIO.read(fileObj, "clustal")
  
  print('\nClustalW alignment\n')
  
  print("Alignment length %i" % alignment.get_alignment_length())
  for record in alignment:
    print(record.seq, record.id)

 
  alignments = [alignment,]
  outputHandle = open("test2.phylip", "w")
  AlignIO.write(alignments, outputHandle, "phylip")




