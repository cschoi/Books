from re import sub
from subprocess import call, Popen, PIPE
from xml.etree import ElementTree

DNA_1 = {'G': { 'G':1, 'C':0, 'A':0, 'T':0 },
         'C': { 'G':0, 'C':1, 'A':0, 'T':0 },
         'A': { 'G':0, 'C':0, 'A':1, 'T':0 },
         'T': { 'G':0, 'C':0, 'A':0, 'T':1 }}

          
REV_COMP = {'G': { 'G':-1, 'C': 1, 'A':-1, 'T':-1 },
            'C': { 'G': 1, 'C':-1, 'A':-1, 'T':-1 },
            'A': { 'G':-1, 'C':-1, 'A':-1, 'T': 1 },
            'T': { 'G':-1, 'C':-1, 'A': 1, 'T':-1 }}

          
DNA_2 = {'G': { 'G': 1, 'C':-3, 'A':-3, 'T':-3, 'N':0 },
         'C': { 'G':-3, 'C': 1, 'A':-3, 'T':-3, 'N':0 },
         'A': { 'G':-3, 'C':-3, 'A': 1, 'T':-3, 'N':0 },
         'T': { 'G':-3, 'C':-3, 'A':-3, 'T': 1, 'N':0 },
         'N': { 'G': 0, 'C': 0, 'A': 0, 'T': 0, 'N':0 }}  

BLOSUM62 = {'A':{'A': 4,'R':-1,'N':-2,'D':-2,'C': 0,'Q':-1,'E':-1,'G': 0,'H':-2,'I':-1,
                 'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S': 1,'T': 0,'W':-3,'Y':-2,'V': 0,'X':0},
            'R':{'A':-1,'R': 5,'N': 0,'D':-2,'C':-3,'Q': 1,'E': 0,'G':-2,'H': 0,'I':-3,
                 'L':-2,'K': 2,'M':-1,'F':-3,'P':-2,'S':-1,'T':-1,'W':-3,'Y':-2,'V':-3,'X':0},
            'N':{'A':-2,'R': 0,'N': 6,'D': 1,'C':-3,'Q': 0,'E': 0,'G': 0,'H': 1,'I':-3,
                 'L':-3,'K': 0,'M':-2,'F':-3,'P':-2,'S': 1,'T': 0,'W':-4,'Y':-2,'V':-3,'X':0},
            'D':{'A':-2,'R':-2,'N': 1,'D': 6,'C':-3,'Q': 0,'E': 2,'G':-1,'H':-1,'I':-3,
                 'L':-4,'K':-1,'M':-3,'F':-3,'P':-1,'S': 0,'T':-1,'W':-4,'Y':-3,'V':-3,'X':0},
            'C':{'A': 0,'R':-3,'N':-3,'D':-3,'C': 9,'Q':-3,'E':-4,'G':-3,'H':-3,'I':-1,
                 'L':-1,'K':-3,'M':-1,'F':-2,'P':-3,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1,'X':0},
            'Q':{'A':-1,'R': 1,'N': 0,'D': 0,'C':-3,'Q': 5,'E': 2,'G':-2,'H': 0,'I':-3,
                 'L':-2,'K': 1,'M': 0,'F':-3,'P':-1,'S': 0,'T':-1,'W':-2,'Y':-1,'V':-2,'X':0},
            'E':{'A':-1,'R': 0,'N': 0,'D': 2,'C':-4,'Q': 2,'E': 5,'G':-2,'H': 0,'I':-3,
                 'L':-3,'K': 1,'M':-2,'F':-3,'P':-1,'S': 0,'T':-1,'W':-3,'Y':-2,'V':-2,'X':0},
            'G':{'A': 0,'R':-2,'N': 0,'D':-1,'C':-3,'Q':-2,'E':-2,'G': 6,'H':-2,'I':-4,
                 'L':-4,'K':-2,'M':-3,'F':-3,'P':-2,'S': 0,'T':-2,'W':-2,'Y':-3,'V':-3,'X':0},
            'H':{'A':-2,'R': 0,'N': 1,'D':-1,'C':-3,'Q': 0,'E': 0,'G':-2,'H': 8,'I':-3,
                 'L':-3,'K':-1,'M':-2,'F':-1,'P':-2,'S':-1,'T':-2,'W':-2,'Y': 2,'V':-3,'X':0},
            'I':{'A':-1,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-3,'E':-3,'G':-4,'H':-3,'I': 4,
                 'L': 2,'K':-3,'M': 1,'F': 0,'P':-3,'S':-2,'T':-1,'W':-3,'Y':-1,'V': 3,'X':0},
            'L':{'A':-1,'R':-2,'N':-3,'D':-4,'C':-1,'Q':-2,'E':-3,'G':-4,'H':-3,'I': 2,
                 'L': 4,'K':-2,'M': 2,'F': 0,'P':-3,'S':-2,'T':-1,'W':-2,'Y':-1,'V': 1,'X':0},
            'K':{'A':-1,'R': 2,'N': 0,'D':-1,'C':-3,'Q': 1,'E': 1,'G':-2,'H':-1,'I':-3,
                 'L':-2,'K': 5,'M':-1,'F':-3,'P':-1,'S': 0,'T':-1,'W':-3,'Y':-2,'V':-2,'X':0},
            'M':{'A':-1,'R':-1,'N':-2,'D':-3,'C':-1,'Q': 0,'E':-2,'G':-3,'H':-2,'I': 1,
                 'L': 2,'K':-1,'M': 5,'F': 0,'P':-2,'S':-1,'T':-1,'W':-1,'Y':-1,'V': 1,'X':0},
            'F':{'A':-2,'R':-3,'N':-3,'D':-3,'C':-2,'Q':-3,'E':-3,'G':-3,'H':-1,'I': 0,
                 'L': 0,'K':-3,'M': 0,'F': 6,'P':-4,'S':-2,'T':-2,'W': 1,'Y': 3,'V':-1,'X':0},
            'P':{'A':-1,'R':-2,'N':-2,'D':-1,'C':-3,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-3,
                 'L':-3,'K':-1,'M':-2,'F':-4,'P': 7,'S':-1,'T':-1,'W':-4,'Y':-3,'V':-2,'X':0},
            'S':{'A': 1,'R':-1,'N': 1,'D': 0,'C':-1,'Q': 0,'E': 0,'G': 0,'H':-1,'I':-2,
                 'L':-2,'K': 0,'M':-1,'F':-2,'P':-1,'S': 4,'T': 1,'W':-3,'Y':-2,'V':-2,'X':0},
            'T':{'A': 0,'R':-1,'N': 0,'D':-1,'C':-1,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-1,
                 'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S': 1,'T': 5,'W':-2,'Y':-2,'V': 0,'X':0},
            'W':{'A':-3,'R':-3,'N':-4,'D':-4,'C':-2,'Q':-2,'E':-3,'G':-2,'H':-2,'I':-3,
                 'L':-2,'K':-3,'M':-1,'F': 1,'P':-4,'S':-3,'T':-2,'W':11,'Y': 2,'V':-3,'X':0},
            'Y':{'A':-2,'R':-2,'N':-2,'D':-3,'C':-2,'Q':-1,'E':-2,'G':-3,'H': 2,'I':-1,
                 'L':-1,'K':-2,'M':-1,'F': 3,'P':-3,'S':-2,'T':-2,'W': 2,'Y': 7,'V':-1,'X':0},
            'V':{'A': 0,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-2,'G':-3,'H':-3,'I': 3,
                 'L': 1,'K':-2,'M': 1,'F':-1,'P':-2,'S':-2,'T': 0,'W':-3,'Y':-1,'V': 4,'X':0},
            'X':{'A': 0,'R': 0,'N': 0,'D': 0,'C': 0,'Q': 0,'E': 0,'G': 0,'H': 0,'I': 0,
                 'L': 0,'K': 0,'M': 0,'F': 0,'P': 0,'S': 0,'T': 0,'W': 0,'Y': 0,'V': 0,'X':0}}

def calcSeqIdentity(seqA, seqB):

  numPlaces = min(len(seqA), len(seqB))
  
  score = 0.0
  
  for i in range(numPlaces):
  
    if seqA[i] == seqB[i]:
       score += 1.0

  return 100.0 * score/numPlaces


def calcSeqSimilarity(seqA, seqB, simMatrix):

  numPlaces = min(len(seqA), len(seqB))
  
  totalScore = 0.0
  
  for i in range(numPlaces):
    
    residueA = seqA[i]
    residueB = seqB[i]
  
    totalScore += simMatrix[residueA][residueB]

  return totalScore


def pairAlignScore(alignA, alignB, simMatrix, insert=8, extend=4):

  totalScore = 0.0
  
  n = min( len(alignA), len(alignB) )
  
  for i in range(n):
    
    residueA = alignA[i]
    residueB = alignB[i]
    
    if '-' not in (residueA, residueB):
      simScore = simMatrix[residueA][residueB]
      
    elif (i > 0) and ('-' in (alignA[i-1], alignB[i-1])):
      simScore = -extend
    else:
      simScore = -insert
    
    totalScore += simScore
  
  return totalScore


def sequenceAlign(seqA, seqB, simMatrix=DNA_2, insert=8, extend=4):

  numI = len(seqA) + 1
  numJ = len(seqB) + 1

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
        
      similarity = simMatrix[ seqA[i-1] ][ seqB[j-1] ]
      
      paths = [scoreMatrix[i-1][j-1] + similarity, # Route 0
               scoreMatrix[i-1][j] - penalty1, # Route 1
               scoreMatrix[i][j-1] - penalty2] # Route 2                     
      
      best = max(paths)
      route = paths.index(best)           

      scoreMatrix[i][j] = best
      routeMatrix[i][j] = route
      
  alignA = []
  alignB = []
  
  i = numI-1
  j = numJ-1
  score = scoreMatrix[i][j]

    
  while i > 0 or j > 0:
    route = routeMatrix[i][j]    

    if route == 0: # Diagonal
      alignA.append( seqA[i-1] )
      alignB.append( seqB[j-1] )
      i -= 1
      j -= 1

    elif route == 1: # Gap in seqB
      alignA.append( seqA[i-1] )
      alignB.append( '-' )
      i -= 1      

    elif route == 2: # Gap in seqA
      alignA.append( '-' )
      alignB.append( seqB[j-1] ) 
      j -= 1
  
  alignA.reverse()
  alignB.reverse()
  alignA = ''.join(alignA)
  alignB = ''.join(alignB)

  return score, alignA, alignB 
  
"""  
def makeBlastDatabase(fastaFile, databaseName,
                      formatdbExe=None, isProtein=True):

  if not formatdbExe:
    formatdbExe = 'formatdb'

  dbType = 'T' if isProtein else 'F'

  cmdArgs = [formatdbExe,
             '-p', dbType,
             '-i', fastaFile,
             '-n', databaseName]
    
  print('Making BLAST database %s...' % databaseName)
  
  try:
    call(cmdArgs)
  
  except Exception as err:
    print('BLAST database creation failed')
    print('Command used: "%s"' % ' '.join(cmdArgs))
    print(err)
    return

  print(' ...done')
"""

def makeBlastDatabase(fastaFile, databaseName,
                     formatdbExe=None, isProtein=True):
  """Make a BLAST sequence database from a FASTA file"""
 
  if not formatdbExe:
    formatdbExe = 'makeblastdb'
  
  if call(['which', formatdbExe]) != 0:
    print 'Executable file %s not found' % formatdbExe
  
  if isProtein:
    molType = 'prot'
  else:
    molType = 'nucl'

  cmdArgs = [formatdbExe,
             '-dbtype', molType,
             '-in', fastaFile,
             '-out', databaseName]
    
  print 'Making BLAST database %s...' % databaseName
  
  try:
    call(cmdArgs)
  
  except Exception as err:
    print('BLAST database creation failed')
    print('Command used: "%s"' % ' '.join(cmdArgs))
    print(err)
    return
  
  print ' ...done'

"""
def blastSearch(seq, database, blastExe=None, eCutoff=10,
                matrix='BLOSUM62', cpuCores=None,
                proteinQuery=True, proteinDatabase=True):
  
  if not blastExe:
    blastExe = 'blast2'
   
  if not cpuCores:
    import multiprocessing
    cpuCores = multiprocessing.cpu_count()

  if proteinQuery:
    if proteinDatabase:
      program = 'blastp'
    else:
      program = 'tblastn'
  
  else:
    if proteinDatabase:
      program = 'blastx'
    else:
      program = 'blastn'  

  querySeq = seq.upper()
  querySeq = sub('(.{60})(.)',r'\1\n\2', querySeq)
  querySeq = '>UserQuery\n%s\n' % querySeq 
  
  cmdArgs = [blastExe,
             '-p', program,
             '-m', '7',             # => XML output
             '-a', str(cpuCores),
             '-d', database,
             '-e', str(eCutoff),
             '-M', matrix]
  
  print(' '.join(cmdArgs))
  
  try:
    proc = Popen(cmdArgs, stdin=PIPE, stdout=PIPE)
    
  except Exception as err:
    print('BLAST command failed')
    print('Command used: "%s"' % ' '.join(cmdArgs))
    print(err)
    return []
  
  stdOutData, stdErrData = proc.communicate(querySeq)
 
  results = []
  root = ElementTree.fromstring(stdOutData)

  iteration = root.find('BlastOutput_iterations').find('Iteration')

  ihit = iteration.find('Iteration_hits')
  
  if ihit:
    hits = ihit.findall('Hit')
  else:
    return []

  for hit in hits:
    hsp = hit.find('Hit_hsps').find('Hsp')

    name = hit.findtext('Hit_def')

    hitDict = {}

    for tag in ('Hsp_bit-score', 'Hsp_evalue'):
      key = tag[4:]
      hitDict[key] = float(hsp.findtext(tag))

    for tag in ('Hsp_score', 'Hsp_query-from',
                'Hsp_query-to', 'Hsp_hit-from',
                'Hsp_hit-to', 'Hsp_query-frame',
                'Hsp_hit-frame', 'Hsp_identity',
                'Hsp_positive', 'Hsp_gaps',
                'Hsp_align-len'):
      key = tag[4:]
      hitDict[key] = int(hsp.findtext(tag, '0'))

    for tag in ('Hsp_qseq', 'Hsp_hseq', 'Hsp_midline'):
      key = tag[4:]
      hitDict[key] = hsp.findtext(tag)

    results.append( (name, hitDict) )
 
  return results
"""  
  
def blastSearch(seq, database, bastExe=None,
                eCutoff=10, matrix='BLOSUM62',
                cpuCores=None): 
  """Run BLAST from Python"""
  
  if not bastExe:
    bastExe = 'blastp'
    
  if not cpuCores:
    import multiprocessing
    cpuCores = multiprocessing.cpu_count()

  # make input
  
  querySeq = seq.upper()
  querySeq = sub('(.{60})(.)',r'\1\n\2', querySeq)
  querySeq = '>UserQuery\n%s\n' % querySeq
    
  # run BLAST
  
  cmdArgs = [bastExe,
             '-outfmt', '5', # => XML outupt
             '-num_threads', str(cpuCores),
             '-db', database, #'-o', outFile,
             '-evalue', str(eCutoff),
             '-matrix', matrix]
  
  #print ' '.join(cmdArgs)
  
  try:
    proc = Popen(cmdArgs, stdin=PIPE, stdout=PIPE)
    
  except Exception as err:
    print('BLAST command failed')
    print('Command used: "%s"' % ' '.join(cmdArgs))
    print(err)
    return []
  
  stdOutData, stdErrData = proc.communicate(querySeq)

  results = []
  root = ElementTree.fromstring(stdOutData)

  iteration1 = root.find('BlastOutput_iterations').find('Iteration')
  if iteration1 is None:
    return []
  
  ihits = iteration1.find('Iteration_hits')
  if ihits is None:
    return []

  hits = ihits.findall('Hit')

  for hit in hits:
    hsp = hit.find('Hit_hsps').find('Hsp')
 
    hitDict = {}
    hitDict['def'] = hit.findtext('Hit_def')
    hitDict['len'] = int(hit.findtext('Hit_len'))
 
    for tag in ('Hsp_bit-score', 'Hsp_evalue'):
      key = tag[4:]
      hitDict[key] = float(hsp.findtext(tag))
 
    for tag in ('Hsp_score', 'Hsp_query-from',
                'Hsp_query-to', 'Hsp_hit-from',
                'Hsp_hit-to', 'Hsp_query-frame',
                'Hsp_hit-frame', 'Hsp_identity',
                'Hsp_positive', 'Hsp_gaps',
                'Hsp_align-len'):
      key = tag[4:]
      hitDict[key] = int(hsp.findtext(tag, '0'))

    for tag in ('Hsp_qseq', 'Hsp_hseq', 'Hsp_midline'):
      key = tag[4:]
      hitDict[key] = hsp.findtext(tag)
    
    results.append( hitDict )
 
  return results
  
if __name__ == '__main__':

  seq1 = 'ALIGNMENTS'
  seq2 = 'ALIGDVENTS'
  seq3 = 'ALIGDPVENTS'
  seq4 = 'ALIGN-MENTS'

  print(calcSeqIdentity(seq1, seq2)) # 80.0%
  print(calcSeqIdentity(seq1, seq3)) # 40.0%
  print(calcSeqIdentity(seq4, seq3)) # 72.7%

  # Test with pre-defined substitution matrices
  # DNA example
  print(calcSeqSimilarity('AGCATCGCTCT', 'AGCATCGTTTT', DNA_2))

  # Protein example
  print(calcSeqSimilarity('ALIGNMENT', 'AYIPNVENT', BLOSUM62))

  # Test
  print(pairAlignScore('ALIGDPPVENTS', 'ALIGN--MENTS', BLOSUM62)) # 28
  print(pairAlignScore('ALIGDPPVENTS', '--ALIGNMENTS', BLOSUM62)) # -3


  seqA = 'WFSEPEIST'
  seqB = 'FSRPAVVIST'

  score, alignA, alignB = sequenceAlign(seqA, seqB, BLOSUM62)

  print(score)  # 17
  print(alignA) # WFSEPE--IST
  print(alignB) # -FSRPAVVIST

  fileName = 'examples/EcoliProteome.fasta'

  makeBlastDatabase(fileName, 'ECOLI_PROT')
 
  seq = 'NGTISYTNEAGKIYQLKPNPAVLICRVRGLHLPEKHVTWRGEAIPGSLFDFA' \
        'LYFFHNYQALLAKGSGPYFYLPKTQSWQEAAWWSEVFSYAEDRFNLPRGTIK' \
        'ATLLIETLPAVFQMDEILHALRDHIVGLNCGRWDYIFSYIKTLKNYPDRVLP'
        
  database = 'ECOLI_PROT'
  
  # "blastp" is location of the program executable
  results = blastSearch(seq, database, 'blastp', eCutoff=1.0, cpuCores=4)
  
  for hitDict in results:
    print(hitDict['def'], hitDict['score'], hitDict['evalue'])
    print(hitDict['qseq'])
    print(hitDict['midline'])
    print(hitDict['hseq'])
