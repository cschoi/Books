import sys, os, pickle, csv
from xml.etree import cElementTree as ElementTree


def readFastaFile(fileName):

  fileObj = open(fileName, 'rU')
  sequences = []
  seqFragments = []

  for line in fileObj:
    if line.startswith('>'):
      # found start of next sequence
      if seqFragments:
        sequence = ''.join(seqFragments)
        sequences.append(sequence)
      seqFragments = []

    else:
      # found more of existing sequence
      seq = line.rstrip() # remove newline character
      seqFragments.append(seq)

  if seqFragments:
    # should be the case if file is not empty
    sequence = ''.join(seqFragments)
    sequences.append(sequence)

  fileObj.close()

  return sequences


def readFastaFile(fileName):

  fileObj = open(fileName, "rU")
  sequences = []
  seq = ''

  for line in fileObj:
    if line.startswith('>'):
      if seq:
        sequences.append(seq)
      seq = ''

    else:
      seq += line.rstrip() 

  if seq:
    sequences.append(seq)

  fileObj.close()

  return sequences


def calcCentroid(pdbFile):

  fileObj = open(pdbFile, "rU")

  natoms = 0
  xsum = ysum = zsum = 0

  for line in fileObj:
    if line[:6] == 'ATOM  ':
      natoms = natoms + 1
      x = float(line[30:38])
      y = float(line[38:46])
      z = float(line[46:54])
      xsum += x
      ysum += y
      zsum += z
    
  fileObj.close()

  if natoms == 0:
    xavg = yavg = zavg = 0
  else:
    xavg = xsum / natoms
    yavg = ysum / natoms
    zavg = zsum / natoms

  return (natoms, xavg, yavg, zavg)



def printPubmedAbstracts(xmlFiles):
  
  for xmlFile in xmlFiles:
    tree = ElementTree.parse(xmlFile)
    root = tree.getroot()

    citationElem = root.find('PubmedArticle/MedlineCitation')
    pmid = citationElem.findtext('PMID')
    articleElem = citationElem.find('Article')
    journalElem = articleElem.find('Journal')
    journalTitle = journalElem.findtext('Title')
    journalYear = journalElem.findtext('JournalIssue/PubDate/Year')
    articleTitle = articleElem.findtext('ArticleTitle')
    articleAbstract = articleElem.findtext('Abstract/AbstractText')

    print('PMID = %s' % pmid)
    print('journalYear = %s' % journalYear)
    print('journalTitle = "%s"' % journalTitle)
    print('articleTitle = "%s"' % articleTitle)
    print('articleAbstract = "%s"' % articleAbstract)
    print('')


def writeSingleFastaSequence(comment, sequence, fastaFile):

  fileObj = open(fastaFile, "w")

  fileObj.write('> %s\n' % comment)

  for (n, code) in enumerate(sequence):
    if n > 0 and n % 60 == 0:
      fileObj.write('\n')

    fileObj.write(code)

  if n % 60 != 0:
    fileObj.write('\n')

  fileObj.close()


def writeFastaSeqs(comments, sequences, fastaFile, width=60):

  fileObj = open(fastaFile, "w")

  for i, seq in enumerate(sequences):
      
    numLines = 1 + (len(seq)-1)//width
    seqLines = [seq[width*x:width*(x+1)] for x in range(numLines)]
      
    seq = '\n'.join(seqLines)
    fileObj.write('> %s\n%s\n' % (comments[i], seq))

  fileObj.close()


def writeListFile(fileName, data, headings, formats, separator='\t'):

  if len(data[0]) != len(headings):
    print("Headings length does not match input list")
    return

  fileObj = open(fileName, 'w')
  
  line = separator.join(headings)
  fileObj.write('%s\n' % line)

  format = separator.join(formats)
  for row in data:
   line = format % tuple(row)
   fileObj.write('%s\n' % line)

  fileObj.close()


def writeChromosomeRegions(fileName, data):
   
  headings = ['chromo', 'start', 'end', 'value']
  formats = ['%s', '%d', '%d', '%.3f']
  writeListFile(fileName, data, headings, formats, ' ')


def readListFile(fileName, converters, separator='\t'):
  
  dataList = []
  fileObj = open(fileName, 'rU')
  header = fileObj.readline()     # Extract first line

  for line in fileObj:           # Loop through remaining lines
    line = line.rstrip()

    data = line.split(separator)

    for index, datum in enumerate(data):
      convertFunc = converters[index]
     
      if convertFunc:
        data[index] = convertFunc(datum)

    dataList.append(data)
   
  return dataList


def readChromosomeRegions(fileName):

  converters = [None, int, int, float]
  dataList = readListFile(fileName, converters, ' ')

  return dataList


def writeCsvFile(fileName, data, headings, separator='\t'):

  if sys.version_info.major > 2:
    fileObj = open(fileName, 'w', newline='')
  else:
    fileObj = open(fileName, 'wb')
    
  writer = csv.writer(fileObj, delimiter=separator)
  writer.writerow(headings)

  for row in data:
    writer.writerow(row)
  fileObj.close()
  
  
def readCsvFile(fileName, converters, separator='\t'):

  dataList = []
  if sys.version_info.major > 2:
    fileObj = open(fileName, 'r', newline='')
  else:
    fileObj = open(fileName, 'rb')
  reader = csv.reader(fileObj, delimiter=separator)

  for n, row in enumerate(reader):
    if n > 0:  # n = 0 is the header, which we ignore
      for index, datum in enumerate(row):
        convertFunc = converters[index]

        if convertFunc:
          row[index] = convertFunc(datum)

      dataList.append(row)

  fileObj.close()

  return dataList
  

def findFiles(directory, suffix):

  files = []
  dirfiles = os.listdir(directory)
    
  for dirfile in dirfiles:
    fullfile = os.path.join(directory, dirfile)

    if os.path.isdir(fullfile):
      # fullfile is a directory, so recurse into that
      files.extend(findFiles(fullfile, suffix))

    elif dirfile.endswith(suffix):
      # fullfile is a normal file, and with correct suffix
      files.append(fullfile)

  return files


def removeFiles(directory, suffix):

  dirfiles = os.listdir(directory)

  for dirfile in dirfiles:
    fullfile = os.path.join(directory, dirfile)

    if os.path.isdir(fullfile):
      # fullfile is a directory, so recurse into that
      removeFiles(fullfile, suffix)

    elif dirfile.endswith(suffix):
      # fullfile is a normal file, and with correct suffix
      os.remove(fullfile)



if __name__ == '__main__':

  # Open and close a file object
  
  path = 'examples/dataFile.txt'
  fileObj = open(path)
  fileObj.close()
  
  
  
  # Reading as a data block or as a line
  
  fileObj = open(path, 'r')

  data = fileObj.read()
  line = fileObj.readline()
  
  
  
  # Parsing all lines - explicit line reads
  
  fileObj = open(path, "rU")
  line = fileObj.readline()
  while line:
    # process line
    line = fileObj.readline()
  fileObj.close()
  
  
  
  # Parsing all lines - iterator looping through file object
  
  fileObj = open(path, "rU")
  for line in fileObj:
    pass # process line

  fileObj.close()


  
  # Reading all lines at once
  
  fileObj = open(path, "rU")
  lines = fileObj.readlines()
  fileObj.close()
  for line in lines:
    pass # process line
  
  # or shorhand
  
  for line in open(path, "rU").readlines():
    pass # process line
  
  
  
  # Iterator for old Python versions
  
  import fileinput

  iterator = fileinput.input(path)

  for line in iterator:
    print(line)

  fileinput.close()

  
  
  # Getting command line options/arguments
  
  pyScriptName = sys.argv[0] # 'programFile.py'
  
  """
  dataFileName = sys.argv[1] # 'data/inputFile.txt'


  def workFunction(*args):
    pass

  workFunction(dataFileName) # Use the file name for something?
  """


  # Reading whitespace separated file
  
  fileObj = open('examples/chromoData.tsv')
  values = []
  header = fileObj.readline() # Don't need this first line

  for line in fileObj:
    data = line.split()

    chromosome, position, value = data
    position = int(position)
    value = float(value)

    values.append(value)

  mean = sum(values)/len(values)
  print('Mean value', mean)
  
  
  
  # Reading FASTA sequence file
  
  fileName = 'examples/demoSequences.fasta'
  sequences = readFastaFile(fileName)
  print(sequences)
  
  
  
  # Reading and working with data from column file (PDB)
  
  print(calcCentroid('examples/protein.pdb'))

  
  
  # Basic XML reading operation
  
  xmlFile = 'examples/protein.xml'

  tree = ElementTree.parse(xmlFile)
  root = tree.getroot()


  """
  # Working with XML nodes
  text = node.findtext(pattern)

  element = node.find(pattern)
  if element:
    text = element.text
  else:
    text = None
  """
  

  # Get PUbMed abstracts from XML file
  
  #xmlFiles = ['examples/pubmed_24009655.xml', 'examples/pubmed_24009636.xml']
  #printPubmedAbstracts(xmlFiles)
  
  
  
  # Creating an open file object for writing 'w' or appending 'a'
  
  path = 'output.txt'
  fileObj = open(path, "w")
  fileObj = open(path, "a")
  fileObj.close()
  
  
  
  # Writing FASTA sequence fils

  comment = 'test sequence'
  sequence = sequences[0]
  fastaFile = 'testOut1.fasta'
  writeSingleFastaSequence(comment, sequence, fastaFile)

  comments = len(sequences)*[comment]
  fastaFile = 'testOut2.fasta'
  writeFastaSeqs(comments, sequences, fastaFile)
  
  
  
  # Read and write simple delimited text files

  fileName = 'testFile.txt'
  headings = ['x', 'y']
  data = [[1, 2], [3, 4], [5, 6]]
  formats = ['%d', '%d']
  writeListFile(fileName, data, headings, formats)

  converters = 2*[int]
  dataList = readListFile(fileName, converters)
  print(dataList)

  
  
  # Read and write character separated files 

  fileName = 'testChromoRegions.txt'

  data = [['chr1',195612601,196518584,0.379],
          ['chr1',52408393,196590488,0.361],
          ['chr1',193237929,196783789,0.473],
          ['chr1',181373059,6104731,0.104],
          ['chr2',7015693,7539562,0.508],
          ['chr2',9097449,9108209,0.302]]
        
  writeChromosomeRegions(fileName, data)

  dataList = readChromosomeRegions(fileName)
  print(dataList)

  
  # CSV files
  
  headings = ['chromo', 'start', 'end', 'value']
  converters = [None, int, int, float]
  
  fileName = 'testChromoRegions.csv'
  writeCsvFile(fileName, data, headings)
  dataList = readCsvFile(fileName, converters)
  print(dataList)
  
  
  # Saving serialised Python objects as pickle files

  data = [1,2,3,4,'a','b','c']
 
  fileObj = open('saveFile.pickle', 'wb')

  pickle.dump(data, fileObj)

  fileObj.close()



  # Reading Python objects from serialised pickle files

  fileObj = open('saveFile.pickle', 'rb')

  data = pickle.load(fileObj)

  fileObj.close()
  


  # Recursion into directories to find files
  
  directory = '.'
  suffix = '.py'
  for fileName in findFiles(directory, suffix):
    print fileName



  # Recursion into directories to remove files
  
  removeFiles(directory, '.pickle')
