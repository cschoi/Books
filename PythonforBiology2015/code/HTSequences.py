FTP_ROOT = 'ftp://ftp.ncbi.nlm.nih.gov/genomes'
ALIGNER_PATH = '/home/user/programs/bowtie-0.12.7/'

import os, gzip
from fnmatch import fnmatch
from subprocess import call
from matplotlib import pyplot
from numpy import array

try:
  # Python 3
  from urllib.request import urlopen
except ImportError:
  # Python 2
  from urllib2 import urlopen

try:
  from HTSeq import FastqReader, SAM_Reader, GFF_Reader
  from HTSeq import GenomicArray, GenomicInterval
  HAVE_HTSEQ = True
except ImportError:
  print('HTSeq library inaccessible or not installed. Some examples will not work')
  HAVE_HTSEQ = False

def downloadFile(remoteFile, localFile):

  print('Downloading to %s ...' % localFile)
  
  response = urlopen(remoteFile)
  data = response.read()
  fileObj = open(localFile, 'wb')
  fileObj.write(data)
  fileObj.close()
  
  print(' ...done')


def downloadGenomeFiles(remoteDir, localDir, fileType='*.fna', url=FTP_ROOT):
  
  if remoteDir[0] != '/':
    remoteDir = '/' + remoteDir

  if remoteDir[-1] != '/':
    remoteDir = remoteDir + '/'
  
  remotePath = url + remoteDir
  print("Reading %s" % remotePath)
  
  response = urlopen(remotePath)
  data = response.read()

  fileNames = []
  filePaths = []
  chromosomeDirs = []
  lines = data.split('\n')

  for line in lines:
    line = line.strip()
    
    if not line:
      continue

    fileName = line.split()[-1]
    
    if fileName.startswith('CHR'):
      chromosomeDirs.append(fileName + '/')

    elif fnmatch(fileName, fileType):
      fileNames.append(fileName)
      continue

  for fileName in fileNames:
    filePath = os.path.join(localDir, fileName)
    filePath = os.path.abspath(filePath)
    filePaths.append(filePath)
    downloadFile(url + remoteDir + fileName, filePath)

  for chromosomeDir in chromosomeDirs:
    subDir = remoteDir + chromosomeDir
    filePaths += downloadGenomeFiles(subDir, localDir,  fileType, url)

  return filePaths


def uncompressGzFile(fileName):

  if fileName.endswith('.gz'):
    inFileObj = gzip.open(fileName, 'rb')
 
    fileName = fileName[:-3]
    outFileObj = open(fileName, 'w')

    for line in inFileObj:
      outFileObj.write(line)
    
    # Faster alternative, given sufficient memory:
    # outFileObj.write(infileObj.read())

    inFileObj.close()
    outFileObj.close()

  return fileName


def indexGenome(genomeName, fileNames, outputDir,
                tableSize=10, quiet=True, pack=True):

  fastaFiles = []
  for fileName in fileNames:
    fileName = uncompressGzFile(fileName)
    fastaFiles.append(fileName)

  fastaFileStr= ','.join(fastaFiles)

  cmdArgs = [ALIGNER_PATH+'bowtie-build', '-f']
  
  if quiet:
    cmdArgs.append('-q')
  
  if pack:
    cmdArgs.append('-p')

  cmdArgs += ['-t', str(tableSize), fastaFileStr, genomeName]

  print(' '.join(cmdArgs))
  call(cmdArgs, cwd=outputDir)

def genomeAlign(genomeName, genomeDir, readFiles, pairedReadFiles=None, 
                outFile=None, leftTrim=0, rightTrim=0, useSOLiD=False, 
                qualType='sanger', maxMismatches=2, highQualLen=28, 
                pairSepRange=(0,250), showHits=1, maxHits=2,
                enforceBest=True, cpuCores=None):

  os.environ['BOWTIE_INDEXES'] = genomeDir
  
  if not outFile:
    outFile = genomeName + '.sam'

  cmdArgs = [ALIGNER_PATH+'bowtie', '-S']
  readFilesStr = ','.join(readFiles)

  foreName, fileType = os.path.splitext(readFiles[0])

  if fileType in ('fa','fna','mfa','fasta'):
    cmdArgs.append('-f')
    
  elif fileType in ('fq','fastq'):
    cmdArgs.append('-q')
  
  else:
    cmdArgs.append('-r')

  if useSOLiD:
    cmdArgs.append('-C')
  
  if qualType == 'illumina1.3': # Phred+64
    cmdArgs.append('--phred64-quals')
  
  elif qualType == 'solexa': # Phred+59
    cmdArgs.append('--solexa-quals')
  
  else: # Sanger, current illumina:  Phred + 33
    cmdArgs.append('--phred33-quals')  

  if enforceBest:
    cmdArgs.append('--best')  
  
  if not cpuCores:
    import multiprocessing
    cpuCores = multiprocessing.cpu_count()

  cmdArgs += ['-5', str(leftTrim),
              '-3', str(rightTrim),
              '-n', str(maxMismatches),
              '-l', str(highQualLen),
              '-k', str(showHits),
              '-m', str(maxHits),
              '-p', str(cpuCores),
              '--chunkmbs', '256']

  if pairedReadFiles:
    pairedReadFiles = ','.join(pairedReadFiles)
    minSep, maxSep = pairSepRange
    cmdArgs += [genomeName,
                '--minins', str(minSep),
                '--maxins', str(maxSep),
                '-1', readFilesStr,
                '-2', pairedReadFiles,
                outFile]
    
  else:
    cmdArgs += [genomeName, readFilesStr, outFile]

  call(cmdArgs)
  
  return outFile

if __name__ == '__main__':

  print('\n Download E.coli genome\n')
  
  remoteGenome = FTP_ROOT+'/Bacteria/Escherichia_coli_536_uid58531/NC_008253.fna'
  downloadFile(remoteGenome, 'EcoliGenome.fasta')
  
  # Human genome (big)
  #filePaths = downloadGenomeFiles('H_sapiens','examples','hs_ref*.fa.gz')


  print('\n Indexing E.coli genome, prior to read mapping\n')

  filePath = os.path.abspath('examples/EcoliGenome.fasta')

  genomeName = 'E_coli'

  indexGenome(genomeName, [filePath,], 'examples')


  print('\n Map sequence reads to E.coli genome\n')

  fastqFile = 'examples/EcoliReads.fastq'

  genomeAlign(genomeName, 'examples', [fastqFile], qualType='sanger')


  if not HAVE_HTSEQ:
    import sys
    print('\n Exiting early. HTSeq library not available.')
    sys.exit(0)


  print('\n Using HTSeq library with FASTQ reads\n')

  fileObj = FastqReader(fastqFile)

  for seqRead in fileObj:
    print(seqRead.name)
    print(seqRead.seq)
    print(seqRead.get_reverse_complement()[::-1])

  numReads = 0.0
  meanQual = None

  for seqRead in fileObj:
    print(seqRead.qual)

    if meanQual is None:
      meanQual = array(seqRead.qual)
    else:
      meanQual += seqRead.qual
    numReads += 1.0

  if numReads:
    pyplot.plot(meanQual/numReads)
    pyplot.show()


  print('\n Using HTSeq library with SAM alignments\n')

  alignFile = 'examples/EcoliGenomeAlign.sam'
  chromosomes = set()

  for alignment in SAM_Reader(alignFile):
 
    if alignment.aligned:
      seqRead = alignment.read
      print(seqRead.name)
      print(seqRead.seq)

      genomeRegion = alignment.iv
      chromo = genomeRegion.chrom
      strand = genomeRegion.strand
      start = genomeRegion.start
      end = genomeRegion.end
    
      chromosomes.add(chromo)
      print(chromo, start, end, strand)


  print('\n HTSeq genomic arrays and intervals\n')

  chromosomes = list(chromosomes)
  hitMap = GenomicArray(chromosomes, stranded=True, typecode='i')

  for alignment in SAM_Reader(alignFile):
 
    if alignment.aligned:
      genomeRegion = alignment.iv
 
      if genomeRegion.strand == '+':
        hitMap[genomeRegion] = 1
      else:
        hitMap[genomeRegion] = -1

  chromo = chromosomes[0]
  endPoint = 2000000
  plusStrand  = GenomicInterval(chromo, 0, endPoint, '+')
  minusStrand = GenomicInterval(chromo, 0, endPoint, '-')
  bothStrands = GenomicInterval(chromo, 0, endPoint, '.')

  pyplot.plot(list(hitMap[plusStrand]))
  pyplot.plot(list(hitMap[minusStrand]))
  pyplot.show()


  print('\n Using HTSeq to access GFF genome features\n')

  remoteFileName = '/Bacteria/Escherichia_coli_536_uid58531/NC_008253.gff'
  gffFile = 'examples/EcoliGenomeFeatures.gff'
  downloadFile(FTP_ROOT+remoteFileName, gffFile)

  fileObj = GFF_Reader(gffFile)

  for genomeFeature in fileObj:

    genomeRegion = genomeFeature.iv

    data = (genomeRegion.chrom,
            genomeRegion.start,
            genomeRegion.end,
            genomeRegion.strand)

    print('%s %s - %s (%s)' % data)

    data = (genomeFeature.name,
          genomeFeature.type,
          genomeFeature.source)
 
    print('%s %s (%s)' % data)

    print(genomeFeature.attr)

    
  print('\n Using HTSeq to plot features in genomic arrays\n')

  geneMap = GenomicArray('auto', stranded=False, typecode='O')
  genePlot = GenomicArray('auto', stranded=False, typecode='i')

  for genomeFeature in fileObj:
 
    if genomeFeature.type == 'gene':
 
      genomeRegion = genomeFeature.iv
      geneMap[genomeRegion] = genomeFeature
      genePlot[genomeRegion] = 1

  for region, feature in geneMap.steps():
 
    if feature:
 
      data = (feature.name,
              region.start,
              region.end,
              feature.iv.strand)
 
      print('%s: %s - %s (%s)' % data)

  chromosome = genomeRegion.chrom
  region = GenomicInterval(chromosome, 0, 40000, '.')
  pyplot.plot(list(genePlot[region]))
  pyplot.show()

