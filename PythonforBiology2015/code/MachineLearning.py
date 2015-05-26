from time import time 

from matplotlib import pyplot
from PIL import Image

import numpy
from numpy import exp, power, array, zeros, sqrt 
from numpy import nonzero, random, abs, sum, dot
from numpy import tanh, ones, append


def getFeatureDistance(vector1, vector2):
  
  distance = 0.0  
  for a, b in zip(vector1, vector2):
    delta = a-b
    distance += delta * delta 
  
  return distance


def kNearestNeighbour(knowns, query, k=7):

  if k >= len(knowns):
    raise Exception('Length of training data must be larger than k')

  dists = []
  for vector, cat in knowns[:k]:
    dist = getFeatureDistance(vector, query)
    dists.append( (dist, cat) ) # Vector could be included
  
  dists.sort()
  closest = dists[:k]
  
  counts = {}
  for dist, cat in closest:
    counts[cat] = counts.get(cat, 0) + 1

  bestCount = max(counts.values())
  bestCats = [cat for cat in counts if counts[cat] == bestCount]

  for dist, cat in closest:
    if cat in bestCats:
      return cat

def kNearestNeighbour(knowns, query, k=7):
  # Numpy version
  
  if k >= len(knowns):
    raise Exception('Length of training data must be larger than k')
  
  vectors, cats = zip(*knowns)
  deltas = array(vectors) - array(query)
  dists = (deltas * deltas).sum(axis=1)
  sortedIndices = dists.argsort()
  closeCats = [cats[i] for i in sortedIndices[:k]]
  
  counts = {}
  for cat in closeCats:
    counts[cat] = counts.get(cat, 0) + 1

  bestCount = max(counts.values())
  bestCats = [cat for cat in counts if counts[cat] == bestCount]

  for cat in closeCats:
    if cat in bestCats:
      return cat


def selfOrganisingMap(inputs, spread, size, steps=1000):

  nRows, nCols = size
  vecLen = len(inputs[0])
  somap = numpy.random.rand(nRows, nCols, vecLen)
  	
  influence = numpy.dstack([spread]*vecLen) # One for each feature
  infWidth = (len(spread)-1) // 2
  makeMesh = numpy.ix_   # Ugly

  for s in range(steps):  # xrange in Python 2
    
    decay = exp(-s/float(steps))
    
    for vector in inputs:
      
      diff =  somap-vector
      diff2 = diff*diff
      dist2 = diff2.sum(axis=2)

      index = dist2.argmin()
      row = index // nRows
      col = index % nRows

      rows = [x % nRows for x in range(row-infWidth, row+1+infWidth)]
      cols = [y % nCols for y in range(col-infWidth, col+1+infWidth)]

      mesh = makeMesh(rows,cols)
      somap[mesh] -= diff[mesh] * influence * decay
          

  return somap


def neuralNetPredict(inputVec, weightsIn, weightsOut):

  signalIn = append(inputVec, 1.0) # input layer

  prod = signalIn * weightsIn.T
  sums = sum(prod, axis=1)
  signalHid = tanh(sums)  # hidden  layer

  prod = signalHid * weightsOut.T
  sums = sum(prod, axis=1)
  signalOut = tanh(sums)  # output  layer

  return signalIn, signalHid, signalOut
 
 
def neuralNetTrain(trainData, numHid, steps=100, rate=0.5, momentum=0.2):

  numInp = len(trainData[0][0])
  numOut = len(trainData[0][1])
  numInp += 1
  minError = None

  sigInp = ones(numInp)
  sigHid = ones(numHid)
  sigOut = ones(numOut)

  wInp = random.random((numInp, numHid))-0.5
  wOut = random.random((numHid, numOut))-0.5
  bestWeightMatrices = (wInp, wOut)

  cInp = zeros((numInp, numHid))
  cOut = zeros((numHid, numOut))

  for x, (inputs, knownOut) in enumerate(trainData):
    trainData[x] = (array(inputs), array(knownOut))
 
  for step in range(steps):  # xrange in Python 2
    random.shuffle(trainData) # Important
    error = 0.0
 
    for inputs, knownOut in trainData:
      sigIn, sigHid, sigOut = neuralNetPredict(inputs, wInp, wOut)

      diff = knownOut - sigOut
      error += sum(diff * diff)

      gradient = ones(numOut) - (sigOut*sigOut)
      outAdjust = gradient * diff 

      diff = sum(outAdjust * wOut, axis=1)
      gradient = ones(numHid) - (sigHid*sigHid)
      hidAdjust = gradient * diff 

      # update output 
      change = outAdjust * sigHid.reshape(numHid, 1)
      wOut += (rate * change) + (momentum * cOut)
      cOut = change
 
      # update input 
      change = hidAdjust * sigIn.reshape(numInp, 1)
      wInp += (rate * change) + (momentum * cInp)
      cInp = change
 
    if (minError is None) or (error < minError):
      minError = error
      bestWeightMatrices = (wInp.copy(), wOut.copy())
      print("Step: %d Error: %f" % (step, error))
  
  return bestWeightMatrices


def convertSeqToVector(seq, indexDict):
  
  numLetters = len(indexDict)
  vector = [0.0] * len(seq) * numLetters

  for pos, letter in enumerate(seq):
    index = pos * numLetters + indexDict[letter]  
    vector[index] = 1.0

  return vector


def kernelGauss(vectorI, vectorJ, sigma=1.0):
  
  sigma2 = sigma * sigma
  diff = vectorI - vectorJ
  dotProd = dot(diff,diff)
  
  return exp( -0.5 * dotProd / sigma2 )


def kernelLinear(vectorI, vectorJ, mean): 

  diffI = vectorI - mean
  diffJ = vectorJ - mean
  
  return dot(diffI, diffJ)


def svmTrain(knowns, data, kernelFunc, kernelParams,
             limit=1.0, maxSteps=500, relax=1.3):

  m, n = data.shape
  supports = zeros(m, float)
  change = 1.0 # arbitrary but big start

  kernelArray = zeros((m,m), float)
  for i in range(m):      # xrange in Python 2
    for j in range(i+1):  # xrange in Python 2
      coincidence = kernelFunc(data[i], data[j], *kernelParams)
      kernelArray[i,j] = kernelArray[j,i] = coincidence
  
  kernelArray += 1
  
  steps = 0   
  while (change > 1e-4) and (steps < maxSteps):
    prevSupports = supports.copy()

    sortSup = [(val,i) for i, val in enumerate(supports)]
    sortSup.sort(reverse=True)

    #random.shuffle(sortSup) - also possible

    for support, i in sortSup:
      pull = sum( supports * kernelArray[i,:] * knowns )

      adjust = knowns[i] * pull - 1.0
      supports[i] -= adjust * relax / kernelArray[i,i]
      supports[i] = max(0.0, min(limit, supports[i]))

    nonZeroSup = [(val,i) for i, val in enumerate(supports) if val > 0]
    
    if not nonZeroSup:
      continue

    nonZeroSup.sort()

    inds = [x[1] for x in nonZeroSup]
    niter = 1 + int(sqrt(len(inds)))

    for i in range(niter):  # xrange in Python 2
      for j in inds:
        pull = sum(kernelArray[j,inds] * knowns[inds] * supports[inds])          
        adjust = knowns[j] * pull - 1.0
        supports[j] -= adjust * relax / kernelArray[j,j]
        supports[j] = max(0.0, min(limit, supports[j]))

    diff = supports - prevSupports
    change = sqrt( sum(diff * diff) )
    steps += 1

  return supports, steps, kernelArray


def svmPredict(query, data, knowns, supports, kernelFunc, kernelParams):

  prediction = 0.0
  for j, vector in enumerate(data):
    support = supports[j]
    
    if support > 0:
      coincidence = kernelFunc(vector, query, *kernelParams) + 1.0
      prediction += coincidence * support * knowns[j] 

  return prediction


def svmSeparation(knowns, supports, kernelArray):

  score = 0.0
  nz = [i for i, val in enumerate(supports) if val > 0]

  for i, known in enumerate(knowns):
    prediction = sum(supports[nz] * knowns[nz] * kernelArray[nz, i] )
  
    if known * prediction > 0.0: # same sign
      score += 1.0

  return 100.0 * score / len(knowns)



if __name__ == '__main__':


  print("\nK-nearest neghbour\n")

  knownClasses = [((1.0, 0.0, 0.0), 'warm'), # red
                  ((0.0, 1.0, 0.0), 'cool'), # green
                  ((0.0, 0.0, 1.0), 'cool'), # blue
                  ((0.0, 1.0, 1.0), 'cool'), # cyan
                  ((1.0, 1.0, 0.0), 'warm'), # yellow
                  ((1.0, 0.0, 1.0), 'warm'), # magenta
                  ((0.0, 0.0, 0.0), 'cool'), # black
                  ((0.5, 0.5, 0.5), 'cool'), # grey
                  ((1.0, 1.0, 1.0), 'cool'), # white
                  ((1.0, 1.0, 0.5), 'warm'), # light yellow
                  ((0.5, 0.0, 0.0), 'warm'), # maroon
                  ((1.0, 0.5, 0.5), 'warm'), # pink
                  ]

  result = kNearestNeighbour(knownClasses, (0.7,0.7,0.2), k=3)

  print('Colour class:', result)


  print("\nSelf-organising map\n")

  spread = numpy.array([[0.0, 0.10, 0.2, 0.10, 0.0],
                        [0.1, 0.35, 0.5, 0.35, 0.1],
                        [0.2, 0.50, 1.0, 0.50, 0.2],
                        [0.1, 0.35, 0.5, 0.35, 0.1],
                        [0.0, 0.10, 0.2, 0.10, 0.0]])
 
  rows, cols = 20, 20
  testInput = numpy.random.rand(rows * cols, 3)


  som = selfOrganisingMap(testInput, spread, (rows, cols), 100)

  colors = som*255
  colors = colors.astype(numpy.uint8)
  img1 = Image.fromarray(colors, 'RGB')
  img1.save('som.png', 'PNG')
  img1.show()


  print("\nFeed-forward neural network simple test\n")

  data = [[[0,0], [0]],
          [[0,1], [1]],
          [[1,0], [1]],
          [[1,1], [0]]]

  wMatrixIn, wMatrixOut = neuralNetTrain(data, 2, 1000)

  for inputs, knownOut in data:
    sIn, sHid, sOut =  neuralNetPredict(array(inputs), wMatrixIn, wMatrixOut)
    print(knownOut, sOut[0])

  
  print("\nFeed-forward neural network sequence training\n")

  seqSecStrucData = [('ADTLL','E'),
                     ('DTLLI','E'),
                     ('TLLIL','E'),
                     ('LLILG','E'),
                     ('LILGD','E'),
                     ('ILGDS','E'),
                     ('LGDSL','C'),
                     ('GDSLS','H'),
                     ('DSLSA','H'),
                     ('SLSAG','H'),
                     ('LSAGY','H'),
                     ('SAGYR','C'),
                     ('AGYRM','C'),
                     ('GYRMS','C'),
                     ('YRMSA','C'),
                     ('RMSAS','C')]


  aminoAcids = 'ACDEFGHIKLMNPQRSTVWY'
  aaIndexDict = {}
  for i, aa in enumerate(aminoAcids):
    aaIndexDict[aa] = i


  ssIndexDict = {}
  ssCodes = 'HCE'
  for i, code in enumerate(ssCodes):
    ssIndexDict[code] = i


  trainingData = []
  for seq, ss in seqSecStrucData:
 
    inputVec = convertSeqToVector(seq, aaIndexDict)
    outputVec = convertSeqToVector(ss, ssIndexDict)
 
    trainingData.append( (inputVec, outputVec) )


  wMatrixIn, wMatrixOut = neuralNetTrain(trainingData, 3, 1000)

  
  print("\nFeed-forward neural network sequence prediction\n")

  testSeq = 'DLLSA'
  testVec = convertSeqToVector(testSeq, aaIndexDict)
  testArray = array( [testVec,] )

  sIn, sHid, sOut =  neuralNetPredict(testArray, wMatrixIn, wMatrixOut)
  index = sOut.argmax()
  print("Test prediction: %s" % ssCodes[index])

  
  print("\nSupport vector machine training\n")

  random.seed(int(time()))

  numPoints = 20
  catData = []

  for x in range(1,6):
    for y in range(1,6):
      xNorm = x/6.0      # Normalise range [0,1]
      yNorm = y/6.0

      if (x == 3) and (y == 3):
        category = -1.0
 
      elif (x%2) == (y%2):
        category = 1.0

      else:
        category = -1.0

      xvals = random.normal(xNorm, 0.05, numPoints)
      yvals = random.normal(yNorm, 0.05, numPoints)

      for i in range(numPoints):  # xrange in Python 2
        catData.append( (xvals[i], yvals[i], category) )

  catData = array(catData)
  random.shuffle(catData)

  knowns = catData [:,-1]
  data = catData [:,:-1]

  params = (0.1,)
  supports, steps, kernelArray = svmTrain(knowns, data, kernelGauss, params)

  score = svmSeparation(knowns, supports, kernelArray)
  print('Known data: %5.2f%% correct' % ( score ))
  

  print("\nSupport vector machine prediction boundaries\n")

  ds1x = []
  ds1y = []
  ds2x = []
  ds2y = []

  x = 0.0
  while x < 1.0:

    y = 0.0
    while y < 1.0:
      query = array( (x,y) )

      prediction = svmPredict(query, data, knowns, supports,
                              kernelGauss, params)

      if prediction > 0:
        ds1x.append(x)
        ds1y.append(y)
      else:
        ds2x.append(x)
        ds2y.append(y)
 
      y += 0.02
    x += 0.02

  pyplot.scatter( ds1x, ds1y, color='grey' )
  pyplot.scatter( ds2x, ds2y, color='black' )
  pyplot.show()

