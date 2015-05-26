from matplotlib import pyplot

from random import randint, sample

from numpy import array, cov, diag, dot, linalg, ones
from numpy import outer, random, sqrt, vstack, zeros


def euclideanDist(vectorA, vectorB):
  
  diff = vectorA-vectorB
  
  return sqrt(dot(diff,diff))


def findNeighbours(data, distFunc, threshold):
  
  neighbourDict = {}
  
  n = len(data)
  for i in range(n):
    neighbourDict[i] = []

  for i in range(0,n-1):
    for j in range(i+1,n):
      dist = distFunc(data[i], data[j])
      
      if dist < threshold:
        neighbourDict[i].append(j)
        neighbourDict[j].append(i)

  return neighbourDict


def simpleCluster(data, threshold, distFunc=euclideanDist):
  
  neighbourDict = findNeighbours(data, distFunc, threshold)

  clusters = []
  pool = set(range(len(data)))

  while pool:
    i = pool.pop()
    neighbours = neighbourDict[i]
    cluster = set()
    cluster.add(i)
 
    pool2 = set(neighbours)
    while pool2:
      j = pool2.pop()
      
      if j in pool:
        pool.remove(j)
        cluster.add(j)
        neighbours2 = neighbourDict[j]
        pool2.update(neighbours2)
 
    clusters.append(cluster)

  clusterData = []
  for cluster in clusters:
    clusterData.append( [data[i] for i in cluster] )
  
  return clusterData


def dbScanCluster(data, threshold, minNeighbour, distFunc=euclideanDist):     
  
  neighbourDict = findNeighbours(data, distFunc, threshold)

  clusters = []
  noise = set()
  pool = set(range(len(data)))

  while pool:
    i = pool.pop()
    neighbours = neighbourDict[i]
    
    if len(neighbours) < minNeighbour:
      noise.add(i)
    
    else:
      cluster = set()
      cluster.add(i)

      pool2 = set(neighbours)
      while pool2:
        j = pool2.pop()
        
        if j in pool:
          pool.remove(j)
          neighbours2 = neighbourDict.get(j, [])
 
          if len(neighbours2) < minNeighbour:
            noise.add(j)
          
          else:  
            pool2.update(neighbours2)
            cluster.add(j)
 
      clusters.append(cluster)
 
  noiseData = [data[i] for i in noise]
  
  clusterData = []
  for cluster in clusters:
    clusterData.append( [data[i] for i in cluster] )
      
  return clusterData, noiseData
  
  
def kMeans(data, k, centers=None):
  
  if centers is None:
    centers = array( sample(list(data), k) )  # list() not needed in Python 2

  change = 1.0
  prev = []

  while change > 1e-8:

    clusters = [[] for x in range(k)]
    for vector in data:
      diffs = centers - vector
      dists = (diffs * diffs).sum(axis=1)
      closest = dists.argmin()
      clusters[closest].append(vector)
     
    change = 0
    for i, cluster in enumerate(clusters):
      cluster = array(cluster)
      center = cluster.sum(axis=0)/len(cluster)
      diff = center - centers[i]
      change += (diff * diff).sum()
      centers[i] = center
    
  return centers, clusters


def kMeansSpread(data, k):

  n = len(data)
  index = randint(0, n-1)
  indices = set([index])
  
  influence = zeros(n)
  while len(indices) < k:
    diff = data - data[index]
    sumSq = (diff * diff).sum(axis=1) + 1.0
    influence += 1.0 / sumSq
    index = influence.argmin()
    
    while index in indices:
      index = randint(0, n-1)
    
    indices.add(index)    
  
  centers = vstack([data[i] for i in indices])
    
  return kMeans(data, k, centers)


def jumpMethodCluster(data, kRange=None, cycles=10):

  n, dims = data.shape

  if kRange is None:
    start, limit = (2, n+1)
  else:
    start, limit = kRange

  power = dims/2.0
  distortions = {}
  invCovMat = linalg.pinv(cov(data.T))
  
  for k in range(start, limit):
    meanDists = zeros(cycles)

    for c in range(cycles):
      sumDist = 0.0
      centers, clusters = kMeans(data, k)

      for i, cluster in enumerate(clusters):
        size = len(cluster)
        diffs = array(cluster) - centers[i]
 
        for j, diff in enumerate(diffs):
          dist = dot(diff.T, dot(diff, invCovMat))
          sumDist += dist / size

      meanDists[c] = sumDist / (dims * k)
    
    distortions[k] = min(meanDists) ** (-power)
  
  maxJump = None
  bestK = None
  
  for k in range(start+1, limit):
    jump = distortions[k] - distortions[k-1]
    
    if (maxJump is None) or (jump > maxJump):
      maxJump = jump
      bestK = k

  return bestK


def principleComponentAnalysis(data, n=2):

  samples, features = data.shape

  meanVec = data.mean(axis=0)
  dataC = (data - meanVec).T

  covar = cov(dataC)
  evals, evecs = linalg.eig(covar)

  indices = evals.argsort()[::-1]

  evecs = evecs[:,indices]

  basis = evecs[:,:n]
  energy = evals[:n].sum()

  # norm wrt to variance
  #sd = sqrt(diag(covar))
  #zscores = dataC.T / sd

  return basis, energy


def extractPrincipleComponent(data, precision=1e-9):

  samples, features = data.shape
  meanVec = data.mean(axis=0)
  dataC = data - meanVec
  
  pc1 = random.random(features)
  pc0 = pc1 - 1.0
   
  while abs((pc0-pc1).sum()) > precision:

    t = zeros(features)    
    for datum in dataC:
      t += dot(datum, pc1) * datum

    pc0 = pc1
    pc1 = t / sqrt(dot(t,t))
    
  return pc1


def twoClassLda(dataA, dataB):

  meanA = dataA.mean(axis=0)
  meanB = dataB.mean(axis=0)
    
  covA  = cov(dataA.T) 
  covB  = cov(dataB.T)
  
  nA = len(dataA)-1.0
  nB = len(dataB)-1.0
  
  scatterWithin = nA * covA + nB * covB
  scatterBetween = meanA - meanB

  discrim = dot(linalg.inv(scatterWithin),scatterBetween)

  transfA = dot(dataA, discrim.T)
  transfB = dot(dataB, discrim.T)
  
  divide = dot(discrim,(meanA+meanB))/2.0
  
  return transfA, transfB, divide


if __name__ == '__main__':


  print("\nSimple associative clustering\n")

  spread = 0.12
  sizeDims = (100,2)
  data = [random.normal(( 0.0, 0.0), spread, sizeDims),
          random.normal(( 1.0, 1.0), spread, sizeDims),
          random.normal(( 1.0, 0.0), spread, sizeDims)]

  data = vstack(data)
  random.shuffle(data) # Randomise order

  clusters = simpleCluster(data, 0.10)

  colors = ['#F0F0F0','#A0A0A0','#505050',
            '#D0D0D0','#808080','#202020']

  markers = ['d','o','s','>','^']

  i = 0
  for cluster in clusters:
     allX, allY = zip(*cluster)

     if len(cluster) > 3:
       color = colors[i % len(colors)]
       marker = markers[i % len(markers)]
       pyplot.scatter(allX, allY, s=30, c=color, marker=marker)
       i += 1
 
     else:
       pyplot.scatter(allX, allY, s=5, c='black', marker='o')
 
  pyplot.show()


  print("\nDensity based associative clustering\n")

  clusters, noise = dbScanCluster(data, 0.10, 2)

  i = 0
  for cluster in clusters:
     allX, allY = zip(*cluster)

     if len(cluster) > 3:
       color = colors[i % len(colors)]
       marker = markers[i % len(markers)]
       pyplot.scatter(allX, allY, s=30, c=color, marker=marker)
       i += 1
 
     else:
       pyplot.scatter(allX, allY, s=5, c='black', marker='o')
 
  pyplot.show()


  print("\nK-means clustering\n")

  testDataA = random.random((1000,2)) # No clumps

  centers, clusters = kMeans(testDataA, 3)


  testDataB1 = random.normal(0.0, 2.0, (100,2))
  testDataB2 = random.normal(7.0, 2.0, (100,2))
  testDataB = vstack([testDataB1, testDataB2]) # Two clumps

  centers, clusters = kMeans(testDataB, 2)


  colors = ['#FF0000','#00FF00','#0000FF',
            '#FFFF00','#00FFFF','#FF00FF']

  for i, cluster in enumerate(clusters):
     x, y = zip(*cluster)
     color = colors[i % len(colors)]
     pyplot.scatter(x, y, c=color, marker='o')

  x, y = zip(*centers)
  pyplot.scatter(x, y, s=40, c='black', marker='o')
  pyplot.show()


  print("\nJump method to determine number of k-means clusters\n")

  data = [random.normal(( 0.0, 0.0), spread, sizeDims),
          random.normal(( 1.0, 1.0), spread, sizeDims),
          random.normal(( 1.0, 0.0), spread, sizeDims)]

  data = vstack(data)
  random.shuffle(data)

  k = jumpMethodCluster(data, (2, 10), 20)

  print('Number of clusters:', k)

 

  print("\nPrinciple component extraction\n")

  testData = random.normal(0.0, 2.0, (100,2))

  shear = array([[2,1],[1,0]])

  testData = dot(testData, shear)

  pc1 = extractPrincipleComponent(testData)
  print('Quick PC1:', pc1)
 

  print("\nFull principle component analysis\n")

  basis, energy = principleComponentAnalysis(testData, n=2)
  print('Full PCA:', basis, energy)

  x,y = zip(*testData)
  pyplot.scatter(x, y, s=20, c='#F0F0F0', marker='o')

  x,y = zip(-10*pc1, 10*pc1)
  pyplot.plot(x, y)


  transformed = dot(testData, basis)

  x,y = zip(*transformed)
  pyplot.scatter(x, y, s=10, c='#000000', marker='^')

  pyplot.show()

  testData1 = random.normal(0.0, 2.0, (100,2)) + array([-10.0,5.0])
  testData2 = random.normal(0.0, 6.0, (100,2))


  print("\nTwo-class linear discriminant analysis\n")

  x, y = zip(*testData1)
  pyplot.scatter(x, y, s=25, c='#404040', marker='o')

  x, y = zip(*testData2)
  pyplot.scatter(x, y, s=25, c='#FFFFFF', marker='^')

  meanA = testData1.mean(axis=0)
  meanB = testData2.mean(axis=0)
  
  x, y = zip(meanA, meanB)
  pyplot.plot(x, y)
  

  pyplot.show()

  proj1, proj2, div = twoClassLda(testData1, testData2)
  print(div)

  x = proj1
  y = [0.5] * len(x)
  pyplot.scatter(x, y, s=35, c='#404040', marker='o')

  x = proj2
  y = [-0.5] * len(x)
  pyplot.scatter(x, y, s=35, c='#FFFFFF', marker='^')

  pyplot.plot((div, div), (1.0, -1.0))

  pyplot.show()


