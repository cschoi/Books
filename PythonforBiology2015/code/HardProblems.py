import time
import sys

from numpy import array, random, zeros, sqrt, dot
from math import exp, acos, cos, sin, atan2, radians
from random import shuffle, randint

uniform = random.uniform
normal = random.normal

try:
  from PySide import QtCore, QtGui
  HAVE_QT = True
  
except ImportError:
  print("* * * * PySide module not available for display of chemical structure * * * * ")
  HAVE_QT = False


def monteCarlo(numSteps, testFunc, spread=0.1, nDims=2):

  bestPoint = uniform(-1.0, 1.0, nDims)
  prevPoint = bestPoint

  bestValue = testFunc(bestPoint)
  prevValue = bestValue

  for i in range(numSteps): # could use xrange in Python 2
    testPoint = normal(prevPoint, spread, nDims) 
    value = testFunc(testPoint)

    prob = exp(prevValue-value)

    if prob > uniform():
      prevPoint = testPoint
      prevValue = value
 
      if value < bestValue:
        bestPoint = testPoint
        bestValue = value
 
        coordinates = ', '.join(['%.3f' % v for v in testPoint])
        print('%5d [%s] value:%e' % (i, coordinates, value))
  
  return bestValue, bestPoint
  

def getRouteLength(distanceData, route):

  distance = 0.0

  for i, pointA in enumerate(route[:-1]):
    pointB = route[i+1]
    key = frozenset((pointA, pointB))  
    distance += distanceData[key]

  return distance


def travellingSalesman(distanceData, cities, numSteps=10000):
  
  n = len(cities)
  bestRoute = cities[:]
  shuffle(bestRoute)
  
  dists = list(distanceData.values()) # list() not needed in Python 2
  scale = 0.5 * array(dists).std()

  bestDistance = getRouteLength(distanceData, bestRoute)
  prevRoute = bestRoute
  prevDistance = bestDistance
  
  for i in range(numSteps): # could use xrange in Python 2
    
    a = randint(0, n-1)
    b = randint(0, n-1)
    
    route = prevRoute[:]
    route[a] = prevRoute[b]
    route[b] = prevRoute[a]
    
    distance = getRouteLength(distanceData, route)
    score = exp((prevDistance-distance)/scale)

    if score > uniform():
      prevRoute = route
      prevDistance = distance

      if distance < bestDistance:
        bestRoute = route[:]
        bestDistance = distance
        print('%5d %.5f' % (i, distance))
  
  return bestDistance, bestRoute


def calcCityDistances(coordDict):

  cities = list(coordDict.keys()) # list() not needed in Python 2
  n = len(cities)
  distances = {}

  for i in range(n-1):
    cityA = cities[i]
    latA, longA = coordDict[cityA]
    latA = radians(latA)
    longA = radians(longA)    

    for j in range(i+1, n):
      cityB = cities[j]
      latB, longB = coordDict[cityB]
      latB = radians(latB)
      longB = radians(longB)
    
      dLong = abs(longA - longB)
      angle = acos(sin(latA)*sin(latB) + cos(latA)*cos(latB)*cos(dLong))
      dist = angle * 6371.1 # Mean Earth radius (km)
    
      key = frozenset((cityA, cityB)) 
      distances[key] = dist

  return distances


def travellingSalesmanSimAnneal(distanceData, cities, numIter=10000):
  
  n = len(cities)
  bestRoute = cities[:]
  shuffle(bestRoute)
  
  dists = list(distanceData.values()) # list() not needed in Python 2
  scale = 0.5 * array(dists).std()

  bestDistance = getRouteLength(distanceData, bestRoute)
  prevRoute = bestRoute
  prevDistance = bestDistance

  m = float(numIter)             # Use to calculate cool
  for i in range(numIter):       # could use xrange in Python 2
    cool = exp(-i/m)             # The 'temperature' the annealing schedule

    a = randint(0, n-1)
    b = randint(0, n-1)
    
    route = prevRoute[:]
    route[a] = prevRoute[b]
    route[b] = prevRoute[a]
    
    distance = getRouteLength(distanceData, route)
    score = exp( (prevDistance - distance) / (scale*cool) ) # Adjusted score
    
    if score > uniform():
      prevRoute = route
      prevDistance = distance

      if distance < bestDistance:
        bestRoute = route[:]
        bestDistance = distance

        print('%5d Dist:%.5f Temp:%.5f' % (i, distance, cool))
  
  return bestDistance, bestRoute


def simAnneal(numIter, testFunc, spread=0.1, nDims=2):

  n = float(numIter)
  bestPoint = uniform(-1.0, 1.0, nDims)
  bestValue = testFunc(bestPoint)
  prevPoint = bestPoint
  prevValue = bestValue
  
  for i in range(numIter): # could use xrange in Python 2
    cool = exp(-i/n)
        
    testPoint = normal(prevPoint, spread, nDims)
 
    value = testFunc(testPoint)

    prob = exp( (prevValue-value) / cool )  # Adjusted acceptance score.
 
    if prob > uniform():
      prevPoint = testPoint
      prevValue = value
 
    if value < bestValue:
      bestPoint = testPoint
      bestValue = value
      
      pointStr = ' '.join(['%.3f' % p for p in testPoint])
      print('%5d T:%.3f %s value:%e' % (i, cool, pointStr, value))

  return bestValue, bestPoint


def chemParticleDynamics(bondDict, numSteps=5000, bondLen=1.0,
                         timeStep=0.01, updateFunc=None):

  atoms = list(bondDict.keys())  # list() not needed in Python 2
  numAtoms = len(atoms)
  atomCoords = uniform(-10.0, 10.0, (numAtoms, 3))

  indices = range(numAtoms)
  n = float(numSteps)

  for step in range(numSteps): # could use xrange in Python 2
    temp = exp(-step/n)
    
    if updateFunc: # Extra for graphical display
      print("Step:", step)
      updateFunc(atomCoords)
    
    for i in indices[1:]:
      atom = atoms[i]
      coords = atomCoords[i]
      velocity = zeros(3, float)
 
      for j in indices:
        if i == j:
          continue

        delta = coords - atomCoords[j]
        delta2 = delta * delta
        dist2 = delta2.sum()

        bound = bondDict[atoms[j]]
        if atom in bound:
          force = bondLen - sqrt(dist2)
 
        else:
          force = 1.0 / (dist2*dist2)
       
        force = min(max(-200.0, force), 200.0)
        velocity += delta * force * temp * timeStep
       
      atomCoords[i] += velocity
     
  center = atomCoords.mean(axis=0)
  atomCoords = atomCoords-center
  
  return atomCoords

# The following code is extra to the book, to display chemical structures
# using Qt graphics

def getRotationMatrix(axis, angle):

  vLen = sqrt( sum([xyz*xyz for xyz in axis]) )

  x, y, z = [xyz/vLen for xyz in axis]

  c = cos(angle)
  d = 1-c
  s = sin(angle)

  R = [[c+d*x*x,   d*x*y-s*z, d*x*z+s*y],
       [d*y*x+s*z, c+d*y*y,   d*y*z-s*x],
       [d*z*x-s*y, d*z*y+s*x, c+d*z*z  ]]

  return R
    
if HAVE_QT:
  class ChemView(QtGui.QWidget):

    def __init__(self, parent, topology):
  
      QtGui.QWidget.__init__(self, parent)
    
      self.topology = topology
      self.coords = None
      self.depthOfField = 5.0
      self.atomRadius = 30.0
    
      atoms = list(topology.keys())  # list() not needed in Python 2
      self.atoms = atoms
      iDict = dict([(x, i) for i, x in enumerate(atoms)])
    
      bondDict = {}
      for atom in atoms:
        bound = [iDict[a] for a in topology[atom]]
        bondDict[iDict[atom]] = bound
    
      self.bondDict = bondDict
      self.setMinimumHeight(600)
      self.setMinimumWidth(800)
      self.movePos = None

    def mousePressEvent(self, event):
   
      QtGui.QWidget.mousePressEvent(self, event)
      self.movePos = event.pos()
    
    def mouseMoveEvent(self, event):
    
      pos = event.pos()
    
      delta = self.movePos - pos
      dx = delta.x()
      dy = delta.y()
    
      rX = array(getRotationMatrix((0.0, 1.0, 0.0), dx*-0.01))
      rY = array(getRotationMatrix((1.0, 0.0, 0.0), dy*0.01))
    
      c = self.coords
      c = dot(c, rX)
      c = dot(c, rY)
    
      self.movePos = pos
      self.updateCoords(c)
    
    def runDynamics(self):
  
      chemParticleDynamics(self.topology, updateFunc=self.updateCoords)
  
    def updateCoords(self, coords):
    
      self.coords = coords
      self.repaint() # Calls paintEvent()
    
    def paintEvent(self, event):  
    
      if self.coords is None:
        return
    
      scale = 50.0
    
      painter = QtGui.QPainter()
      painter.begin(self)

      dof = self.depthOfField
      rad = self.atomRadius
    
      cx = self.width() / 2.0
      cy = self.height() / 2.0
    
      p1 = QtGui.QColor(0, 0, 0)
      p2 = QtGui.QColor(64, 64, 64)  
    
      setPen = painter.setPen
    
      painter.setBrush(QtGui.QColor(128, 128, 128, 128))
    
      drawEllipse = painter.drawEllipse
      drawText = painter.drawText
      setPen(p2)
    
      for i in self.bondDict.keys():
        x, y, z =  self.coords[i]
        perpective = dof / (z-dof)
        xView = cx + scale*x * perpective
        yView = cy + scale*y * perpective
      
        for j in self.bondDict[i]:
          x2, y2, z2 =  self.coords[j]
          perpective2 = dof / (z2-dof)
          xView2 = cx + scale*x2 * perpective2
          yView2 = cy + scale*y2 * perpective2
        
          painter.drawLine(xView, yView, xView2, yView2)
    
      sortCoords = [(z, i, x, y) for i, (x, y, z) in enumerate(self.coords)]
      sortCoords.sort()
    
      QPointF = QtCore.QPointF
      for z, i, x, y in sortCoords:
        perpective = dof / (z-dof)
        xView = cx + scale*x * perpective
        yView = cy + scale*y * perpective
        rView = rad * perpective
      
        point = QPointF(xView, yView)
        setPen(p2)
        drawEllipse(point,  rView, rView)
        setPen(p1)
        drawText(point,  self.atoms[i])
        
      painter.end()
    
    
if __name__ == '__main__':

  print("Monte Carlo intergration")
  
  uniform = random.uniform
  numSamples = 100000
  numInside = 0
  for i in range(numSamples): # could use xrange in Python 2
    x, y = uniform(-1.0, 1.0, 2)

    if (x * x) + (y * y) < 1.0:
      numInside += 1

  pi = 4.0 * numInside / float(numSamples)
  print(pi)


  print("\nMonte Carlo function minimisation - uniform")

  def testFunc(point):
    x, y = point
    a = 1.0 - x
    b = y - (x * x)
    return (a * a) + (100 * b * b)


  bestPoint = uniform(-5, 5, 2)
  bestValue = testFunc(bestPoint)

  numSteps = 100000
  for i in range(numSteps):
    point = uniform(-5, 5, 2)
    value = testFunc(point)
 
    if value < bestValue:
      bestPoint = point
      bestValue = value
 
      x, y =  point
      print('%5d x:%.3f y:%.3f value:%.3f' % (i, x, y, value))


  print("\nMonte Carlo function minimisation - Markov Chain")

  normal = random.normal
  mumSteps = 100000
  for i in range(numSteps): # could use xrange in Python 2
    point = normal(bestPoint, 0.5, 2)
    value = testFunc(point)
 
    if value < bestValue:
      bestPoint = point
      bestValue = value
 
      x, y = point
      print('%5d x:%.3f y:%.3f value:%e' % (i, x, y, value))


  print("\nMonte Carlo function minimisation - Metropolis-Hastings")

  monteCarlo(100000, testFunc, 2)
  
  
  print("\nTravelling Salesman solver - Markov Chain swaps")

  cityCoords = {'Paris':(48.856667, 2.350833),
                'Marseille':(43.296386, 5.369954),
                'Lyon':(45.759723, 4.842223),
                'Toulouse':(43.604503, 1.444026),
                'Nice':(43.703393, 7.266274),
                'Strasbourg':(48.584445, 7.748612),
                'Nantes':(47.21806, -1.55278),
                'Bordeaux':(44.838611, -0.578333),
                'Montpellier':(43.61194, 3.87722),
                'Rennes':(48.114722, -1.679444),
                'Lille':(50.637222, 3.063333),
                'Le Havre':(49.498889, 0.121111),
                'Reims':(49.26278, 4.03472),
                'Saint-Etienne':(45.434722, 4.390278),
                'Toulon':(43.125, 5.930556)}

  distances = calcCityDistances(cityCoords)
  cities = list(cityCoords.keys())          # Use all the cities

  dist, route = travellingSalesman(distances, cities, 1000000)
  print('%.3f %s' % (dist, ', '.join(route)))

  
  print("\nTravelling Saleman solver - Simulated annealing Markov Chain swaps")

  distances = calcCityDistances(cityCoords)
  cities = list(cityCoords.keys())  # list() not needed in Python 2

  dist, route = travellingSalesmanSimAnneal(distances, cities, 1000000)
  print('%.3f %s' % (dist, '-'.join(route)))

  
  
  print("\nFunction minimisation by simulated annealing -")

  numSteps = 100000
  
  simAnneal(numSteps, testFunc)
  monteCarlo(numSteps, testFunc)

  print("\nParticle dynamics by simulated annealing")

  chemBonds = {'H1': ['O1',], 'O1': ['H1', 'C1'],  'C1': ['O1', 'C2', 'C6'],
               'C2': ['C1', 'H2', 'C3'], 'H2': ['C2'],  'C3': ['C2', 'H3', 'C4'],
               'H3': ['C3'], 'C4': ['C3', 'N7', 'C5'], 'C5': ['C4', 'H5', 'C6'],
               'H5': ['C5'], 'C6': ['C5', 'H6', 'C1'], 'H6': ['C6'],
               'N7': ['C4', 'H7', 'C8'], 'H7': ['N7'], 'C8': ['N7', 'O8', 'C9'],
               'O8': ['C8'], 'C9': ['C8', 'H9a', 'H9b', 'H9c'], 'H9a': ['C9',],
               'H9b': ['C9',], 'H9c': ['C9',]}


  coords = chemParticleDynamics(chemBonds)
  print(coords)

  
  if HAVE_QT:
  
    print("\nStart Qt based gui to display particle dynamics structure")
 
    app = QtGui.QApplication(['Qt Chemical View'])
 
    window = QtGui.QWidget()
    layout = QtGui.QVBoxLayout()
 
    viewPanel = ChemView(window, chemBonds)
    layout.addWidget(viewPanel)
 
    button = QtGui.QPushButton('Run annealing', window)
    button.clicked.connect(viewPanel.runDynamics)
    layout.addWidget(button)
    
    viewPanel.updateCoords(coords)
 
    window.show()
 
    sys.exit(app.exec_())
