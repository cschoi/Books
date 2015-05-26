x = abs(-3.0)
print x       # Result: 3.0

def calcDistance(vec1, vec2):

  from math import sqrt

  d2 = 0

  for i in range(len(vec1)):
    delta = vec1[i] - vec2[i]
    d2 += delta * delta

  dist = sqrt(d2)

  return dist

w1 = (23.1, 17.8, -5.6)
w2 = (8.4, 15.9, 7.7)
distance = calcDistance(w1, w2)

distance = calcDistance((23.1, 17.8, -5.6), (8.4, 15.9, 7.7))

calcDistance(w1, w2)

def getSign(value):

  if value > 0:
    return 'Positive'
  elif value < 0:
    return 'Negative'

  return

  print 'Hello world'  # This line is ignored

print getSign(33.6) # Result: 'Positive'
print getSign(-7)   # Result: 'Negative'
print getSign(0)    # Result: None


def myFunction(value1, value2):
  
  total = value1 + value2
  difference = value1 - value2
  product = value1 * value2

  return total, difference, product

values = myFunction(3, 7) # Grab output as a whole tuple

print values  # Result: (10, -4, 21)

x, y, z = myFunction(3, 7) # Grab individual values

print x  # Result: 10
print y  # Result: -4
print z  # Result: 21



def runSimulation(numberSteps=1000):
  pass

runSimulation(500)

runSimulation()

def myFunction(parameters=[]):

  parameters.append(100)
  print parameters

myFunction() # Result: [100]
myFunction() # Result: [100, 100]
myFunction() # Result: [100, 100, 100]


def myFunction(parameters=None):

  if parameters is None:
    parameters = []

  parameters.append(100)
  print parameters

myFunction() # Result: [100]
myFunction() # Result: [100]

def runSimulation(initialTemperature, numberSteps=1000):
  pass

def runSimulation(numberSteps=1000, initialTemperature=300.0):
  pass

runSimulation(500, 400.0)

runSimulation(500)

runSimulation()

runSimulation(initialTemperature=400.0)

runSimulation(500, initialTemperature=400.0)

runSimulation(numberSteps=500, initialTemperature=400.0)

runSimulation(initialTemperature=400.0, numberSteps=500)

runSimulation(500, 400.0)


def reverseComplement(sequence, isDna=True):
    
  from string import maketrans

  if isDna:
    sequence = sequence.replace('U','T')
    transTable = maketrans('ATGC', 'TACG')

  else:
    sequence = sequence.replace('T','U')
    transTable = maketrans('AUGC', 'UACG')

  complement = sequence.translate(transTable)
  reverseComp = complement[::-1]

  return reverseComp


seq1 = 'GATTACA'
seq2 = "AUGGUG"

print  reverseComplement(seq1)			# TGTAATC
print  reverseComplement(seq1, isDna=False)	# UGUAAUC
print  reverseComplement(seq2, False)		# CACCAU

def draw(molecule, *args, **kw):
  pass

myMolecule = None

draw(myMolecule, 1, 100, color="red", speed=1.0)

draw(myMolecule, color="red", speed=1.0)

draw(myMolecule)

def draw(molecule, *irrelevantArgs, **irrelevantKw):
  pass

def setStyle(*args, **kw):
  pass

def setStyle(fgColor="red", bgColor="black", linestyle="plain"):
  pass

def draw(molecule, *args, **kw):
  # some code
  setStyle(*args, **kw)
  # some more code
 
draw(myMolecule, bgColor="green")

(value1, value2, value3) = 1,2,3

tupleArgs = (value1, value2, value3)

dictKw = {'color':'blue', 'depth':3, 'gamma':0.271728}

def someFunc(*args, **kw):
  pass
 
someFunc(*tupleArgs, **dictKw)


def mathFunction(x, y):
  z = (x+y)*(x-y)
  return z

answer = mathFunction(4, 7)

print answer # Fine

print z # Fails; z does not exist outside function definition

KILOS_PER_POUND = 0.45359237

def poundsToKilos(pounds):

  kilos = pounds *  KILOS_PER_POUND
  return  kilos

def kilosToPounds(kilos):

  pounds = kilos /  KILOS_PER_POUND
  return  pounds


counter = 0

def someFunction(argument):
  
  global counter
  counter += 1
  performOperation(argument)


counter = 0
performOperation = str

def someFunction(argument, counterVal):

  counterVal += 1
  performOperation(argument)    

  return counterVal

counter = someFunction(input, counter)

def drawAtoms(atoms):

  def compareAtoms(atom1, atom2):
    return cmp(atom1.z, atom2.z)

  atoms.sort(compareAtoms)

cube = lambda x: x*x*x
 
print cube(3)  # Result: 27

def cube(x):
  return x*x*x

print cube(3)  # Result: 27

def jobFunc(arg1, errorFunc):
  pass
  
"""
def jobFunc(arg1, errorFunc('Warning', color='Red')): # Wrong
  pass
"""
lmb = lambda: errorFunc('Warning', color='Red')

def jobFunc(arg1, lmb): 
  pass
  
def drawAtoms(atoms):

  atoms.sort(lambda atom1,atom2: cmp(atom1.z, atom2.z))


def getSumSquares(n):
  result = sum([i*i for i in range(n)])
  return result

def addEmphasis(inputFunc):
  def modifyFunc(*args):
    result = inputFunc(*args)
    result = '****%d****' % result
    return result  
  return modifyFunc


print getSumSquares(127)  # Gives: 674751

getSumSquares =  addEmphasis(getSumSquares) # Redefine

print getSumSquares(127)  # Gives: ****674751****


@addEmphasis
def getSumSquares(n):
  result = sum([i*i for i in range(n)])
  return result


def usingTimer(function):
  from time import time
  
  def timer(*args, **kw):
    start = time()
    output = function(*args, **kw)
    end = time()
    print 'Function call took %f seconds' % (end-start)
    
    return output
  
  return timer

@usingTimer 
@addEmphasis
def getSumSquares(n):
  result = sum([i*i for i in range(n)])
  return result

print getSumSquares(127)

