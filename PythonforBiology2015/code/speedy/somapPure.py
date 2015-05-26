import math
import random

def updateSomap(somap, vector, size, spread, decay):

  nrows, ncols = size
  depth = len(vector)
  width = len(spread)
  height = len(spread[0])

  imin = jmin = -1
  diff2min = 0
  for i in range(nrows):
    for j in range(ncols):
      diff2 = 0.0
      for k in range(depth):
        diff = somap[i][j][k] - vector[k]
        diff2 += diff * diff

      if (imin == -1) or (diff2 < diff2min):
        imin = i
        jmin = j
        diff2min = diff2

  halfWidth = (width-1) // 2
  halfHeight = (height-1) // 2

  for k in range(width):
    i = (imin + k - halfWidth + nrows) % nrows
    for l in range(height):
      j = (jmin + l - halfHeight + ncols) % ncols
      alpha = decay * spread[k][l]
      for m in range(depth):
        somap[i][j][m] = (1.0-alpha) * somap[i][j][m] + alpha * vector[m]

def selfOrganisingMap(inputs, spread, size, steps=1000):

  nrows, ncols = size
  depth = len(inputs[0])

  somap = nrows * [0]
  for i in range(nrows):
    somap[i] = ncols * [0]
    for j in range(ncols):
      somap[i][j] = depth * [0]
      for k in range(depth):
        somap[i][j][k] = random.random()

  for step in range(steps):
    decay = math.exp(-step / float(steps))
    for vector in inputs:
      updateSomap(somap, vector, size, spread, decay)

  return somap

if __name__ == '__main__':

  import numpy
  import time
  #from PIL import Image

  spread = [[0.0, 0.10, 0.2, 0.10, 0.0],
            [0.1, 0.35, 0.5, 0.35, 0.1],
            [0.2, 0.50, 1.0, 0.50, 0.2],
            [0.1, 0.35, 0.5, 0.35, 0.1],
            [0.0, 0.10, 0.2, 0.10, 0.0]]

  rows, cols = size = (50, 50)
  depth = 3
  testInput = (rows * cols) * [0]
  for i in range(rows*cols):
    testInput[i] = depth*[0]
    for j in range(depth):
      testInput[i][j] = random.random()

  nsteps = 10

  t0 = time.time()

  som = selfOrganisingMap(testInput, spread, size, nsteps)

  t1 = time.time()
  print 'time taken = %.3f' % (t1-t0)

  """
  colors = numpy.ndarray(shape=(rows,cols,depth), dtype='uint8')
  for i in range(rows):
    for j in range(cols):
      for k in range(depth):
        colors[i][j][k] = int(255*som[i][j][k])

  img = Image.fromarray(colors, 'RGB')
  img.save('som.png', 'PNG')
  img.show()
"""

