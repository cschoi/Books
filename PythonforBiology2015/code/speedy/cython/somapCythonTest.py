from somapCython import selfOrganisingMap

if __name__ == '__main__':

  import numpy
  import time
  #from PIL import Image

  spread = [[0.0, 0.10, 0.2, 0.10, 0.0],
            [0.1, 0.35, 0.5, 0.35, 0.1],
            [0.2, 0.50, 1.0, 0.50, 0.2],
            [0.1, 0.35, 0.5, 0.35, 0.1],
            [0.0, 0.10, 0.2, 0.10, 0.0]]
  spread = numpy.array(spread)

  rows, cols = 100, 100
  depth = 3
  testInput = numpy.random.random(size=(rows*cols, depth))

  nsteps = 100

  t0 = time.time()

  size = (rows, cols)
  somap = selfOrganisingMap(testInput, spread, size, nsteps)

  t1 = time.time()
  print 'time taken = %.3f' % (t1-t0)

  """
  colors = numpy.ndarray(shape=(rows,cols,depth), dtype='uint8')
  for i in range(rows):
    for j in range(cols):
      for k in range(depth):
        colors[i][j][k] = int(255*somap[i][j][k])

  img = Image.fromarray(colors, 'RGB')
  img.save('som.png', 'PNG')
  img.show()
"""

