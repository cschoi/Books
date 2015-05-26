import numpy
import time

#from PIL import Image

import somap

if __name__ == '__main__':

  spread = numpy.array([[0.0, 0.10, 0.2, 0.10, 0.0],
                      [0.1, 0.35, 0.5, 0.35, 0.1],
                      [0.2, 0.50, 1.0, 0.50, 0.2],
                      [0.1, 0.35, 0.5, 0.35, 0.1],
                      [0.0, 0.10, 0.2, 0.10, 0.0]])

  #rows, cols = 100, 100
  rows, cols = 50, 50
  inputs = numpy.random.rand(rows*cols, 3)
 
  nsteps = 10
  
  t0 = time.time()

  result = somap.somap(inputs, spread, rows, cols, nsteps)

  t1 = time.time()
  print('time taken = %.3f' % (t1-t0))

  #colors = result*255
  #colors = colors.astype(numpy.uint8)
  #img = Image.fromarray(colors, 'RGB')
  #img.save('som.png', 'PNG')
  #img.show()

