from numpy cimport ndarray
from numpy import random, exp
import cython

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef void updateSomap(ndarray[double, ndim=3] somap,
                      ndarray[double, ndim=1] inputs,
                      ndarray[double, ndim=2] spread,
                      int nrows, int ncols, int depth,
                      int width, int height, double decay):

  cdef int halfWidth, halfHeight, i, j, k, l, m, imin, jmin
  cdef double alpha, diff, diff2, diff2min
  """below is for NumPy alternative
  cdef ndarray[double, ndim=3] diff
  cdef ndarray[double, ndim=2] diff2
"""

  halfWidth = (width-1) // 2
  halfHeight = (height-1) // 2

  imin = jmin = -1
  diff2min = 0.0
  for i in range(nrows):
    for j in range(ncols):
      diff2 = 0.0
      for k in range(depth):
        diff = somap[i,j,k] - inputs[k]
        diff2 += diff * diff

      if ((imin == -1) or (diff2 < diff2min)):
        imin = i
        jmin = j
        diff2min = diff2
  """below is for NumPy alternative
  diff = somap - inputs[i]
  diff2 = (diff*diff).sum(axis=2)
  cdef int index = dist2.argmin()
  imin = index // nrows
  jmin = index % nrows
"""

  for k in range(width):
    i = (imin + k - halfWidth + nrows) % nrows
    for l in range(height):
      j = (jmin + l - halfHeight + ncols) % ncols
      alpha = decay * spread[k,l]
      for m in range(depth):
        somap[i,j,m] = (1.0-alpha) * somap[i,j,m] + alpha * inputs[m]

def selfOrganisingMap(ndarray[double, ndim=2] inputs,
                      ndarray[double, ndim=2] spread,
                      size, nsteps):
  
  cdef int nrows, ncols, ninputs, depth, width, height
  cdef int step, i, j
  cdef double decay
  cdef ndarray[double, ndim=3] somap

  nrows, ncols = size
  ninputs = len(inputs)
  depth = len(inputs[0])
  width = len(spread)
  height = len(spread[0])

  somap = random.random((nrows, ncols, depth))

  for step in range(nsteps):
    decay = exp(-step / float(nsteps))
    for i in range(ninputs):
      updateSomap(somap, inputs[i], spread, nrows, ncols, depth,
                  width, height, decay)

  return somap

