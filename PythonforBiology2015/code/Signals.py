from Images import imageToPixmapRGB,  pixmapToImage

from matplotlib import pyplot
from PIL import Image
from numpy import arange, zeros, exp, pi, sqrt, array, outer, argwhere
from numpy.random import standard_normal as normal
from scipy.fftpack import fft
from scipy.ndimage.filters import maximum_filter
from scipy import optimize

I = complex(0, 1)  # Square root of -1; 0+1j

def createSignal(parameters, tStep, nPoints, noise):

  sig = zeros(nPoints, dtype=complex)
  t = tStep * arange(nPoints, dtype=float)

  for amplitude, frequency, decay in parameters:
    sig += amplitude * exp(2*pi*I*frequency*t) * exp(-decay*t)

  noise *= sqrt(0.5)
  sig += noise*(normal(nPoints) + I*normal(nPoints))

  return sig


def savePlot(x, y, xlabel, ylabel):

  pyplot.plot(x, y, color='k') # k means BlacK
  pyplot.xlabel(xlabel)
  pyplot.ylabel(ylabel)
  fileName = ylabel.replace(' ', '')
  pyplot.savefig(fileName)
  pyplot.close()


class Peak:

  def __init__(self, position, data, dataHeight=None, 
                 linewidth=None):

    self.position = tuple(position)
    self.data = data

    self.point = tuple([int(round(x)) for x in position])

    if dataHeight is None:
      dataHeight = data[self.point]

    self.dataHeight = dataHeight

    if linewidth is None:
      linewidth = self._calcHalfHeightWidth()
      
    self.linewidth = linewidth

    self.fitAmplitude = None
    self.fitPosition = None
    self.fitLinewidth = None

  def _calcHalfHeightWidth(self):

    dimWidths = []
    
    for dim in range(self.data.ndim):
      posA, posB = self._findHalfPoints(dim)
      width = posB - posA
      dimWidths.append(width)

    return dimWidths

  def _findHalfPoints(self, dim):

    height = abs(self.dataHeight)
    halfHt = 0.5 * height
    data = self.data
    point = self.point
    
    testPoint = list(point)
    posA = posB = point[dim]
    
    prevValue = height
    while posA > 0: # Search backwards
      posA -= 1
      testPoint[dim] = posA
      value =  abs(data[tuple(testPoint)])
      
      if value <= halfHt:
        posA += (halfHt-value)/(prevValue-value)  # linear interpolation
        break 
      
      prevValue = value

    lastPoint = data.shape[dim] - 1   
    prevValue = height
    while posB < lastPoint-1: # Search forwards
      posB += 1 
      testPoint[dim] = posB
      value = abs(data[tuple(testPoint)])
      
      if value <= halfHt:
        posB -= (halfHt-value)/(prevValue-value)
        break     
      
      prevValue = value
    
    return posA, posB


  def fit(self, fitWidth=2):

    region = []
    numPoints = self.data.shape
    
    for dim, point in enumerate(peak.position):
      start = max(point-fitWidth, 0)
      end = min(point+fitWidth+1, numPoints[dim])
      region.append( (start, end) )

    self.fitData = self._getRegionData(region) / self.dataHeight

    amplitudeScale = 1.0
    offset = 0.0
    linewidthScale = 1.0

    ndim = self.data.ndim
    params = [amplitudeScale]
    params += ndim * [offset]
    params += ndim * [linewidthScale]

    fitFunc = lambda params: self._fitFunc(region, params)
    result = optimize.fmin(fitFunc, params, xtol=0.01)

    amplitudeScale = result[0]
    offset = result[1:ndim+1]
    linewidthScale = result[ndim+1:]

    peak.fitAmplitude = float(amplitudeScale * peak.dataHeight)
    peak.fitPosition  = list(peak.position + offset)
    peak.fitLinewidth = list(linewidthScale * peak.linewidth)


  def _getRegionData(self, region):

    slices = tuple([slice(start, end) for start, end in region])
    
    return self.data[slices]


  def _fitFunc(self, region, params):

    ndim = self.data.ndim

    amplitudeScale = params[0]
    offset = params[1:1+ndim]
    linewidthScale = params[1+ndim:]
    sliceData = ndim * [0]

    for dim in range(ndim):
      linewidth = linewidthScale[dim] * self.linewidth[dim]
      testPos = offset[dim] + self.position[dim]
      (start, end) = region[dim]
      
      if linewidth > 0:
        x = array(range(start, end))
        x = (x - testPos) / linewidth
        slice1d = 1.0 / (1.0 + 4.0*x*x)
      else:
        slice1d = zeros(end-start)

      sliceData[dim] = slice1d

    heights = amplitudeScale * self._outerProduct(sliceData)
    diff2 = ((heights-self.fitData)**2).mean()
    return sqrt(diff2)


  def _outerProduct(self, data):

    size = [d.shape[0] for d in data]
    product = data[0]

    for dim in range(1, len(size)):
      product = outer(product, data[dim])

    product = product.reshape(size)

    return product


def findPeaks(data, threshold, size=3, mode='wrap'):

  peaks = []

  if (data.size == 0) or (data.max() < threshold):
    return peaks
  
  boolsVal = data > threshold

  maxFilter = maximum_filter(data, size=size, mode=mode) 
  boolsMax = data == maxFilter

  boolsPeak = boolsVal & boolsMax
  
  indices = argwhere(boolsPeak) # True positional indices

  for position in indices:
    position = tuple(position)
    height = data[position]
    peak = Peak(position, data, height)
    peaks.append(peak)

  return peaks

if __name__ == '__main__':

  sigParams = ((1.0, 0.1, 0.01), # Amplitude, frequency, decay
               (2.5, 0.7, 0.05))

  nPoints = 100
  tStep = 1.0
  noise = 0.5

  sig = createSignal(sigParams, tStep, nPoints, noise)


  times = [i*tStep for i in range(nPoints)]
  savePlot(times, sig, 'time', 'signal')

  freqs = fft(sig)

  freqReal = [f.real for f in freqs]
  savePlot(times, freqReal, 'freq', 'FT real')

  freqImag = [f.imag for f in freqs]
  savePlot(times, freqImag, 'freq', 'FT imag')

  powerSpec = [abs(f)**2 for f in freqs]
  savePlot(times, powerSpec, 'freq', 'FT power')


  img = Image.open('examples/Gel2D.png')
  pixmap = imageToPixmapRGB(img)

  data = pixmap.mean(axis=2)
  data -= data.min()
  data /= data.max()
  data = 1.0 - data


  threshold = 0.3 * data.max()
  peaks = findPeaks(data, threshold, size=7, mode='wrap')

  color = (255.0, 0.0, 0.0)
  for peak in peaks:
    x,y = peak.position
    xIndices = (x-1, x, x, x, x+1)
    yIndices = (y, y-1, y, y+1, y)
    pixmap[xIndices, yIndices,:] = color

  img2 = pixmapToImage(pixmap, mode='RGB')
  img2.show()
  img2.save('PickedGel.png')

  for peak in peaks:
    peak.fit(fitWidth=3)
    print(peak.fitAmplitude, peak.fitAmplitude, peak.fitLinewidth)


