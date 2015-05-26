from PIL import Image, ImageEnhance

from matplotlib import pyplot
from numpy import array, dot, dstack, ones, random, uint8, zeros, exp, mgrid, sqrt
from scipy import ndimage 

def pixmapToImage(pixmap, mode='RGB'):
  
  if pixmap.max() > 255:
    pixmap *= 255.0 / pixmap.max()

  pixmap = array(pixmap, uint8)
  img = Image.fromarray(pixmap, mode)
  
  return img
  
  
def imageToPixmapRGB(img):
  
  img2 = img.convert('RGB')
  w, h = img2.size  
  data = img2.getdata()

  pixmap = array(data, float)
  pixmap = pixmap.reshape((h,w,3))
  
  return pixmap


def gammaAdjust(pixmap, gamma=1.0):

  pixmap = array(pixmap, float)/255.0
  pixmap = pixmap ** gamma
  pixmap *= 255
  
  return pixmap


def normalisePixmap(pixmap):
  
  pixmap -= pixmap.min()
  maxVal = pixmap.max()
  
  
  if maxVal > 0:
    pixmap *= 255.0 / maxVal
  
  return pixmap


def clipPixmapValues(pixmap, minimum=0, maximum=255):
  
  pixmap2 = pixmap.copy()
  grey = pixmap2.mean(axis=2)
  
  boolArray = grey < minimum
  indices = boolArray.nonzero()
  pixmap2[indices] = minimum
  
  boolArray = grey > maximum
  indices = boolArray.nonzero()
  pixmap2[indices] = maximum
  
  return pixmap2


def showHistogram(pixmap):

  grey = pixmap.mean(axis=2)
  values = grey.flatten().tolist()
  
  pyplot.hist(values, 256)
  pyplot.show()
  
  
def convolveMatrix2D(pixmap, matrix, mode='reflect'):
  
  matrix = array(matrix)

  if matrix.ndim != 2:
    raise Exception('Convolution matrix must be 2D')

  if pixmap.ndim not in (2,3):
    raise Exception('Pixmap must be 2D or 3D')
  
  if pixmap.ndim == 2:
    pixmap2 = ndimage.convolve(pixmap, matrix, mode=mode)

  else:
    layers = []
    
    for i in range(3):
      layer = ndimage.convolve(pixmap[:,:,i], matrix, mode=mode)
      layers.append(layer)
      
    pixmap2 = dstack(layers)
    
  
  return pixmap2


def sharpenPixmap(pixmap):

  matrix = [[-1,-1,-1],
            [-1, 8,-1],
            [-1,-1,-1]]
  
  grey = pixmap.mean(axis=2)
  
  pixmapEdge = convolveMatrix2D(grey, matrix)
  normalisePixmap(pixmapEdge)
  
  pixmapEdge -= pixmapEdge.mean()
  pixmapEdge = dstack([pixmapEdge, pixmapEdge, pixmapEdge])
  
  pixmapSharp = pixmap + pixmapEdge 
  pixmapSharp = pixmapSharp.clip(0, 255)
  
  return pixmapSharp


def gaussFilter(pixmap, r=2, sigma=1.4):

  x, y = mgrid[-r:r+1, -r:r+1]

  s2 = 2.0 * sigma * sigma
  
  x2 = x * x / s2
  y2 = y * y / s2
  
  matrix = exp( -(x2 + y2))
  matrix /= matrix.sum()
  
  pixmap2 = convolveMatrix2D(pixmap,  matrix)
  
  return pixmap2


def sobelFilter(pixmap):

  matrix = array([[-1, 0, 1],
                  [-2, 0, 2],
                  [-1, 0, 1]])

  grey = pixmap.mean(axis=2)
  edgeX = convolveMatrix2D(grey, matrix)
  edgeY = convolveMatrix2D(grey, matrix.T)
 
  pixmap2 = sqrt(edgeX * edgeX + edgeY * edgeY)
  normalisePixmap(pixmap2) # Put min, max at 0, 255

  return pixmap2 

  
def getNeighbours(point, points):

  i, j = point
  check = [(i-1, j), (i+1, j),
           (i, j-1), (i, j+1)]
  
  neighbours = [p for p in check if p in points]
  
  return neighbours


def brightPixelCluster(pixmap, threshold=60):

  boolArray = pixmap > threshold
  indices = array(boolArray.nonzero()).T
  points = set([tuple(point) for point in indices])

  clusters = []
  pool = set(points)
  clustered = set()

  while pool:
    pointA = pool.pop()
    neighbours = getNeighbours(pointA, points)
    
    cluster = []
    cluster.append(pointA)
    clustered.add(pointA)
    
    pool2 = set(neighbours)
    while pool2:
      pointB = pool2.pop()
      
      if pointB in pool:
        pool.remove(pointB)
        neighbours2 = getNeighbours(pointB, points)
        pool2.update(neighbours2)
        cluster.append(pointB)

    clusters.append(cluster)
  
  return clusters


if __name__ == '__main__':

  # Make an Image object from file, show its atttributes and methods
  img = Image.open('examples/Cells.jpg')
 
  print(img.size)
  print(img.mode)

  img.show()

  img.save('Cells.png', 'PNG')

  img.convert('CMYK') #  Cyan, Magenta, Yellow, blacK
  img.convert('L')    #  Luminance = greyscale

  w, h = img.size
  img = img.crop((10, 10, w-10, h-10))

  img = img.rotate(270)
  img.save('CellsAdj.png', 'PNG')

  img2 = img.resize((w/2, h/2), Image.ANTIALIAS)
  img2.save('CellsHalfSize.png', 'PNG')

  img2 = img.copy()
  img2.thumbnail((50, 50), Image.ANTIALIAS)
  img2.save('CellsThumb.png', 'PNG')
  
  
  
  # Create images from array pixmaps

  h = 200
  w = 200
  d = 3

  pixmap = zeros((h, w, d)) # black
  pixmap = ones((h, w, d)) # white
  pixmap = random.random((h, w, d)) # random colours
  pixmap *= 255
  
  img1 = pixmapToImage(pixmap)
  img1.show()

  size = (h,w)
  redMatrix   = ones(size)
  greenMatrix = ones(size)
  blueMatrix  = zeros(size)

  pixmap = dstack([redMatrix, greenMatrix, blueMatrix])
  pixmap *= 255

  img1 = pixmapToImage(pixmap)
  img1.show()
  
  
  # Pixmap colour transformations
    
  img = Image.open('examples/CellNucleusRedGreen.png')
  pixmap = imageToPixmapRGB(img)
  
  # Alternative scipy loading
  from scipy import ndimage
  pixmap = ndimage.imread('examples/CellNucleusRedGreen.png')
  pixmap = pixmap[:,:,:3]

  transform = array([[1.0, 1.0, 0.0],
                     [0.0, 0.0, 1.0],
                     [0.0, 0.0, 0.0]])

  pixmap2 = dot(pixmap, transform)

  img2 = pixmapToImage(pixmap2)
  img2.show()

  pixmap = array(pixmap, float)
  pixmap /= 255.0
  pixmap = pixmap ** 2
  pixmap = array(255*pixmap, uint8)

  img2 = pixmapToImage(pixmap)
  img2.show()

  # Image adjustments and filters

  img = Image.open('examples/Cells.jpg')

  processObj = ImageEnhance.Contrast(img)
  img2 = processObj.enhance(2.0)
  img2.show()

  processObj = ImageEnhance.Sharpness(img)

  img2Sharp = processObj.enhance(4.0)
  img2Sharp.show()

  imgBlur = processObj.enhance(0.5)
  imgBlur.show()
 
  imgBrighter = ImageEnhance.Brightness(img).enhance(0.5)
  imgDull = ImageEnhance.Color(img).enhance(0.1)


  # Image adjustments and histograms with NumPy

  img = Image.open('examples/Cells.jpg')
  pixmap = imageToPixmapRGB(img)

  showHistogram(pixmap)

  pixmap2 = gammaAdjust(pixmap, 0.5)
  pixmap3 = clipPixmapValues(pixmap2, 0, 145)
  pixmap4 = normalisePixmap(pixmap3)

  pixmapToImage(pixmap4).show()

  # Image filters with NumPy/SciPy

  img = Image.open('examples/Cells.jpg')
  pixmap = imageToPixmapRGB(img)
  
  matrix = [[ 1, 1, 1],
            [ 1, 8, 1],
            [ 1, 1, 1]]

  pixmapBlur = convolveMatrix2D(pixmap, matrix)
  pixmapBlur /= array(matrix).sum()
  
  pixmapToImage(pixmapBlur).show()

  img = Image.open('examples/Cells.jpg')
  imageToPixmapRGB(img)

  pixmap = sharpenPixmap(pixmap)
  pixmapToImage(pixmap).show()

  pixmap = gaussFilter(pixmap)
  pixmapToImage(pixmap).show()

  pixmapGrey = sobelFilter(pixmap)
  pixmapToImage(pixmapGrey, mode='L').show()
  
  
  # Example feature detection - counting cells in an image
  
  img = Image.open('examples/Cells.jpg')

  pixmap =  imageToPixmapRGB(img)
  pixmap2 = gaussFilter(pixmap)
  pixmap2 = sobelFilter(pixmap2)

  normalisePixmap(pixmap2)

  pixmapToImage(pixmap2, mode='L').show()

  clusters = brightPixelCluster(pixmap2)
  sizes = [len(c) for c in clusters]
  pyplot.hist(sizes, 40)
  pyplot.show()

  smallBlobs = []
  mediumBlobs = []
  bigBlobs = []

  for cluster in clusters:
    n = len(cluster)
 
    if n < 80:
      smallBlobs.append(cluster)
    elif n < 320:
      mediumBlobs.append(cluster)
    else:
      bigBlobs.append(cluster)

  print('Found %d small blobs' % len(smallBlobs))
  print('Found %d medium blobs' % len(mediumBlobs))
  print('Found %d big blobs' % len(bigBlobs))

  grey = pixmap.mean(axis=2)
  colorMap = dstack([grey, grey, grey])

  colors = [(255, 0, 0), (255, 255, 0), (0, 0, 255)]
  categories = [smallBlobs, mediumBlobs, bigBlobs]

  for i, blobs in enumerate(categories):
    color = colors[i]
 
    for cluster in blobs:
      x, y = zip(*cluster)

      colorMap[x,y] = color


  Image.fromarray(array(colorMap, uint8), 'RGB').show()

  numCells = len(mediumBlobs)  # initial guess
  cellAreas = [len(blob) for blob in mediumBlobs]
  meanCellArea = sum(cellAreas) / float(numCells)

  for blob in bigBlobs:
    numCells += int( len(blob) // meanCellArea )

  print('Estimated number of cells: %d' %  numCells)
  
