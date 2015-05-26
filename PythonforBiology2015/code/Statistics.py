from matplotlib import pyplot
from numpy import random

mean = 0.0
stdDev = 1.0

for nPoints in (10, 100, 1000, 10000, 100000):
  sample = random.normal(mean, stdDev, nPoints)
  pyplot.hist(sample, bins=20, range=(-4,4), normed=True)
  pyplot.show()

values = [1,2,2,3,2,1,4,2,3,1,0]
counts = [(values.count(val), val) for val in set(values)]
count, mode = max(counts)

print( mode )

from scipy import stats
from numpy import array

valArray = array(values, float)
mode, count = stats.mode(valArray)

print('Mode:', mode[0] ) # Result is 2

def getMedian(values):

  vSorted = sorted(values)
  nValues = len(values)

  if nValues % 2 == 0: # even number
    index = nValues//2
    median = sum(vSorted[index-1:index+1])/2.0
  
  else:
    index = (nValues-1)//2
    median = vSorted[index]

  return median

med = getMedian(values)

from numpy import median

med = median(valArray)
print('Median:', med)   # Result is 2

values = [1,2,2,3,2,1,4,2,3,1,0]
mean = sum(values)/float(len(values))

from numpy import array, mean

valArray = array(values, float)

m = valArray.mean()
# or 
m = mean(valArray)

print( 'Mean', m)    # Result is 1.909

valArray2 = array([[7,9,5],
                   [1,4,3]])

print(  valArray2.mean() )
# All elements - result is 4.8333

print(  valArray2.mean(axis=0) )
# Column means - result is [4.0, 6.5, 4.0]

print(  valArray2.mean(axis=1) )
# Row means - result is [7.0, 2.6667]

values = [1,2,2,3,2,1,4,2,3,1,0]

n = float(len(values))
mean = sum(values)/n
diffs = [v-mean for v in values]
variance = sum([d*d for d in diffs])/(n-1) # Unbiased estimate

from numpy import array
valArray = array(values)
 
variance = valArray.var() # Biased estimate
print('Var 1', variance)  # Result is 1.1736

variance = valArray.var(ddof=1) # Unbiased estimate
print('Var 2', variance)        # Result is 1.2909

from numpy import std, sqrt

stdDev = sqrt(variance)
stdDev = std(valArray)        # Biased estimate   - 1.0833
stdDev = valArray.std(ddof=1) # "Unbiased" estimate - 1.1362

print('Std:', stdDev)

stdErrMean = valArray.std(ddof=1)/sqrt(len(valArray))

from scipy.stats import sem

stdErrMean = sem(valArray, ddof=1) # Result is 0.3426

from scipy.stats import skew
from numpy import random

samples = random.gamma(3.0, 2.0, 100) # Example data

skewness = skew(samples)
print( 'Skew', skewness )      # Result is 1.0537

from scipy.stats import binom_test

count, nTrials, pEvent = 530, 1000, 0.5
result = binom_test(count, nTrials, pEvent)
print( 'Binomial two tail', result)

from scipy.stats import binom
from numpy import array, zeros

def binomialTailTest(counts, nTrials, pEvent, oneSided=True):
  
  counts = array(counts)
  
  mean = nTrials * pEvent
  
  if oneSided:
    result = zeros(counts.shape)
    isAboveMean = counts > mean
    aboveIdx = isAboveMean.nonzero()
    belowIdx = (~isAboveMean).nonzero()
    result[aboveIdx] = binom.sf(counts[aboveIdx]-1, nTrials, pEvent)
    result[belowIdx] = binom.cdf(counts[belowIdx], nTrials, pEvent)
    
  else:
    diffs = abs(counts-mean)
    result = binom.cdf(mean-diffs, nTrials, pEvent)
    result += binom.sf(mean+diffs-1, nTrials, pEvent)
    
  return result

counts = [530]
result = binomialTailTest(counts, 1000, 0.5, oneSided=True)
print( 'Binomial one tail', result)

result = binomialTailTest(counts, 1000, 0.5, oneSided=False)
print( 'Binomial two tail', result)

from scipy.stats import poisson
from numpy import abs, array

def poissonTailTest(counts, eventRate, oneSided=True):
     
  counts = array(counts)
  diffs = abs(counts-eventRate)
  result = poisson.cdf(eventRate-diffs, eventRate)
  if not oneSided:
    result *= 2

  return result

counts = [2300, 2400, 2550]
result = poissonTailTest(counts, 2500, oneSided=False)
print( 'Poisson two tail', result)
# Result is [0.00005310, 0.045492, 0.3222]

result = poissonTailTest(counts, 2500, oneSided=True)
print( 'Poisson one tail', result)

from matplotlib import pyplot
from scipy.stats import norm
from numpy import arange

mean = 1.76
stdDev = 0.075
stepSize = 0.0001
normRandVar = norm(mean, stdDev)  # Random variable object

xVals = arange(1.42, 2.1, stepSize)  # Graph range
yVals = normRandVar.pdf(xVals)       # Note PDF not PMF

pyplot.plot(xVals, yVals, color='black')

xVals = arange(mean-2*stdDev, mean+2*stdDev, stepSize)
yVals = normRandVar.pdf(xVals)
pyplot.fill_between(xVals, yVals, 0, color='lightgrey')

xVals = arange(mean-stdDev, mean+stdDev, stepSize) 
yVals = normRandVar.pdf(xVals)
pyplot.fill_between(xVals, yVals, 0, color='grey')

pyplot.show()

from scipy.stats import norm

def normalTailTest(values, meanVal, stdDev, oneSided=True):
  
  normRandVar = norm(meanVal, stdDev)
  
  diffs = abs(values-meanVal)
  result = normRandVar.cdf(meanVal-diffs) # Distrib is symmetric
  if not oneSided:
    result *= 2

  return result

mean = 1.76
stdDev = 0.075
values = array([1.8, 1.9, 2.0])

result = normalTailTest(values, mean, stdDev, oneSided=True)
print( 'Normal one tail', result)
# Result is: [0.297, 0.03097, 0.000687]

from numpy import abs

mean = 1.76
stdDev = 0.075
values = array([1.8, 1.9, 2.0])
zScores = abs(values - mean)/stdDev

print('Z scores', zScores)  

from scipy.stats import zscore, norm

samples = norm.rvs(mean, stdDev, size=25)  # Values for testing
zScores = zscore(samples, ddof=1)          # Unbiased estimators

print('Est. Z scores ', zScores) 

from numpy import sqrt
from scipy.special import erf

def zTestMean(sMean, nSamples, normMean, stdDev, oneSided=True):
  
  zScore = abs(sMean - normMean) / (stdDev / sqrt(nSamples))
  prob = 1-erf(zScore/sqrt(2))
  
  if oneSided:
    prob *= 0.5

  return prob

samples = array([1.752, 1.818, 1.597, 1.697, 1.644, 1.593,
                 1.878, 1.648, 1.819, 1.794, 1.745, 1.827])
mean = 1.76
stDev = 0.075
result = zTestMean(samples.mean(), len(samples),
                   mean, stdDev, oneSided=True)

print( 'Z-test', result) # Result is 0.1179

result = zTestMean(0.59, 100, 0.61, 0.1)

from scipy.stats import ttest_1samp

trueMean = 1.76
samples = array([1.752, 1.818, 1.597, 1.697, 1.644,  1.593,
                 1.878, 1.648, 1.819, 1.794, 1.745,  1.827])

tStat, twoTailProb = ttest_1samp(samples, trueMean)
# Result is: -0.918, 0.378

from scipy.stats import ttest_ind

samples1 = array([1.752, 1.818, 1.597, 1.697, 1.644,  1.593])
samples2 = array([1.878, 1.648, 1.819, 1.794, 1.745,  1.827])

tStat, twoTailProb = ttest_ind(samples1, samples2)
# Result is: -2.072, 0.0650

tStat, twoTailProb = ttest_ind(samples1, samples2, equal_var=False)
# Result is: # -2.072 0.0654  

from numpy import mean, std, sqrt
from scipy.stats import t
 
def tConfInterval(samples, confidence, isOneSided=True):

  n = len(samples)
  sampleMean = mean(samples)
  sampleStdDev = std(samples, ddof=1) # Unbiased estimate

  if not isOneSided:
    confidence = 0.5 * (1+confidence)

  interval = t(n-1).ppf(confidence) * sampleStdDev / sqrt(n)

  return sampleMean, interval

from numpy import array

samples = array([1.752, 1.818, 1.597, 1.697, 1.644,  1.593,
                 1.878, 1.648, 1.819, 1.794, 1.745,  1.827])

sMean, intvl = tConfInterval(samples, 0.95, isOneSided=False)

print('Sample mean: %.3f, 95%% interval:%.4f' % (sMean, intvl))

from scipy.stats import chisquare

obs = array([530, 470])
exp = array([500, 500])
chSqStat, pValue = chisquare(obs, exp)

print('DNA Chi-square:', chSqStat, pValue) # 3.6, 0.05778

from scipy.stats import norm
from numpy import array

bins = array([1.65, 1.7, 1.75, 1.8, 1.85,])
obsd = array([ 14, 15, 33, 22, 8,])

mean = 1.76
std = 0.075

nObs = obsd.sum()
expd = norm.pdf(bins, mean, std)
expd *= nObs / expd.sum()

# Expected counts: 9.196, 19.576, 26.720, 23.385, 13.122

chSqStat, pValue = chisquare(obsd, expd)

print('Chi square A', chSqStat, pValue)

from scipy.stats import chi2

degFree = len(obsd)-1
pValue = chi2.sf(chSqStat, degFree) # Result is 0.129

bins = array([1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.90])
obsd = array([   3,   4,  14,   15,   33,  22,    8,    1])

mean = 1.76
std = 0.075

nObs = obsd.sum()
expd = norm.pdf(bins, mean, std)
expd *= nObs/expd.sum()

#Expected: 0.535, 2.769, 9.194, 19.571, 26.713, 23.379, 13.119, 4.720

from numpy import log

g = 2.0 * sum(obsd * log(obsd/expd))
degFree = len(obsd)-1
pValue = chi2.sf(g, degFree)

print('G test', g, pValue)  # Result is: 17.34, 0.015

chSqStat, pValue = chisquare(obsd, expd)
print('Chi-square', chSqStat, pValue) # Result: 21.98, 0.00256

from numpy import random, cov

xVals = random.normal(0.0, 1.0, 100)
yVals1 = random.normal(0.0, 1.0, 100) # Random, independent of xVals

deltas = random.normal(0.0, 0.75, 100)
yVals2 = 0.5 + 2.0 * (xVals - deltas) # Derived from xVals

cov1 = cov(xVals, yVals1)
# The exact value depends on the random numbers
# Cov 1: [[0.848, 0.022]
#         [0.022, 1.048]]

cov2 = cov(xVals, yVals2)
# The exact value depends on the random numbers
# Cov 2: [[0.848, 1.809]
#         [1.809, 5.819]]

from numpy import corrcoef

r1 = corrcoef(xVals, yVals1)[0, 1] # Result is: 0.0231
r2 = corrcoef(xVals, yVals2)[0, 1] # Result is: 0.8145

from numpy import std

cov2 = cov2[0,1] # X-Y element
stdDevX = std(xVals, ddof=1)
stdDevY = std(yVals2, ddof=1)

r2 = cov2 / (stdDevX*stdDevY)

from numpy import sqrt
from scipy.stats import t
 
nVals = range(5, 101)
rVals = []

for n in nVals:
  tVal = t(n-2).ppf(0.975)
  tVal2 = tVal * tVal
  rVal = sqrt(tVal2/(n-2+tVal2))
  rVals.append(rVal)

pyplot.plot(nVals, rVals, color='black')
pyplot.show()

from numpy import cov, var, mean, random

xVals = random.normal(0.0, 1.0, 100)
yVals = 2.0 + -0.7 * xVals + random.normal(0.0, 0.2, 100)

grad = cov(xVals, yVals)/var(xVals, ddof=1)
yInt = mean(yVals) - grad * mean(xVals)
print('LR 1:', grad, yInt) # Result for one run was: -0.711 2.04

from numpy import sqrt
from scipy.stats import t
 
nVals = range(5, 101)
rVals = []

for n in nVals:
  tVal = t(n-2).ppf(0.975)
  tVal2 = tVal * tVal
  rVal = sqrt(tVal2/(n-2+tVal2))
  rVals.append(rVal)

pyplot.plot(nVals, rVals, color='black')
pyplot.show()

from scipy.stats import linregress
from matplotlib import pyplot

grad, yInt, corrCoeff, pValue, stdErr = linregress(xVals, yVals)

print('LR 2:', grad, yInt, corrCoeff, pValue, stdErr)
# Result for one run was: -0.712, 2.04, -0.949, 9.639e-51, 0.0240

xValsFit = [xVals.min(),xVals.max()]
yValsFit = [yInt + x*grad for x in xValsFit]

pyplot.plot(xVals, yVals, 'o')
pyplot.plot(xValsFit, yValsFit)
pyplot.show()
