import rpy2.robjects as R

def binomialTailTest(count, nTrials, pEvent, oneSided=True):
  
  alt = 'greater' if oneSided else 'two.sided'
  
  func = R.r['binom.test']
  result = func(x=count, n=nTrials, p=pEvent, alternative=alt)
    
  return result[2][0]

count = 530
nTrials = 1000
pEvent = 0.5

result = binomialTailTest(count, nTrials, pEvent, oneSided=True)
print( 'Binomial one tail', result)

result = binomialTailTest(count, nTrials, pEvent, oneSided=False)
print( 'Binomial two tail', result)

def tTest(x, y, sameVariance=False):
  
  func = R.r['t.test']
  argDict = {'var.equal': sameVariance}
  result = func(x=R.FloatVector(x), y=R.FloatVector(y), **argDict)
   
  return result[0][0], result[2][0]
  
from numpy import array

samples1 = array([1.752, 1.818, 1.597, 1.697, 1.644,  1.593])
samples2 = array([1.878, 1.648, 1.819, 1.794, 1.745,  1.827])

print('Same variance result', tTest(samples1, samples2, sameVariance=True))
# Result is: -2.072, 0.0650

print('Not same variance result', tTest(samples1, samples2, sameVariance=False))
# Result is: # -2.072 0.0654


  