from multiprocessing import Process, Queue, Pool

def calcFunc(n, m):

  print("Running %d %d" % (n,m)) 
  result = sum([x*x for x in range(n) if x % m == 0])
  print("Result %d %d : %d " % (n, m, result))
  
  return result

if __name__ == '__main__':

  job1 = Process(target=calcFunc, args=(8745678, 2))
  job2 = Process(target=calcFunc, args=(2359141, 3))

  job1.start()
  job2.start()

  job1.join()
  job2.join()

def calcFuncWithQ(queue, n, m):

  result = sum([x*x for x in range(n) if x % m == 0])
  queue.put( (n, m, result) )

if __name__ == '__main__':

  queue = Queue()

  job1 = Process(target=calcFuncWithQ, args=(queue, 8745676, 2) )
  job2 = Process(target=calcFuncWithQ, args=(queue, 2359461, 3) )

  job1.start()
  job2.start()

  job1.join()
  job2.join()

  print("Result", queue.get())
  print("Result", queue.get())

  queue.close()

if __name__ == '__main__':

  inputList = [37645, 8374634, 3487584, 191981, 754967, 12345]
  pool = Pool()

  jobs = []
  for value in inputList:
    inputArgs = (value, 2)
    job = pool.apply_async(calcFunc, inputArgs)
    jobs.append(job)

  results = []
  for job in jobs:
    result = job.get()
    results.append(result)

  pool.close()
  pool.join()

  print(results)

from ctypes import cdll

import sys
if sys.platform[:3] == "win":
  libc = cdll.msvcrt
else:
  from ctypes.util import find_library
  fileName = find_library("c") # "c" for C runtime library
  libc = cdll.LoadLibrary(fileName)

print("time = %d" % libc.time(None))

from ctypes import c_double
x = 3.14159
libc.printf(b"x = %.3f\n", c_double(x))

from ctypes import Structure, c_int

class TimeStruct(Structure):
  _fields_ = [ \
    ('tm_sec', c_int),  # seconds
    ('tm_min', c_int),  # minutes
    ('tm_hour', c_int), # hours
    ('tm_mday', c_int), # day of the month
    ('tm_mon', c_int),  # month
    ('tm_year', c_int), # year
    ('tm_wday', c_int), # day of the week
    ('tm_yday', c_int), # day in the year
    ('tm_isdst', c_int) # daylight saving time
  ]

from ctypes import POINTER, c_long, byref
libc.localtime.restype = POINTER(TimeStruct)

t = libc.time(None)
t = c_long(t)
resultPtr = libc.localtime(byref(t))
result = resultPtr[0]
#print("year = %d, month = %d, day = %d, hour = %d, mins = %d, secs = %d" %
print("day = %04d %02d %02d, time = %02d:%02d:%02d" %
      (result.tm_year+1900, result.tm_mon+1, result.tm_mday,
       result.tm_hour, result.tm_min, result.tm_sec))

