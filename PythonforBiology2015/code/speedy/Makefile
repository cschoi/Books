PYTHON_DIR=/usr

all: somap.so somapTest

somap.so: py_somap.o somap.o randomNumber.o
	cc -fPIC -shared -o $@ py_somap.o somap.o randomNumber.o -lpython2.7 -lm

somapTest: somapTest.o somap.o randomNumber.o
	cc -o $@ somapTest.o somap.o randomNumber.o -lm

py_somap.o: py_somap.c
	cc -fPIC -c py_somap.c -I${PYTHON_DIR}/include/python2.7 -I${PYTHON_DIR}/lib/python2.7/dist-packages/numpy/core/include/numpy

.c.o:; cc -fPIC -c $<

