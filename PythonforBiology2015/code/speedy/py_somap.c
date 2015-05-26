#include "somap.h"

#include "Python.h"
#include "arrayobject.h"
#include <time.h>

static PyObject *ErrorObject;   /* locally-raised exception */

#define RETURN_ERROR(message) \
	{ PyErr_SetString(ErrorObject, message); return NULL; }

static void freeArray2(double **array, int m)
{
    int i;

    for (i = 0; i < m; i++)
	free(array[i]);

    free(array);
}

static void freeArray3(double ***array, int m, int n)
{
    int i;

    for (i = 0; i < m; i++)
	freeArray2(array[i], n);

    free(array);
}

static double **allocArray2(int m, int n)
{
    int i;
    double **array;

    array = (double **) malloc(m*sizeof(double *));
    if (array)
    {
	for (i = 0; i < m; i++)
	    array[i] = (double *) malloc(n*sizeof(double));
    }

    return array;
}

static double **copyArray2(PyArrayObject *array_obj)
{
    int i, j, m, n;
    double **array;

    m = PyArray_DIM(array_obj, 0);
    n = PyArray_DIM(array_obj, 1);
    array = allocArray2(m, n);

    if (array)
    {
	for (i = 0; i < m; i++)
	{
	    for (j = 0; j < n; j++)
          	array[i][j] = *((double *) PyArray_GETPTR2(array_obj, i, j));
	}
    }

    return array;
}

static PyObject *somap(PyObject *self, PyObject *args)
{
    int ninputs, nsteps, nrows, ncols, depth, width, height, i, j, k;
    PyArrayObject *inputs_obj, *spread_obj;
    PyObject *somap_obj;
    double **inputs, **spread, ***somap;
    time_t t0, t1;
    npy_intp dims[3];

    if (!PyArg_ParseTuple(args, "O!O!iii", &PyArray_Type, &inputs_obj, &PyArray_Type, &spread_obj, &nrows, &ncols, &nsteps))
	RETURN_ERROR("need five arguments: inputs, spread, nrows, ncols, nsteps");

    if (!PyArray_Check(inputs_obj))
	RETURN_ERROR("inputs needs to be NumPy array");

    if (PyArray_NDIM(inputs_obj) != 2)
	RETURN_ERROR("inputs needs to be NumPy array with ndim 2");

    if (!PyArray_Check(spread_obj))
	RETURN_ERROR("spread needs to be NumPy array");

    if (PyArray_NDIM(spread_obj) != 2)
	RETURN_ERROR("spread needs to be NumPy array with ndim 2");

    if (PyArray_TYPE(inputs_obj) != NPY_DOUBLE)
	RETURN_ERROR("inputs needs to be array of doubles");

    if (PyArray_TYPE(spread_obj) != NPY_DOUBLE)
	RETURN_ERROR("spread needs to be array of doubles");

    ninputs = PyArray_DIM(inputs_obj, 0);
    depth = PyArray_DIM(inputs_obj, 1);

    width = PyArray_DIM(spread_obj, 0);
    height = PyArray_DIM(spread_obj, 1);

    //if (PyArray_AsCArray((PyObject**) &inputs_obj, (void *) &inputs, PyArray_DIMS(inputs_obj), PyArray_NDIM(inputs_obj), PyArray_DescrFromType(NPY_DOUBLE)) < 0)
    if (!(inputs = copyArray2(inputs_obj)))
	RETURN_ERROR("getting inputs as C array");

    //if (PyArray_AsCArray((PyObject**) &spread_obj, (void *) &spread, PyArray_DIMS(spread_obj), PyArray_NDIM(spread_obj), PyArray_DescrFromType(NPY_DOUBLE)) < 0)
    if (!(spread = copyArray2(spread_obj)))
    {
    	//PyArray_Free((PyObject*) inputs_obj, (void *) inputs);
        freeArray2(inputs, ninputs);
	RETURN_ERROR("getting spread as C array");
    }

    t0 = time(NULL);

    somap = selfOrganisingMap(inputs, ninputs, depth, nrows, ncols, spread, width, height, nsteps);

    t1 = time(NULL);

    printf("Time for just somap = %ld\n", t1-t0);

    //somap_obj = PyArray_NewFromDescr(&PyArray_Type, PyArray_DescrFromType(NPY_DOUBLE), PyArray_NDIM(inputs_obj), PyArray_DIMS(inputs_obj), NULL, (void *) somap, 0, NULL);
    // below does not work because data not contiguous
    //somap_obj = PyArray_SimpleNewFromData(PyArray_NDIM(inputs_obj), PyArray_DIMS(inputs_obj), NPY_DOUBLE, (void *) somap);
    dims[0] = nrows;
    dims[1] = ncols;
    dims[2] = depth;
    somap_obj = PyArray_SimpleNew(3, dims, NPY_DOUBLE);
    for (i = 0; i < nrows; i++)
      for (j = 0; j < ncols; j++)
        for (k = 0; k < depth; k++)
          *((double *) PyArray_GETPTR3(somap_obj, i, j, k)) = somap[i][j][k];

    //PyArray_Free((PyObject*) inputs_obj, (void *) inputs);
    //PyArray_Free((PyObject*) spread_obj, (void *) spread);
    freeArray3(somap, nrows, ncols);
    freeArray2(inputs, ninputs);
    freeArray2(spread, width);

    return somap_obj;
}

static char somap_doc[] = "Creates a self organising map";

static struct PyMethodDef Somap_type_methods[] =
{
    { "somap",      (PyCFunction) somap,	METH_VARARGS,	somap_doc },
    { NULL,         NULL,          		0,		NULL }
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    "somap",             /* m_name */
    NULL,                /* m_doc */
    -1,                  /* m_size */
    Somap_type_methods,  /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
};
#endif

#if PY_MAJOR_VERSION >= 3
PyObject *PyInit_somap(void)
#else
void initsomap(void)
#endif
{
    PyObject *module;

    /* create the module and add the functions */

#if PY_MAJOR_VERSION >= 3
    module = PyModule_Create(&module_def);
#else
    module = Py_InitModule("somap", Somap_type_methods);
#endif

    import_array();  /* needed for numpy, otherwise it crashes */

    /* create exception object and add to module */
    ErrorObject = PyErr_NewException("somap.error", NULL, NULL);
    Py_INCREF(ErrorObject);
    PyModule_AddObject(module, "error", ErrorObject);

    /* check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module somap");

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

