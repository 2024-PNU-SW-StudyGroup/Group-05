#define PY_SSIZE_T_CLEAN // 크기 형식을 파이썬 기본 형식으로 지정
#include <Python.h>
#include "numpy/npy_math.h"
#include "numpy/arrayobject.h"

#include <cmath> // sqrt

using namespace std;

double euclidean_distance(double *a, double *b, size_t n) { // unsigned int = size_t (음수 아닌 int)
  double rtn = 0;
  for (size_t i = 0; i < n; i++) {
    double d = a[i] - b[i];
    rtn += d * d;  // sum += pow(d, 2); 로 거듭제곱 계산 가능 (실수여도 가능)
  }
  return sqrt(rtn);
}

double manhattan_distance(double *a, double *b, size_t n) {
  double rtn = 0;
  for (size_t i = 0; i < n; i++) {
    rtn += abs(a[i] - b[i]);
  }
  return rtn;
}

static PyObject* calc_distances(PyObject* self, PyObject* args) {
  PyObject *arg1 = NULL;
  PyObject *arg2 = NULL;
  PyArrayObject *arr1;
  PyArrayObject *arr2;


  double *a, *b;
  size_t n;
  double ed, md;
  PyObject* rtn;

  if (!PyArg_ParseTuple(args, "OO", &arg1, &arg2)) return NULL;
  if ((arr1 = (PyArrayObject*)PyArray_FROM_OTF(arg1, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) goto fail;
  if ((arr2 = (PyArrayObject*)PyArray_FROM_OTF(arg2, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY)) == NULL) goto fail;
  
  a = (double *)PyArray_DATA(arr1);
  b = (double *)PyArray_DATA(arr2);
  n = (size_t)PyArray_SIZE(arr1);

  ed = euclidean_distance(a, b, n);
  md = manhattan_distance(a, b, n);

  Py_DECREF(arr1);
  Py_DECREF(arr2);

  rtn = (PyObject *)PyTuple_New(2);
  PyTuple_SetItem(rtn, 0, PyFloat_FromDouble(ed));
  PyTuple_SetItem(rtn, 1, PyFloat_FromDouble(md));

  return rtn;

fail:
  Py_XDECREF(arr1);
  Py_XDECREF(arr2);
  return NULL;
}

static PyMethodDef ModuleMethods[] = {
    {"calc_distances", calc_distances, METH_VARARGS, "Calculate the distances of two vectors."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "npmodule_utils",
    NULL,
    -1,
    ModuleMethods
};

PyMODINIT_FUNC PyInit_utils(void) {
    PyObject* module = PyModule_Create(&moduledef);
    import_array();  // 이게 들어가 있는게 중요!!
    return module;
}