// 의생프 - 파이썬 연동
#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyObject* math_cube(PyObject* self, PyObject* args) {
    double x;
    if (!PyArg_ParseTuple(args, "d", &x)) {
        return NULL;
    }
    return PyFloat_FromDouble(x * x * x);
}

static PyMethodDef MathMethods[] = {
    {"cube", math_cube, METH_VARARGS, "Calculate the cube of a number."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mathmodule = {
    PyModuleDef_HEAD_INIT,
    "mathmodule",
    "A simple math module",
    -1,
    MathMethods
};

PyMODINIT_FUNC PyInit_mathmodule(void) {
    return PyModule_Create(&mathmodule);
}
