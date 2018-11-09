#include <Python.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>


typedef enum
{
  CLIB_MSG_TRACE,
  CLIB_MSG_DEVMESSAGE,
  CLIB_MSG_WARNING,
  CLIB_MSG_ERROR,
  CLIB_MSG_FATAL
} CLIB_MSG_IMPORTANCE;


/*
 * The API version must be changed manually each time the API is
 * changed.
 */
static char clib_api_version[] = "0.0.1";


static CLIB_MSG_IMPORTANCE message_importance_threshold = CLIB_MSG_WARNING;


/* #define REFCOUNTDEBUG */

/***** reference count monitoring *********************************************/

#ifdef REFCOUNTDEBUG

#define REFCOUNTDEBUG_STRINGLENGTH 200


typedef struct tag_refcountdebug_record
{
  struct tag_refcountdebug_record *next;
  int refcount;
  PyObject *python_object;
  char name[REFCOUNTDEBUG_STRINGLENGTH];
  int line;
} REFCOUNTDEBUG_RECORD;


typedef struct
{
  REFCOUNTDEBUG_RECORD *record_list;
} REFCOUNTDEBUG_STATE;


REFCOUNTDEBUG_STATE _refcountdebug_state = {NULL};


static REFCOUNTDEBUG_RECORD *_refcountDebug_findRecord(PyObject *o)
{
  REFCOUNTDEBUG_RECORD *r;

  for (r = _refcountdebug_state.record_list; r; r = r->next)
  {
    if (r->python_object == o)
    {
      return (r);
    }
  }
  return (NULL);
}


static void _refcountDebug_setString(char *str, const char *source)
{
  strncpy(str, source, REFCOUNTDEBUG_STRINGLENGTH - 1);
  str[REFCOUNTDEBUG_STRINGLENGTH - 1] = '\0';
}


static void _refcountDebug_setInfo(REFCOUNTDEBUG_RECORD *r, const char *name, int line)
{
  _refcountDebug_setString(r->name, name);
  r->line = line;
}


static REFCOUNTDEBUG_RECORD *_refcountDebug_newRecord(PyObject *o)
{
  REFCOUNTDEBUG_RECORD *r;

  r = (REFCOUNTDEBUG_RECORD *) malloc(sizeof(REFCOUNTDEBUG_RECORD));
  if (r == NULL)
  {
    return (NULL);
  }
  r->next = _refcountdebug_state.record_list;
  r->python_object = o;
  r->refcount = 0;
  _refcountdebug_state.record_list = r;
  return (r);
}


static void _refcountDebug_findOrNewAndIncrement(PyObject *o, const char *name, int line)
{
  REFCOUNTDEBUG_RECORD *r;

  if (o == NULL)
  {
    return;
  }
  r = _refcountDebug_findRecord(o);
  if (r == NULL)
  {
    r = _refcountDebug_newRecord(o);
  }
  if (r != NULL)
  {
    if (r->refcount == 0)
    {
      _refcountDebug_setInfo(r, name, line);
    }
    r->refcount++;
    /* fprintf(stderr, "_refcountDebug_findOrNewAndIncrement: record for %p: refcount = %d\n", (void *) o, r->refcount); */
  }
  else
  {
    fprintf(stderr, "_refcountDebug_findOrNewAndIncrement: newRecord failed\n");
  }
}


static void _refcountDebug_Py_INCREF(PyObject *o, const char *name, int line)
{
  Py_XINCREF(o);
  _refcountDebug_findOrNewAndIncrement(o, name, line);
}


static void _refcountDebug_checkAndDecrement(PyObject *o, const char *funcname, const char *name, int line)
{
  REFCOUNTDEBUG_RECORD *r = _refcountDebug_findRecord(o);

  if (r == NULL)
  {
    fprintf(stderr, "refcountDebug: %s on unknown object at %p\n", funcname, (void *) o);
    fprintf(stderr, "  name: \"%s\"\n", name);
    fprintf(stderr, "  line: %d\n", line);
  }
  else
  {
    if (r->refcount > 0)
    {
      r->refcount--;
      /* fprintf(stderr, "_refcountDebug_Py_DECREF: decremented refcount of %p to %d\n", (void *) o, r->refcount); */
    }
    else
    {
      fprintf(stderr, "refcountDebug: Py_DECREF beyond 0\n");
      fprintf(stderr, "  name: \"%s\"\n", name);
      fprintf(stderr, "  line: %d\n", line);
      fprintf(stderr, "  record name: \"%s\"\n", r->name);
      fprintf(stderr, "  record line: %d\n", r->line);
    }
  }
}


static void _refcountDebug_Py_DECREF(PyObject *o, const char *name, int line)
{
  if (o != NULL)
  {
    _refcountDebug_checkAndDecrement(o, "Py_DECREF", name, line);
  }
  Py_XDECREF(o);
}


static PyObject *_refcountDebug_PyObject_GetAttrString(PyObject *o, char *attr, const char *name, int line)
{
  PyObject *python_object = PyObject_GetAttrString(o, attr);
  char aname[REFCOUNTDEBUG_STRINGLENGTH];

  if (strlen(name) + strlen(attr) + 1 < REFCOUNTDEBUG_STRINGLENGTH)
  {
    sprintf(aname, "%s.%s", name, attr);
  }
  else
  {
    sprintf(aname, "attribute");
  }
  _refcountDebug_findOrNewAndIncrement(python_object, aname, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyTuple_New(int size, const char *name, int line)
{
  PyObject *python_object = PyTuple_New(size);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyList_New(int size, const char *name, int line)
{
  PyObject *python_object = PyList_New(size);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyDict_New(const char *name, int line)
{
  PyObject *python_object = PyDict_New();

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyInt_FromLong(long v, const char *name, int line)
{
  PyObject *python_object = PyInt_FromLong(v);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyFloat_FromDouble(double v, const char *name, int line)
{
  PyObject *python_object = PyFloat_FromDouble(v);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyInstance_New(PyObject *pyClass, PyObject *args, PyObject *kw, const char *name, int line)
{
  PyObject *python_object = PyInstance_New(pyClass, args, kw);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static PyObject *_refcountDebug_PyObject_Call(PyObject *callable_object, PyObject *args, PyObject *kw, const char *name, int line)
{
  PyObject *python_object = PyObject_Call(callable_object, args, kw);

  _refcountDebug_findOrNewAndIncrement(python_object, name, line);
  return (python_object);
}


static int _refcountDebug_PyTuple_SetItem(PyObject *tuple, int i, PyObject *item, const char *name, int line)
{
  _refcountDebug_checkAndDecrement(item, "PyTuple_SetItem", name, line);
  return (PyTuple_SetItem(tuple, i, item));
}


static int _refcountDebug_PyList_SetItem(PyObject *list, int i, PyObject *item, const char *name, int line)
{
  _refcountDebug_checkAndDecrement(item, "PyList_SetItem", name, line);
  return (PyList_SetItem(list, i, item));
}


static int _refcountDebug_numDisplayableRecords(int showall)
{
  REFCOUNTDEBUG_RECORD *r;
  int n = 0;

  for (r = _refcountdebug_state.record_list; r; r = r->next)
  {
    if ((r->refcount > 0) || showall)
    {
      n++;
    }
  }
  return (n);
}


static void _refcountDebug_report(const char *title, int showall, const char *filename, int line)
{
  REFCOUNTDEBUG_RECORD *r;

  fprintf(stderr, "+----- refcountDebug: %-45s -----+\n", title);
  fprintf(stderr, "| file: %-64s |\n", filename);
  fprintf(stderr, "| line: %-6d                                                           |\n", line);
  fprintf(stderr, "+------------------------------------------------------------------------+\n");
  if (_refcountDebug_numDisplayableRecords(showall) == 0)
  {
  fprintf(stderr, "| no records                                                             |\n");
  }
  else
  {
    fprintf(stderr, "| %-41s  %10s  %6s   %6s |\n", "reference name", "PyObject", "line", "count");
    fprintf(stderr, "+------------------------------------------------------------------------+\n");
    for (r = _refcountdebug_state.record_list; r; r = r->next)
    {
      if ((r->refcount > 0) || showall)
      {
	fprintf(stderr, "| %-41s  %10p  %6d   %6d |\n", r->name, (void *) r->python_object, r->line, r->refcount);
      }
    }
  }
  fprintf(stderr, "+------------------------------------------------------------------------+\n\n");
}


static void refcountDebug_init(void)
{
  REFCOUNTDEBUG_RECORD *r = _refcountdebug_state.record_list, *r1;

  while (r)
  {
    r1 = r;
    r = r->next;
    free(r1);
  }
  _refcountdebug_state.record_list = NULL;
}


/* should also intercept the functions that steal references */

#define PyObject_GetAttrString(o, s) _refcountDebug_PyObject_GetAttrString(o, s, #o, __LINE__)
#define PyTuple_New(size) _refcountDebug_PyTuple_New(size, "new tuple", __LINE__)
#define PyList_New(size) _refcountDebug_PyList_New(size, "new list", __LINE__)
#define PyDict_New() _refcountDebug_PyDict_New("new dict", __LINE__)
#define PyInt_FromLong(v) _refcountDebug_PyInt_FromLong(v, #v, __LINE__)
#define PyFloat_FromDouble(v) _refcountDebug_PyFloat_FromDouble(v, #v, __LINE__)
#define PyInstance_New(pyclass, args, kw) _refcountDebug_PyInstance_New(pyclass, args, kw, #pyclass, __LINE__)
#define PyObject_Call(callable_object, args, kw) _refcountDebug_PyObject_Call(callable_object, args, kw, #callable_object, __LINE__)
#define PyTuple_SetItem(tuple, i, item) _refcountDebug_PyTuple_SetItem(tuple, i, item, #item, __LINE__)
#define PyList_SetItem(list, i, item) _refcountDebug_PyList_SetItem(list, i, item, #item, __LINE__)
#undef Py_INCREF
#undef Py_DECREF
#undef Py_XINCREF
#undef Py_XDECREF
#define Py_INCREF(o) _refcountDebug_Py_INCREF(o, #o, __LINE__)
#define Py_DECREF(o) _refcountDebug_Py_DECREF(o, #o, __LINE__)
#define Py_XINCREF(o) _refcountDebug_Py_INCREF(o, #o, __LINE__)
#define Py_XDECREF(o) _refcountDebug_Py_DECREF(o, #o, __LINE__)

#define refcountDebug_report(title, showall) _refcountDebug_report(title, showall, __FILE__, __LINE__)

#else /* REFCOUNTDEBUG */

#define refcountDebug_report(title, showall)
#define refcountDebug_init()

#endif /* REFCOUNTDEBUG */

/******************************************************************************/



static void clib_message(CLIB_MSG_IMPORTANCE importance, const char *format, ...)
{
  va_list arglist;
  int imp = (int) importance;

  if (imp >= ((int) message_importance_threshold))
  {
    va_start(arglist, format);
    vfprintf(stderr, format, arglist);
  }
}


static PyObject *clib_setverbose(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "i", &message_importance_threshold))
  {
    return (NULL);
  }
  clib_message(CLIB_MSG_TRACE, "transsys.clib: verbosity level set to %d\n", message_importance_threshold);
  Py_INCREF(Py_None);
  return (Py_None);
}


static PyObject *clib_dummy(PyObject *self, PyObject *args)
{
  fprintf(stderr, "clib_dummy called\n");
  PyErr_SetString(PyExc_SystemError, "demo system error");
  return (NULL);
  /*
  Py_INCREF(Py_None);
  return (Py_None);
  */
}


static PyMethodDef clib_methods[] = {
  {"dummy", clib_dummy, METH_VARARGS, "dummy test function for clib development"},
  {"setverbose", clib_setverbose, METH_VARARGS, "set verbosity level for transsys.clib module"},
  {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC initclib(void)
{
  PyObject *clib_module;
  clib_module = Py_InitModule("paftol.clib", clib_methods);
  /* FIXME: should not ignore return value */
  PyModule_AddStringConstant(clib_module, "clib_api_version", clib_api_version);
}

/* don't forget to change the clib_api_version */
