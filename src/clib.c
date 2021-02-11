/* Copyright © 2020 The Board of Trustees of the Royal Botanic Gardens, Kew */

#include <Python.h>

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>


#define MAX_LINE_LENGTH 1000


typedef enum
{
  CLIB_MSG_TRACE,
  CLIB_MSG_DEVMESSAGE,
  CLIB_MSG_WARNING,
  CLIB_MSG_ERROR,
  CLIB_MSG_FATAL
} CLIB_MSG_IMPORTANCE;


typedef struct
{
  char *symbol;
  double **score;
} SYMBOL_SCORE_MATRIX;


typedef struct
{
  char *id;
  char *description;
  char *seq;
} BIOSEQUENCE;


typedef struct
{
  BIOSEQUENCE *seq0;
  BIOSEQUENCE *seq1;
  double score;
} PAIRWISE_ALIGNMENT;


typedef struct
{
  int i;
  int j;
  double **m;
  double score;
} BACKTRACK_POSITION;



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





static int find_symbol_index(const SYMBOL_SCORE_MATRIX *symbol_score_matrix, char symbol)
{
  int i;

  /* fprintf(stderr, "symbol list: %s\n", symbol_score_matrix->symbol); */
  for (i = 0; symbol_score_matrix->symbol[i] != '\0'; i++)
  {
    if (symbol == symbol_score_matrix->symbol[i])
    {
      return (i);
    }
  }
  return (-1);
}


static double find_symbol_score(const SYMBOL_SCORE_MATRIX *symbol_score_matrix, char sym0, char sym1)
{
  int i0, i1;

  i0 = find_symbol_index(symbol_score_matrix, (char) toupper(sym0));
  if (i0 == -1)
  {
    fprintf(stderr, "symbol 0 '%c' not in symbol score matrix\n", sym0);
    return (strtod("NaN", NULL));
  }
  i1 = find_symbol_index(symbol_score_matrix, (char) toupper(sym1));
  if (i1 == -1)
  {
    fprintf(stderr, "symbol 1 '%c' not in symbol score matrix\n", sym1);
    return (strtod("NaN", NULL));
  }
  return (symbol_score_matrix->score[i0][i1]);
}


static void write_biosequence_fasta(const BIOSEQUENCE *biosequence, FILE *f)
{
  fprintf(f, ">%s %s\n", biosequence->id, biosequence->description);
  fprintf(f, "%s\n", biosequence->seq);
}


static void write_pairwise_alignment_fasta(const PAIRWISE_ALIGNMENT *pairwise_alignment, FILE *f)
{
  write_biosequence_fasta(pairwise_alignment->seq0, f);
  write_biosequence_fasta(pairwise_alignment->seq1, f);
}


static int biosequence_length(const BIOSEQUENCE *biosequence)
{
  return (strlen(biosequence->seq));
}


static void free_matrix(double **m)
{
  free(m[0]);
  free(m);
}


static double **malloc_matrix(int num_rows, int num_columns)
{
  double **m;
  int r;

  m = (double **) malloc(num_rows * sizeof(double *));
  if (m == NULL)
  {
    return (NULL);
  }
  m[0] = (double *) malloc(num_rows * num_columns * sizeof(double));
  if (m[0] == NULL)
  {
    free(m);
    return (NULL);
  }
  for (r = 1; r < num_rows; r++)
  {
    m[r] = m[0] + r * num_columns;
  }
  return (m);
}


static void free_symbol_score_matrix(SYMBOL_SCORE_MATRIX *symbol_score_matrix)
{
  free_matrix(symbol_score_matrix->score);
  free(symbol_score_matrix->symbol);
  free(symbol_score_matrix);
}


static SYMBOL_SCORE_MATRIX *malloc_symbol_score_matrix(int num_symbols)
{
  SYMBOL_SCORE_MATRIX *symbol_score_matrix = (SYMBOL_SCORE_MATRIX *) malloc(sizeof(SYMBOL_SCORE_MATRIX));
  if (symbol_score_matrix == NULL)
  {
    return (NULL);
  }
  symbol_score_matrix->symbol = (char *) malloc(sizeof(char) * (num_symbols + 1));
  if (symbol_score_matrix->symbol == NULL)
  {
    free(symbol_score_matrix);
    return (NULL);
  }
  symbol_score_matrix->score = malloc_matrix(num_symbols, num_symbols);
  if (symbol_score_matrix->score == NULL)
  {
    free(symbol_score_matrix->symbol);
    free(symbol_score_matrix);
    return (NULL);
  }
  return (symbol_score_matrix);
}


static int get_num_symbols(const SYMBOL_SCORE_MATRIX *symbol_score_matrix)
{
  return (strlen(symbol_score_matrix->symbol));
}


static char *fget_next_line(char *s, int size, FILE *f)
{
  char *fgets_retval;

  fgets_retval = fgets(s, size, f);
  if (fgets_retval == NULL)
  {
    return (NULL);
  }
  while (s[0] == '#')
  {
    fgets_retval = fgets(s, size, f);
    if (fgets_retval == NULL)
    {
      return (NULL);
    }
  }
  return (fgets_retval);
}


static void free_biosequence(BIOSEQUENCE *biosequence)
{
  free(biosequence->seq);
  free(biosequence->description);
  free(biosequence->id);
  free(biosequence);
}


static BIOSEQUENCE *malloc_biosequence(int id_length, int description_length, int sequence_length)
{
  BIOSEQUENCE *biosequence;

  biosequence = (BIOSEQUENCE *) malloc(sizeof(BIOSEQUENCE));
  if (biosequence == NULL)
  {
    return (NULL);
  }
  biosequence->id = (char *) malloc((id_length + 1) * sizeof(char));
  if (biosequence->id == NULL)
  {
    free(biosequence);
    return (NULL);
  }
  biosequence->description = (char *) malloc((description_length + 1) * sizeof(char));
  if (biosequence->description == NULL)
  {
    free(biosequence->id);
    free(biosequence);
    return (NULL);
  }
  biosequence->seq = (char *) malloc((sequence_length + 1) * sizeof(char));
  if (biosequence->seq == NULL)
  {
    free(biosequence->description);
    free(biosequence->id);
    free(biosequence);
    return (NULL);
  }
  return (biosequence);
}


static BIOSEQUENCE *wrap_biosequence(char *id, char *description, char *seq)
{
  BIOSEQUENCE *biosequence;

  biosequence = (BIOSEQUENCE *) malloc(sizeof(BIOSEQUENCE));
  if (biosequence == NULL)
  {
    return (NULL);
  }
  biosequence->id = id;
  biosequence->description = description;
  biosequence->seq = seq;
  return (biosequence);
}


static BIOSEQUENCE *new_biosequence(const char *id, const char *description, const char *seq)
{
  BIOSEQUENCE *biosequence;

  biosequence = malloc_biosequence(strlen(id), strlen(description), strlen(seq));
  if (biosequence == NULL)
  {
    return (NULL);
  }
  strncpy(biosequence->id, id, strlen(id) + 1);
  strncpy(biosequence->description, description, strlen(description) + 1);
  strncpy(biosequence->seq, seq, strlen(seq) + 1);
  return (biosequence);
}


static BIOSEQUENCE *clone_biosequence(const BIOSEQUENCE *biosequence)
{
  return(new_biosequence(biosequence->id, biosequence->description, biosequence->seq));
}


void free_pairwise_alignment(PAIRWISE_ALIGNMENT *pairwise_alignment)
{
  free_biosequence(pairwise_alignment->seq0);
  free_biosequence(pairwise_alignment->seq1);
}


static PAIRWISE_ALIGNMENT *wrap_pairwise_alignment(BIOSEQUENCE *biosequence0, BIOSEQUENCE *biosequence1)
{
  PAIRWISE_ALIGNMENT *pairwise_alignment;

  pairwise_alignment = (PAIRWISE_ALIGNMENT *) malloc(sizeof(PAIRWISE_ALIGNMENT));
  if (pairwise_alignment == NULL)
  {
    return (NULL);
  }
  pairwise_alignment->seq0 = biosequence0;
  pairwise_alignment->seq1 = biosequence1;
  return (pairwise_alignment);
}


static PAIRWISE_ALIGNMENT *new_pairwise_alignment(const BIOSEQUENCE *biosequence0, const BIOSEQUENCE *biosequence1)
{
  PAIRWISE_ALIGNMENT *pairwise_alignment;
  BIOSEQUENCE *clone0, *clone1;

  clone0 = clone_biosequence(biosequence0);
  if (clone0 == NULL)
  {
    return (NULL);
  }
  clone1 = clone_biosequence(biosequence1);
  if (clone1 == NULL)
  {
    free_biosequence(clone0);
    return (NULL);
  }
  pairwise_alignment = wrap_pairwise_alignment(clone0, clone1);
  if (pairwise_alignment == NULL)
  {
    free_biosequence(clone1);
    free_biosequence(clone0);
    return (NULL);
  }
  return (pairwise_alignment);
}


static const char *next_nonspace(const char *s)
{
  while ((*s != '\0') && isspace(*s))
  {
    s++;
  }
  return (s);
}



static const char *next_double(const char *s, double *p)
{
  double d;
  char *e;

  s = next_nonspace(s);
  if (*s == '\0')
  {
    return (s);
  }
  d = strtod(s, &e);
  if (s == e)
  {
    fprintf(stderr, "invalid double: \"%s\"", s);
    return (s);
  }
  *p = d;
  return (e);
}


static int read_matrix_row(char buf[], SYMBOL_SCORE_MATRIX *symbol_score_matrix, int r)
{
  const char *s, *e;
  char row_symbol;
  int num_symbols = get_num_symbols(symbol_score_matrix);
  int j;
  double d;

  s = next_nonspace(buf);
  row_symbol = *s;
  if (row_symbol != symbol_score_matrix->symbol[r])
  {
    fprintf(stderr, "inconsistent row: expected symbol %c but got %c\n", symbol_score_matrix->symbol[r], row_symbol);
  }
  s++;
  for (j = 0; j < num_symbols; j++)
  {
    e = next_double(s, &d);
    if (s != e)
    {
      symbol_score_matrix->score[r][j] = d;
    }
    s = e;
  }
  return (0);
}


static int read_symbol_list(char buf[], char *symbol)
{
  int num_symbols = 0;
  const char *s;

  fprintf(stderr, "read_symbol_list: %s", buf);
  s = next_nonspace(buf);
  while (*s != '\0')
  {
    symbol[num_symbols++] = *s;
    s++;
    s = next_nonspace(s);
  }
  symbol[num_symbols] = '\0';
  return (num_symbols);
}


static SYMBOL_SCORE_MATRIX *read_symbol_score_matrix(FILE *f)
{
  SYMBOL_SCORE_MATRIX *symbol_score_matrix = NULL;
  int num_symbols, r;
  char buf[MAX_LINE_LENGTH], symbol_list[MAX_LINE_LENGTH];

  if (fget_next_line(buf, MAX_LINE_LENGTH, f) == NULL)
  {
    return (NULL);
  }
  num_symbols = read_symbol_list(buf, symbol_list);
  fprintf(stderr, "read_symbol_score_matrix: got %d symbols\n", num_symbols);
  symbol_score_matrix = malloc_symbol_score_matrix(num_symbols);
  for (r = 0; r < num_symbols; r++)
  {
    symbol_score_matrix->symbol[r] = symbol_list[r];
  }
  symbol_score_matrix->symbol[num_symbols] = '\0';
  for (r = 0; r < num_symbols; r++)
  {
    if (fget_next_line(buf, MAX_LINE_LENGTH, f) == NULL)
    {
      free_symbol_score_matrix(symbol_score_matrix);
      return (NULL);
    }
    read_matrix_row(buf, symbol_score_matrix, r);
  }
  return (symbol_score_matrix);
}


static SYMBOL_SCORE_MATRIX *make_ednafull_matrix()
{
  char ednafull_symbol[] = "ATGCSWRYKMBVHDNU";
  double ednafull_score[] = {
    5, -4, -4, -4, -4,  1,  1, -4, -4,  1, -4, -1, -1, -1, -2, -4,
    -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2,  5,
    -4, -4,  5, -4,  1, -4,  1, -4,  1, -4, -1, -1, -4, -1, -2, -4,
    -4, -4, -4,  5,  1, -4, -4,  1, -4,  1, -1, -1, -1, -4, -2, -4,
    -4, -4,  1,  1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1, -4,
    1,  1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1,  1,
    1, -4,  1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1, -4,
    -4,  1, -4,  1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1,  1,
    -4,  1,  1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1,  1,
    1, -4, -4,  1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1, -4,
    -4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1, -1,
    -1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1, -4,
    -1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1, -1,
    -1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1, -1,
    -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,
    -4,  5, -4, -4, -4,  1, -4,  1,  1, -4, -1, -4, -1, -1, -2,  5
  };
  size_t num_symbols = strlen(ednafull_symbol);
  size_t i, j;
  SYMBOL_SCORE_MATRIX *ednafull_matrix = malloc_symbol_score_matrix(num_symbols);

  if (ednafull_matrix == NULL)
  {
    return (NULL);
  }
  strcpy(ednafull_matrix->symbol, ednafull_symbol);
  for (i = 0; i < num_symbols; i++)
  {
    for (j = 0; j < num_symbols; j++)
    {
      ednafull_matrix->score[i][j] = ednafull_score[i * num_symbols + j];
    }
  }
  return (ednafull_matrix);
}


static void write_symbol_score_matrix(FILE *f, const SYMBOL_SCORE_MATRIX *symbol_score_matrix)
{
  int num_symbols = get_num_symbols(symbol_score_matrix);
  int i, j;

  for (i = 0; i < num_symbols; i++)
  {
    fprintf(f, "  %c", symbol_score_matrix->symbol[i]);
  }
  fprintf(f, "\n");
  for (i = 0; i < num_symbols; i++)
  {
    fprintf(f, "%c ", symbol_score_matrix->symbol[i]);
    for (j = 0; j < num_symbols; j++)
    {
      fprintf(f, "  %f", symbol_score_matrix->score[i][j]);
    }
    fprintf(f, "\n");
  }
}


static char *right_padded_extension(const BIOSEQUENCE *biosequence, int extended_length)
{
  char gapchar = '-';
  char *extended_seq;
  int i;

  extended_seq = (char *) malloc(sizeof(char) * (extended_length + 1));
  if (extended_seq == NULL)
  {
    return (NULL);
  }
  strcpy(extended_seq, biosequence->seq);
  for (i = strlen(biosequence->seq); i < extended_length; i++)
  {
    extended_seq[i] = gapchar;
  }
  extended_seq[extended_length] = '\0';
  return (extended_seq);
}


static PAIRWISE_ALIGNMENT *align_by_padding_right(const BIOSEQUENCE *seq0, const BIOSEQUENCE *seq1, const SYMBOL_SCORE_MATRIX *symbol_score_matrix, double gap_creation_penalty, double gap_extension_penalty)
{
  int length0, length1;
  char *extended_seq;
  BIOSEQUENCE *aligned0, *aligned1;;
  PAIRWISE_ALIGNMENT *pairwise_alignment;

  length0 = biosequence_length(seq0);
  length1 = biosequence_length(seq1);
  aligned0 = clone_biosequence(seq0);
  aligned1 = clone_biosequence(seq1);
  if (length0 < length1)
  {
    extended_seq = right_padded_extension(seq0, biosequence_length(seq1));
    if (extended_seq == NULL)
    {
      free_biosequence(aligned0);
      free_biosequence(aligned1);
      return (NULL);
    }
    free(aligned0->seq);
    aligned0->seq = extended_seq;
  }
  else if (length0 > length1)
  {
    extended_seq = right_padded_extension(seq1, biosequence_length(seq0));
    if (extended_seq == NULL)
    {
      free_biosequence(aligned0);
      free_biosequence(aligned1);
      return (NULL);
    }
    free(aligned1->seq);
    aligned1->seq = extended_seq;
  }
  pairwise_alignment = wrap_pairwise_alignment(aligned0, aligned1);
  if (pairwise_alignment == NULL)
  {
    free_biosequence(aligned0);
    free_biosequence(aligned1);
    return (NULL);
  }
  pairwise_alignment->score = 0.0;
  return (pairwise_alignment);
}


static double **print_matrix(FILE *f, double **m, int nrow, int ncol)
{
  int i, j;

  for (i = 0; i < nrow; i++)
  {
    for (j = 0; j < ncol; j++)
    {
      fprintf(f, "  %5.1f", m[i][j]);
    }
    fprintf(f, "\n");
  }
  fprintf(f, "\n");
  return (m);
}


static double **fill_matrix(double **m, double v, int nrow, int ncol)
{
  int i, j;

  for (i = 0; i < nrow; i++)
  {
    for (j = 0; j < ncol; j++)
    {
      m[i][j] = v;
    }
  }
  return (m);
}


static double max3(double x, double y, double z)
{
  double m;

  m = x > y ? x : y;
  return (m > z ? m : z);
}


static char *reverse_string(char *s)
{
  size_t j = strlen(s);
  size_t i = 0;
  char c;

  if (j == 0)
  {
    return (s);
  }
  j--;
  while (i < j)
  {
    /* fprintf(stderr, "reversing, i = %lu, j = %lu\n", (unsigned long) i, (unsigned long) j); */
    c = s[i];
    s[i++] = s[j];
    s[j--] = c;
  }
  return (s);
}


static PAIRWISE_ALIGNMENT *align_semiglobal(const BIOSEQUENCE *seq0, const BIOSEQUENCE *seq1, const SYMBOL_SCORE_MATRIX *symbol_score_matrix, double gap_creation_penalty, double gap_extension_penalty)
{
  char gapchar = '-';
  BIOSEQUENCE *aligned_seq0, *aligned_seq1;
  PAIRWISE_ALIGNMENT *pairwise_alignment = NULL;
  int l0 = strlen(seq0->seq);
  int l1 = strlen(seq1->seq);
  int i, j, k;
  char symbol0, symbol1;
  char *aln0, *aln1;
  double **m, **m0, **m1;
  double symbol_score, alignment_score;
  double nan;
  BACKTRACK_POSITION backtrack_position;

  nan = strtod("NaN", NULL);
  m = malloc_matrix(l0 + 1, l1 + 1);
  if (m == NULL)
  {
    return (NULL);
  }
  m0 = malloc_matrix(l0 + 1, l1 + 1);   if (m0 == NULL)
  {
    free_matrix(m);
    return (NULL);
  }
  m1 = malloc_matrix(l0 + 1, l1 + 1);
  if (m1 == NULL)
  {
    free_matrix(m0);
    free_matrix(m);
    return (NULL);
  }
  aln0 = (char *) malloc(l0 + l1 + 1);
  if (aln0 == NULL)
  {
    free_matrix(m1);
    free_matrix(m0);
    free_matrix(m);
    return (NULL);
  }
  aln1 = (char *) malloc(l0 + l1 + 1);
  if (aln0 == NULL)
  {
    free(aln0);
    free_matrix(m1);
    free_matrix(m0);
    free_matrix(m);
    return (NULL);
  }
  fill_matrix(m, nan, l0 + 1, l1 + 1);
  fill_matrix(m0, nan, l0 + 1, l1 + 1);
  fill_matrix(m1, nan, l0 + 1, l1 + 1);
  for (i = 0; i <= l0; i++)
  {
    m[i][0] = 0.0;
    m0[i][0] = -DBL_MAX;
    m1[i][0] = -DBL_MAX;
  }
  for (j = 0; j <= l1; j++)
  {
    m[0][j] = 0.0;
    m0[0][j] = -DBL_MAX;
    m1[0][j] = -DBL_MAX;
  }
  for (i = 1; i <= l0; i++)
  {
    for (j = 1; j <= l1; j++)
    {
      symbol0 = seq0->seq[i - 1];
      symbol1 = seq1->seq[j - 1];
      symbol_score = find_symbol_score(symbol_score_matrix, symbol0, symbol1);
      m[i][j] = max3(m[i - 1][j - 1], m0[i - 1][j - 1], m1[i - 1][j - 1]) + symbol_score;
      m0[i][j] = max3(m[i][j - 1] - gap_creation_penalty, m0[i][j - 1] - gap_extension_penalty, m1[i][j - 1] - gap_creation_penalty);
      m1[i][j] = max3(m[i - 1][j] - gap_creation_penalty, m0[i - 1][j] - gap_creation_penalty, m1[i - 1][j] - gap_extension_penalty);
    }
  }
  /*
  print_matrix(stderr, m, l0 + 1, l1 + 1);
  print_matrix(stderr, m0, l0 + 1, l1 + 1);
  print_matrix(stderr, m1, l0 + 1, l1 + 1);
  */
  backtrack_position.i = 0;
  backtrack_position.j = l1;
  backtrack_position.m = m;
  for (i = 0; i <= l0; i++)
  {
    if (backtrack_position.m[backtrack_position.i][backtrack_position.j] < m[i][l1])
    {
      backtrack_position.i = i;
      backtrack_position.j = l1;
      backtrack_position.m = m;
    }
    if (backtrack_position.m[backtrack_position.i][backtrack_position.j] < m0[i][l1])
    {
      backtrack_position.i = i;
      backtrack_position.j = l1;
      backtrack_position.m = m0;
    }
    if (backtrack_position.m[backtrack_position.i][backtrack_position.j] < m1[i][l1])
    {
      backtrack_position.i = i;
      backtrack_position.j = l1;
      backtrack_position.m = m1;
    }
  }
  for (j = 0; j <= l1; j++)
  {
    if (backtrack_position.m[backtrack_position.i][backtrack_position.j] < m[l0][j])
    {
      backtrack_position.i = l0;
      backtrack_position.j = j;
      backtrack_position.m = m;
    }
    if (backtrack_position.m[backtrack_position.i][backtrack_position.j] < m0[l0][j])
    {
      backtrack_position.i = l0;
      backtrack_position.j = j;
      backtrack_position.m = m0;
    }
    if (backtrack_position.m[backtrack_position.i][backtrack_position.j] < m1[l0][j])
    {
      backtrack_position.i = l0;
      backtrack_position.j = j;
      backtrack_position.m = m1;
    }
  }
  alignment_score = backtrack_position.m[backtrack_position.i][backtrack_position.j];
  /*
  fprintf(stderr, "l0 = %d, l1 = %d, matrices: m = %p, m0 = %p, m1 = %p\n", l0, l1, m, m0, m1);
  fprintf(stderr, "backtrack_position: i = %d, j = %d, m = %p, score = %f\n", backtrack_position.i, backtrack_position.j, backtrack_position.m, backtrack_position.m[backtrack_position.i][backtrack_position.j]);
  */
  k = 0;
  for (i = l0 - 1; i >= backtrack_position.i; i--)
  {
    aln0[k] = seq0->seq[i];
    aln1[k] = gapchar;
    k++;
  }
  for (j = l1 - 1; j >= backtrack_position.j; j--)
  {
    aln0[k] = gapchar;
    aln1[k] = seq1->seq[j];
    k++;
  }
  while ((backtrack_position.i > 0) && (backtrack_position.j > 0))
  {
    symbol0 = seq0->seq[backtrack_position.i - 1];
    symbol1 = seq1->seq[backtrack_position.j - 1];
    symbol_score = find_symbol_score(symbol_score_matrix, symbol0, symbol1);
    if (backtrack_position.m == m)
    {
      /* m[i][j] = max3(m[i - 1][j - 1], m0[i - 1][j - 1], m1[i - 1][j - 1]) + symbol_score; */
      if (m[backtrack_position.i][backtrack_position.j] == m[backtrack_position.i - 1][backtrack_position.j - 1] + symbol_score)
      {
        /* fprintf(stderr, "(%d, %d): m -> m\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.i--;
        backtrack_position.j--;
        backtrack_position.m = m;
      }
      else if (m[backtrack_position.i][backtrack_position.j] == m0[backtrack_position.i - 1][backtrack_position.j - 1] + symbol_score)
      {
        /* fprintf(stderr, "(%d, %d): m -> m0\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.i--;
        backtrack_position.j--;
        backtrack_position.m = m0;
      }
      else if (m[backtrack_position.i][backtrack_position.j] == m1[backtrack_position.i - 1][backtrack_position.j - 1] + symbol_score)
      {
        /* fprintf(stderr, "(%d, %d): m -> m1\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.i--;
        backtrack_position.j--;
        backtrack_position.m = m1;
      }
      else
      {
        fprintf(stderr, "no backtracking step from m[%d][%d]\n", backtrack_position.i, backtrack_position.j);
      }
      aln0[k] = seq0->seq[backtrack_position.i];
      aln1[k] = seq1->seq[backtrack_position.j];
    }
    else if (backtrack_position.m == m0)
    {
      /* m0[i][j] = max3(m[i][j - 1] - gap_creation_penalty, m0[i][j - 1] - gap_extension_penalty, m1[i][j - 1] - gap_creation_penalty); */
      if (m0[backtrack_position.i][backtrack_position.j] == m[backtrack_position.i][backtrack_position.j - 1] - gap_creation_penalty)
      {
        /* fprintf(stderr, "(%d, %d): m0 -> m\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.j--;
        backtrack_position.m = m;
      }
      else if (m0[backtrack_position.i][backtrack_position.j] == m0[backtrack_position.i][backtrack_position.j - 1] - gap_extension_penalty)
      {
        /* fprintf(stderr, "(%d, %d): m0 -> m0\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.j--;
        backtrack_position.m = m0;
      }
      else if (m0[backtrack_position.i][backtrack_position.j] == m1[backtrack_position.i][backtrack_position.j - 1] - gap_creation_penalty)
      {
        /* fprintf(stderr, "(%d, %d): m0 -> m1\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.j--;
        backtrack_position.m = m1;
      }
      else
      {
        fprintf(stderr, "no backtracking step from m0[%d][%d]\n", backtrack_position.i, backtrack_position.j);
      }
      aln0[k] = gapchar;
      aln1[k] = seq1->seq[backtrack_position.j];
    }
    else if (backtrack_position.m == m1)
    {
      /* m1[i][j] = max3(m[i - 1][j] - gap_creation_penalty, m0[i - 1][j] - gap_creation_penalty, m1[i - 1][j] - gap_extension_penalty); */
      if (m1[backtrack_position.i][backtrack_position.j] == m[backtrack_position.i - 1][backtrack_position.j] - gap_creation_penalty)
      {
        /* fprintf(stderr, "(%d, %d): m1 -> m\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.i--;
        backtrack_position.m = m;
      }
      else if (m1[backtrack_position.i][backtrack_position.j] == m0[backtrack_position.i - 1][backtrack_position.j] - gap_creation_penalty)
      {
        /* fprintf(stderr, "(%d, %d): m1 -> m0\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.i--;
        backtrack_position.m = m0;
      }
      else if (m1[backtrack_position.i][backtrack_position.j] == m1[backtrack_position.i - 1][backtrack_position.j] - gap_extension_penalty)
      {
        /* fprintf(stderr, "(%d, %d): m1 -> m1\n", backtrack_position.i, backtrack_position.j); */
        backtrack_position.i--;
        backtrack_position.m = m1;
      }
      else
      {
        fprintf(stderr, "no backtracking step from m0[%d][%d]\n", backtrack_position.i, backtrack_position.j);
      }
      aln0[k] = seq0->seq[backtrack_position.i];
      aln1[k] = gapchar;
    }
    else
    {
      fprintf(stderr, "internal error: backtracking from unknown matrix %p (m = %p, m0 = %p, m1 = %p)\n", (void *) backtrack_position.m, (void *) m, (void *) m0, (void *) m1);
    }
    /* fprintf(stderr, "aln0[%d] = %c, aln1[%d] = %c\n", k, aln0[k], k, aln1[k]); */
    k++;
    /* fprintf(stderr, "backtrack_position: i = %d, j = %d, m = %p, score = %f\n", backtrack_position.i, backtrack_position.j, backtrack_position.m, backtrack_position.m[backtrack_position.i][backtrack_position.j]); */
  }
  while (backtrack_position.i > 0)
  {
    aln0[k] = seq0->seq[--backtrack_position.i];
    aln1[k] = gapchar;
    k++;
  }
  while (backtrack_position.j > 0)
  {
    aln0[k] = gapchar;
    aln1[k] = seq1->seq[--backtrack_position.j];
    k++;
  }
  aln0[k] = '\0';
  aln1[k] = '\0';
  /*
  fprintf(stderr, "aln0 rev: %s\n", aln0);
  fprintf(stderr, "aln1 rev: %s\n", aln1);
  */
  reverse_string(aln0);
  reverse_string(aln1);
  /*
  fprintf(stderr, "aln0: %s\n", aln0);
  fprintf(stderr, "aln1: %s\n", aln1);
  */
  aligned_seq0 = new_biosequence(seq0->id, seq0->description, aln0);
  aligned_seq1 = new_biosequence(seq1->id, seq1->description, aln1);
  free(aln1);
  free(aln0);
  free_matrix(m);
  free_matrix(m0);
  free_matrix(m1);
  pairwise_alignment = wrap_pairwise_alignment(aligned_seq0, aligned_seq1);
  pairwise_alignment->score = alignment_score;
  if (pairwise_alignment == NULL)
  {
    free_biosequence(aligned_seq0);
    free_biosequence(aligned_seq1);
  }
  return (pairwise_alignment);
}

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


PyObject *clib_setverbose(PyObject *self, PyObject *args)
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
  PyObject *py_helloworld = PyString_FromString("hello world");
  return(py_helloworld);
  /*
  PyErr_SetString(PyExc_SystemError, "demo system error");
  return (NULL);
  */
  /*
  Py_INCREF(Py_None);
  return (Py_None);
  */
}


static PyObject *clib_align_semiglobal(PyObject *self, PyObject *args)
{
  const char *s0, *s1;
  double gap_creation_penalty, gap_extension_penalty;
  BIOSEQUENCE *biosequence0, *biosequence1;
  SYMBOL_SCORE_MATRIX *symbol_score_matrix;
  PAIRWISE_ALIGNMENT *pairwise_alignment;
  PyObject *r;

  /* fprintf(stderr, "starting\n"); */
  if (!PyArg_ParseTuple(args, "ss", &s0, &s1))
  {
    return (NULL);
  }
  /* fprintf(stderr, "got args\n"); */
  /* fprintf(stderr, "s0 = %s, s1 = %s\n", s0, s1); */
  biosequence0 = new_biosequence("seq0", "seq0", s0);
  if (biosequence0 == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "failed to allocate biosequence0");
    return (NULL);
  }
  biosequence1 = new_biosequence("seq1", "seq1", s1);
  if (biosequence1 == NULL)
  {
    free_biosequence(biosequence1);
    PyErr_SetString(PyExc_MemoryError, "failed to allocate biosequence1");
    return (NULL);
  }
  /* fprintf(stderr, "made biosequences\n"); */
  symbol_score_matrix = make_ednafull_matrix();
  if (symbol_score_matrix == NULL)
  {
    free_biosequence(biosequence0);
    free_biosequence(biosequence1);
    PyErr_SetString(PyExc_MemoryError, "failed to allocate symbol score matrix");
    return (NULL);
  }
  /* fprintf(stderr, "made symbol scoring matrix\n"); */
  gap_creation_penalty = 10.0;
  gap_extension_penalty = 0.5;
  pairwise_alignment = align_semiglobal(biosequence0, biosequence1, symbol_score_matrix, gap_creation_penalty, gap_extension_penalty);
  /* fprintf(stderr, "aligned biosequences\n"); */
  free_biosequence(biosequence0);
  free_biosequence(biosequence1);
  free_symbol_score_matrix(symbol_score_matrix);
  r = Py_BuildValue("ssd", pairwise_alignment->seq0->seq, pairwise_alignment->seq1->seq, pairwise_alignment->score);
  /* fprintf(stderr, "made return value tuple\n"); */
  free_pairwise_alignment(pairwise_alignment);
  /* fprintf(stderr, "freed alignment\n"); */
  return (r);
}


static BIOSEQUENCE *extract_biosequence_from_list(PyObject *python_list, Py_ssize_t i)
{
  PyObject *python_seqstr;
  BIOSEQUENCE *biosequence;
  char *s;

  python_seqstr = PySequence_GetItem(python_list, i);
  if (python_seqstr == NULL)
  {
    PyErr_SetString(PyExc_RuntimeError, "failed to retrieve element from list");
    return (NULL);
  }
  if (!PyString_Check(python_seqstr))
  {
    PyErr_SetString(PyExc_RuntimeError, "retrieved non-string element from list");
    Py_DECREF(python_seqstr);
    return (NULL);
  }
  s = PyString_AsString(python_seqstr);
  if (s == NULL)
  {
    PyErr_SetString(PyExc_RuntimeError, "PyString_AsString failed");
    Py_DECREF(python_seqstr);
    return (NULL);
  }
  /* fprintf(stderr, "extracted: #%ld: %s\n", (long) i, s); */
  biosequence = new_biosequence("seq0", "seq0", s);
  if (biosequence == NULL)
  {
    PyErr_SetString(PyExc_MemoryError, "failed to allocate biosequence");
    Py_DECREF(python_seqstr);
    return (NULL);
  }
  Py_DECREF(python_seqstr);
  return (biosequence);
}


static PyObject *clib_semiglobal_alignment_series(PyObject *self, PyObject *args)
{
  double gap_creation_penalty, gap_extension_penalty;
  BIOSEQUENCE *biosequence0, *biosequence1;
  SYMBOL_SCORE_MATRIX *symbol_score_matrix;
  PAIRWISE_ALIGNMENT *pairwise_alignment;
  PyObject *result_list, *alignment_tuple, *python_sequence_list, *python_symbol_score_matrix;
  Py_ssize_t num_sequences, i;

  /* fprintf(stderr, "starting\n"); */
  if (!PyArg_ParseTuple(args, "OddO", &python_sequence_list, &gap_creation_penalty, &gap_extension_penalty, &python_symbol_score_matrix))
  {
    return (NULL);
  }
  /* fprintf(stderr, "gap_creation_penalty = %f, gap_extension_penalty = %f\n", gap_creation_penalty, gap_extension_penalty); */
  if (!PySequence_Check(python_sequence_list))
  {
    PyErr_SetString(PyExc_TypeError, "sequence_list (arg 0) was not a sequence");
    return (NULL);
  }
  if (!PyMapping_Check(python_symbol_score_matrix))
  {
    PyErr_SetString(PyExc_TypeError, "symbol_score_matrix (arg 3) was not a mapping");
    return (NULL);
  }
  num_sequences = PySequence_Length(python_sequence_list);
  if (num_sequences < 2)
  {
    PyErr_SetString(PyExc_TypeError, "sequence_list has too few elements");
    return (NULL);
  }
  fprintf(stderr, "ignoring symbol_score_matrix and using hard-coded EDNAFULL instead\n");
  symbol_score_matrix = make_ednafull_matrix();
  Py_INCREF(python_sequence_list);
  Py_INCREF(python_symbol_score_matrix);
  biosequence1 = extract_biosequence_from_list(python_sequence_list, 0);
  if (biosequence1 == NULL)
  {
    free_symbol_score_matrix(symbol_score_matrix);
    Py_DECREF(python_sequence_list);
    Py_DECREF(python_symbol_score_matrix);
    return (NULL);
  }
  result_list = PyList_New(0);
  if (result_list == NULL)
  {
    free_symbol_score_matrix(symbol_score_matrix);
    Py_DECREF(python_sequence_list);
    Py_DECREF(python_symbol_score_matrix);
    return (NULL);
  }
  for (i = 1; i < num_sequences; i++)
  {
    biosequence0 = biosequence1;
    biosequence1 = extract_biosequence_from_list(python_sequence_list, i);
    /* fprintf(stderr, "aligning: %s, %s\n", biosequence0->seq, biosequence1->seq); */
    if (biosequence1 == NULL)
    {
      free_biosequence(biosequence0);
      free_symbol_score_matrix(symbol_score_matrix);
      Py_DECREF(python_sequence_list);
      Py_DECREF(python_symbol_score_matrix);
      return (NULL);
    }
    pairwise_alignment = align_semiglobal(biosequence0, biosequence1, symbol_score_matrix, gap_creation_penalty, gap_extension_penalty);
    if (pairwise_alignment == NULL)
    {
      free_biosequence(biosequence0);
      free_biosequence(biosequence1);
      free_symbol_score_matrix(symbol_score_matrix);
      Py_DECREF(python_sequence_list);
      Py_DECREF(python_symbol_score_matrix);
      PyErr_SetString(PyExc_MemoryError, "failed to allocate pairwise alignment");
      return (NULL);
    }
    alignment_tuple = Py_BuildValue("ssd", pairwise_alignment->seq0->seq, pairwise_alignment->seq1->seq, pairwise_alignment->score);
    free_pairwise_alignment(pairwise_alignment);
    if (alignment_tuple == NULL)
    {
      free_biosequence(biosequence0);
      free_biosequence(biosequence1);
      free_symbol_score_matrix(symbol_score_matrix);
      Py_DECREF(python_sequence_list);
      Py_DECREF(python_symbol_score_matrix);
      return (NULL);
    }
    if (PyList_Append(result_list, alignment_tuple) != 0)
    {
      free_biosequence(biosequence0);
      free_biosequence(biosequence1);
      free_symbol_score_matrix(symbol_score_matrix);
      Py_DECREF(python_sequence_list);
      Py_DECREF(python_symbol_score_matrix);
      return (NULL);
    }
    free_biosequence(biosequence0);
  }
  free_biosequence(biosequence1);
  free_symbol_score_matrix(symbol_score_matrix);
  Py_DECREF(python_sequence_list);
  Py_DECREF(python_symbol_score_matrix);
  return (result_list);
}


static PyMethodDef clib_methods[] = {
  {"dummy", clib_dummy, METH_VARARGS, "dummy test function for clib development"},
  {"align_semiglobal", clib_align_semiglobal, METH_VARARGS, "compute semiglobal alignment of two sequences"},
  {"semiglobal_alignment_series", clib_semiglobal_alignment_series, METH_VARARGS, "compute consecutive series of semiglobal alignments"},
  {"setverbose", clib_setverbose, METH_VARARGS, "set verbosity level for paftol.clib module"},
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
