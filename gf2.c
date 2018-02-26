#include <Python.h>

#include "gf2.h"

unsigned int gf2bitlength(GF2 a) {
	if(a == 0) return 0;
	unsigned int n = 1;
	while(((a >> (n-1)) >> 1) > 0) {
		n <<= 1;
	}
	n >>= 1;
	unsigned int b = n;
	unsigned int nnew;
	while(b > 0) {
		b >>= 1;
		nnew = n ^ b;
		if((a >> nnew) > 0) n = nnew;
	}
	return n + 1;
}

GF2 gf2mulmod(GF2 a, GF2 b, GF2 m) {
	GF2 b2;
	GF2 c = 0;
	while(a > 0) {
		if((a & 1) != 0){
			c ^= b;
		}
		b <<= 1;
		b2 = b ^ m;
		if(b2 < b) {
			b = b2;
		}
		a >>= 1;
	}
	return c;
}

GF2 gf2powmod(GF2 a, GF2 k, GF2 m) {
	GF2 c = 1;
	while(k > 0) {
		if((k & 1) != 0) {
			c = gf2mulmod(c, a, m);
		}
		a = gf2mulmod(a, a, m);
		k >>= 1;
	}
	return c;
}


GF2 gf2mod(GF2 a, GF2 m) {
	int i = 0;
	while(m < a) {
		m <<= 1;
		++i;
	}

	GF2 b;
	while(i >= 0) {
		b = a ^ m;
		if(b < a)
			a = b;
		m >>= 1;
		--i;
	}
	return a;
}

GF2* gf2divmodvect(GF2 avec[], GF2 dvec[], size_t len) {
	GF2 q, test;
	unsigned int na, nd;
	int i;
	size_t j;
	na = gf2bitlength(avec[0]);
	nd = gf2bitlength(dvec[0]);
	i = na - nd;
	q = 0;
	test = (GF2)1 << (na-1);
	while(i >= 0) {
		if((avec[0] & test) != 0) {
			for(j = 0; j < len; ++j) {
				avec[j] = avec[j] ^ (dvec[j] << i);
			}
			q |= ((GF2)1 << i);
		}
		--i;
		test >>= 1;
	}
	
	return avec;
}

GF2* gf2exteuc(GF2 a, GF2 b, GF2* res) {
	GF2 arowData[3];
	GF2* arow = arowData;
	GF2* brow = res;
	GF2* rrow;
	arow[0] = a; arow[1] = 1; arow[2] = 0;
	brow[0] = b; brow[1] = 0; brow[2] = 1;
	while(1) {
		rrow = gf2divmodvect(arow, brow, 3);
		if(rrow[0] == 0) break;
		arow = brow;
		brow = rrow;
	}
	res[0] = brow[0]; res[1] = brow[1]; res[2] = brow[2];
	return res;
}

GF2 gf2modinv(GF2 x, GF2 m, bool* err) {
	GF2 extEucRes[3];
	gf2exteuc(x, m, extEucRes);
	if(extEucRes[0] != 1) {
		*err = true;
		return -1;
	}
	return extEucRes[1];
}

GF2 gf2gcd(GF2 x, GF2 m) {
	GF2 extEucRes[3];
	gf2exteuc(x, m, extEucRes);
	return extEucRes[0];
}

GF2 gf2divmod(GF2 a, GF2 d, GF2 m, bool* err) {
	return gf2mulmod(a, gf2modinv(d, m, err), m);
}

GF2 gf2add(GF2 a, GF2 b, unsigned int lgsize) {
	return (a ^ b) & (((GF2)1 << lgsize) - (GF2)1);
}

GF2 gf2sub(GF2 a, GF2 b, unsigned int lgsize) {
	return (a ^ b) & (((GF2)1 << lgsize) - (GF2)1);
}

bool gf2CheckIrreduciblePolynomial(GF2 poly, unsigned int lgsize) {
	GF2 i;
	for(i = 2; i < ((GF2)1 << lgsize) - 1; ++i) {
		if(gf2mod(poly, i) == 0) return false;
	}
	return true;
}

bool gf2CheckIrreduciblePolynomialRabin(GF2 poly, unsigned int lgsize) {
	GF2 h = gf2powmod((GF2)2, (GF2)1 << (lgsize/2), poly);
	if(gf2gcd(poly, h ^ 2) != 1) return false;
	
	h = gf2powmod((GF2)2, (GF2)1 << lgsize, poly);
	if(gf2mod(h, poly) != 2) return false;
	
	return true;
}

bool gf2CheckGeneratorPolynomial(GF2 gen, GF2 mod, unsigned int lgsize) {
	GF2* binPows = (GF2*) malloc(lgsize * sizeof(GF2));
	binPows[0] = gen;

	size_t i;
	for(i = 1; i < lgsize; ++i) {
		binPows[i] = gf2mulmod(binPows[i-1], binPows[i-1], mod);
	}

	GF2 j;
	bool bit;
	GF2 product;
	for(j = 2; j < ((GF2)1 << lgsize) - 1; ++j) {
		product = 1;
		for(i = 0; i < lgsize; ++i) {
			bit = (j >> i) & 1;
			if(bit) {
				product = gf2mulmod(product, binPows[i], mod);
			}
		}
		if(product == 1) {
			free(binPows);
			return false;
		}
	}
	free(binPows);
	return true;
}

static PyObject* _gf2bitlength(PyObject *self, PyObject *args) {
	GF2 a;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "K", &a))
		return NULL;

	return Py_BuildValue("i", gf2bitlength(a));
}

static PyObject* _gf2mulmod(PyObject *self, PyObject *args) {
	GF2 a, b, m;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKK", &a, &b, &m))
		return NULL;
	
	return Py_BuildValue("K", gf2mulmod(a, b, m));
}

static PyObject* _gf2powmod(PyObject *self, PyObject *args) {
	GF2 a, k, m;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKK", &a, &k, &m))
		return NULL;
	
	return Py_BuildValue("K", gf2powmod(a, k, m));
}

static PyObject* _gf2mod(PyObject *self, PyObject *args) {
	GF2 a, m;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KK", &a, &m))
		return NULL;
	
	return Py_BuildValue("K", gf2mod(a, m));
}

static PyObject* _gf2modinv( PyObject *self, PyObject *args ) {
	GF2 x, m;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KK", &x, &m))
		return NULL;
	
	bool err = false;
	GF2 c = gf2modinv(x, m, &err);
	if(err) {
		PyErr_SetString(PyExc_ValueError, "Not relatively prime");
		return NULL;
	}
	return Py_BuildValue("K", c);
}

static PyObject* _gf2gcd( PyObject *self, PyObject *args ) {
	GF2 x, m;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KK", &x, &m))
		return NULL;
	
	GF2 c = gf2gcd(x, m);
	return Py_BuildValue("K", c);
}

static PyObject* _gf2add( PyObject *self, PyObject *args ) {
	GF2 a, b;
	unsigned int lgsize;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKI", &a, &b, &lgsize))
		return NULL;
	
	return Py_BuildValue("K", gf2add(a, b, lgsize));
}

static PyObject* _gf2sub( PyObject *self, PyObject *args ) {
	GF2 a, b;
	unsigned int lgsize;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKI", &a, &b, &lgsize))
		return NULL;
	
	return Py_BuildValue("K", gf2sub(a, b, lgsize));
}

static PyObject* _gf2divmod( PyObject *self, PyObject *args ) {
	GF2 a, b, m;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKK", &a, &b, &m))
		return NULL;
	
	bool err = false;
	GF2 c = gf2divmod(a, b, m, &err);
	if(err) {
		PyErr_SetString(PyExc_ValueError, "Division Error");
		return NULL;
	}

	return Py_BuildValue("K", c);
}

static PyObject* _gf2CheckIrreduciblePolynomial( PyObject *self, PyObject *args ) {
	GF2 poly;
	unsigned int lgsize;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KI", &poly, &lgsize))
		return NULL;
	
	if(gf2CheckIrreduciblePolynomial(poly, lgsize)) Py_RETURN_TRUE;
	Py_RETURN_FALSE;
}

static PyObject* _gf2CheckIrreduciblePolynomialRabin( PyObject *self, PyObject *args ) {
	GF2 poly;
	unsigned int lgsize;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KI", &poly, &lgsize))
		return NULL;
	
	if(gf2CheckIrreduciblePolynomialRabin(poly, lgsize)) Py_RETURN_TRUE;
	Py_RETURN_FALSE;
}

static PyObject* _gf2CheckGeneratorPolynomial( PyObject *self, PyObject *args ) {
	GF2 gen, mod;
	unsigned int lgsize;

	// Parse the input tuple
	if (!PyArg_ParseTuple(args, "KKI", &gen, &mod, &lgsize))
		return NULL;
	
	if(gf2CheckGeneratorPolynomial(gen, mod, lgsize)) Py_RETURN_TRUE;
	Py_RETURN_FALSE;
}

static PyMethodDef gf2_funcs[] = {
	{"_gf2mulmod", _gf2mulmod, METH_VARARGS, "Computes ``a * x mod m``."},
	{"_gf2powmod", _gf2powmod, METH_VARARGS, "Computes `a` to the `k`th power, in the quotient ring GF(2)[x]/m(x)"},
	{"_gf2mod", _gf2mod, METH_VARARGS, "Computes ``a mod m``."},
	{"_gf2modinv", _gf2modinv, METH_VARARGS, "Computes the value ``b`` such that ``x * b mod m == 1``."},
	{"_gf2add", _gf2add, METH_VARARGS, "Computes ``a + b of GF(2^lgsize)`."},
	{"_gf2sub", _gf2sub, METH_VARARGS, "Computes ``a - b of GF(2^lgsize)``."},
	{"_gf2divmod", _gf2divmod, METH_VARARGS, "Computes ``a / b of GF(2^lgsize)``."},
	{"gf2CheckIrreduciblePolynomial", _gf2CheckIrreduciblePolynomial, METH_VARARGS, "Checks if poly is an irreducible polynomial of the field GF(2^lgsize)``."},
	{"gf2CheckIrreduciblePolynomialRabin", _gf2CheckIrreduciblePolynomialRabin, METH_VARARGS, "Checks if poly is an irreducible polynomial of the field GF(2^lgsize)``."},
	{"gf2CheckGeneratorPolynomial", _gf2CheckGeneratorPolynomial, METH_VARARGS, "Checks if gen is a generator of the field GF(2^lgsize) %% mod``."},
	{"_gf2bitlength", _gf2bitlength, METH_VARARGS, "Computes the length of a polynomial coefficient bit vector = degree + 1."},
	{"_gf2gcd", _gf2gcd, METH_VARARGS, "Computes the GCD of a and b"},
	{NULL, NULL, 0, NULL}
};

static struct PyModuleDef gf2_definition = { 
	PyModuleDef_HEAD_INIT,
	"gf2c",
	"",
	-1, 
	gf2_funcs
};

PyMODINIT_FUNC PyInit_gf2c(void)
{
	Py_Initialize();

	return PyModule_Create(&gf2_definition);
}

