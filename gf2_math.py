"""
https://bitbucket.org/jason_s/libgf2/src/fc688f8678b5c876603e034b7a1e20651dc9f4e2/src/libgf2/gf2.py?at=default&fileviewer=file-view-default
Objects and functions for working with polynomials in GF(2).

Copyright 2013-2017 Jason M. Sachs

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

# ---------- Basic arithmetic     

def _gf2mod(a,m):
    """
    Computes `a` mod `m`.

    Parameters
    ----------
    a, m : integer
        Polynomial coefficient bit vectors.

    Returns
    -------
    integer
        Polynomial coefficient bit vectors of `a` mod `m`.    
    """
    m2 = m
    i = 0
    while m2 < a:
        m2 <<= 1
        i += 1
    while i >= 0:
        anew = a ^ m2
        if anew < a:
            a = anew
        m2 >>= 1
        i -= 1
    return a    

def _gf2mul(a,b):
    """
    Computes ``a * b``.

    Parameters
    ----------
    a, b : integer
        Polynomial coefficient bit vectors.

    Returns
    -------
    c : integer
        Polynomial coefficient bit vector of ``c = a * b``.    
    """

    c = 0
    while a > 0:
        if (a & 1) != 0:
            c ^= b
        b <<= 1
        a >>= 1
    return c

def _gf2bitlength_linear(a):
    """
    Computes the length of a polynomial coefficient bit vector = degree + 1. 

    Parameters
    ----------
    a : integer
        Polynomial coefficient bit vector.

    Returns
    -------
    n : integer
        length of  polynomial `a`.    
    """     
    n = 0
    while a > 0:
        n += 1
        a >>= 1
    return n

def _gf2bitlength(a):
    """
    Computes the length of a polynomial coefficient bit vector = degree + 1. 

    Parameters
    ----------
    a : integer
        Polynomial coefficient bit vector.

    Returns
    -------
    n : integer
        length of  polynomial `a`.    
    """     
    if a == 0:
        return 0
    n = 1
    while (a >> n) > 0:
        n <<= 1
    n >>= 1
    b = n
    while b > 0:
        b >>= 1
        nnew = n ^ b
        if (a >> nnew) > 0:
            n = nnew
    return n + 1

def _gf2mulmod(a,b,m):
    """
    Computes ``a * b mod m``.
    *NOTE*: Does *not* check whether `a` and `b` are both
    smaller in degree than `m`.

    Parameters
    ----------
    a, b, m : integer
        Polynomial coefficient bit vectors.
        Polynomials `a` and `b` should be smaller degree than `m`.

    Returns
    -------
    c : integer
        Polynomial coefficient bit vector of ``c = a * b mod m``.    
    """
    c = 0
    while a > 0:
        if (a & 1) != 0:
            c ^= b
        b <<= 1
        b2 = b ^ m
        if b2 < b:
            b = b2
        a >>= 1
    return c

def _gf2mulxmod(a,m):
    """
    Computes ``a * x mod m``.
    *NOTE*: Does *not* check whether `a` is smaller in degree than `m`.

    Parameters
    ----------
    a, m : integer
        Polynomial coefficient bit vectors.
        Polynomial `a` should be smaller degree than `m`.

    Returns
    -------
    c : integer
        Polynomial coefficient bit vector of ``c = a * x mod m``.    
    """
    c = a << 1
    c2 = c^m
    if c2 < c:
        c = c2
    return c

def _gf2mulxinvmod(a,m):
    """
    Computes ``a * x^(-1) mod m``.
    *NOTE*: Does *not* check whether `a` is smaller in degree than `m`.

    Parameters
    ----------
    a, m : integer
        Polynomial coefficient bit vectors.
        Polynomial `a` should be smaller degree than `m`.

    Returns
    -------
    c : integer
        Polynomial coefficient bit vector of ``c = a * x^(-1) mod m``.    
    """
    c = (a ^ ((a&1)*m)) >> 1
    return c             

def _gf2power(a,k):
    """
    Computes `a` to the `k`th power.

    Parameters
    ----------
    a : integer
        Polynomial coefficient bit vector.
    k : integer     
        Exponent.

    Returns
    -------
    c : integer
        Polynomial coefficient bit vector of `a` raised to the `k`th power.
    """    
    c = 1
    while k > 0:
        if (k & 1) != 0:
            c = _gf2mul(c,a)
        a = _gf2mul(a,a)
        k >>= 1
    return c

def _gf2powmod(a, k, m):
    """
    Computes `a` to the `k`th power, in the quotient ring GF(2)[x]/m(x).
    *NOTE*: Does *not* check whether k is nonnegative.

    Parameters
    ----------
    a : integer
        Polynomial coefficient bit vector.
    k : integer
        Nonnegative exponent.       
    m : integer
        Polynomial coefficient bit vector. 

    Returns
    -------
    c : integer 
        Polynomial coefficient bit vector of `a` raised to the `k`th power, 
        mod `m`. 
    """
    c = 1
    while k > 0:
        if (k & 1) != 0:
            c = _gf2mulmod(c,a,m)
        a = _gf2mulmod(a,a,m)
        k >>= 1
    return c 

def _gf2divmodvect(avec,dvec):
    """
    Computes `q` = `a`/`d` and `r` = `a`%`d` in a vectorized form.

    Parameters
    ----------
    avec : array_like
        Array of dividends, each of which should be a polynomial coefficient bit vector
    dvec : array_like
        Array of divisors, each of which should be a polynomial coefficient bit vector     

    Returns
    -------
    q : array_like
        Array of quotients, each of which is a polynomial coefficient bit vector
    r : array_like 
        Array of quotients, each of which is a polynomial coefficient bit vector
    """
    na = _gf2bitlength(avec[0])
    nd = _gf2bitlength(dvec[0])
    i = na - nd
    q = 0
    test = 1 << (na-1)
    while i >= 0:
        if (avec[0] & test) != 0:
            avec = [a ^ (d << i) for (a,d) in zip(avec,dvec)]
            q |= (1 << i)
        i -= 1
        test >>= 1
    r = avec    
    return (q,r)

def _gf2exteuc(a,b):
    r"""
    Computes the extended Euclidean algorithm using Blankinship's method.

    Returns :math:`(g,u,v)` such that :math:`g = \gcd(a,b)` and :math:`g = au+bv` in :math:`GF(2)`.

    Parameters
    ----------
    a, b : integer
        Polynomial coefficient bit vectors.

    Returns
    -------
    g : integer
        Polynomial coefficient bit vector :math:`g = \gcd(a,b)`
    u, v : integer    
        Polynomial coefficient bit vectors where :math:`g = au+bv`
    """
    arow = [a,1,0]
    brow = [b,0,1]
    while True:
        (_,rrow) = _gf2divmodvect(arow, brow)
        if rrow[0] == 0:
            break
        arow = brow
        brow = rrow
    return tuple(brow)

def _exteuc(a,b):
    """
    based on Blankenship's algorithm
    return (g,u,v) such that g = gcd(a,b) and g = au+bv
    """
    arow = [a,1,0]
    brow = [b,0,1]
    while True:
        (q,r) = divmod(arow[0],brow[0])
        if r == 0:
            break
        rrow = [r,arow[1]-q*brow[1],arow[2]-q*brow[2]]
        arow = brow
        brow = rrow
    return tuple(brow)

def _gf2modinv(x,m):
    #(r,y,_) = _exteuc(x,m)
    (r,y,_) = _gf2exteuc(x,m)
    if r != 1:
        raise ValueError('%d and %d are not relatively prime but have a common factor of %d' % (x,m,r))
    return y

def _modinv(x,m):
    (r,y,_) = _exteuc(x,m)
    if r != 1:
        raise ValueError('%d and %d are not relatively prime but have a common factor of %d' % (x,m,r))
    return y