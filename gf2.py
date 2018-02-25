from functools import reduce
from collections import defaultdict
import math

try:
    from gf2c import gf2CheckIrreduciblePolynomial as checkIrreduciblePolynomialC
    from gf2c import gf2CheckIrreduciblePolynomialRabin as checkIrreduciblePolynomialRabinC
    from gf2c import gf2CheckGeneratorPolynomial as checkGeneratorPolynomialC
except ImportError:
    pass

try:
    from gf2c import _gf2powmod as gf2powmod
except ImportError:
    from gf2_math import _gf2powmod as gf2powmod
    
try:
    from gf2c import _gf2mulmod as gf2mulmod
except ImportError:
    from gf2_math import _gf2mulmod as gf2mulmod
    
try:
    from gf2c import _gf2mod as gf2mod
except ImportError:
    from gf2_math import _gf2mod as gf2mod
    
try:
    from gf2c import _gf2modinv as gf2modinv
except ImportError:
    from gf2_math import _gf2modinv as gf2modinv
    
try:
    from gf2c import _gf2gcd as gf2gcd
except ImportError:
    pass

class InverseException(Exception):
    pass   

class GF2:
    def __init__(self, *, degree = 0, coefficients = None, size = 8, value = None, mod=None):
        self.mod = mod
        if self.mod is None:
            self.mod = 0x11B
        
        if isinstance(value, GF2):
            value = value.value
            
        self.size = size
        if coefficients is None and value is None:
            self.value = 1 << degree
        
        elif coefficients is not None:
            self.value = reduce(lambda x, y: (x << 1) + y, coefficients)
                
        elif value is not None:
            self.value = value % 2**size
            
        self._coefficients = None
    
    @property
    def coefficients(self):
        if self._coefficients is None:
            self._coefficients = tuple([(self.value >> i) & 1 for i in range(self.value.bit_length())])
            if len(self._coefficients) == 0:
                self._coefficients = (0,)
                
        return self._coefficients
    
    def __repr__(self):
        return bin(self.value)
    
    def __invert__(self):
        return gf2modinv(self.value, self.mod)
    
    def __xor__(self, other):
        if isinstance(other, GF2):
            return GF2(value = self.value ^ other.value, size = self.size, mod = self.mod)
        elif isinstance(other, int):
            return self ^ GF2(value = other, size = self.size, mod = self.mod)
        else:
            return NotImplemented

    def __rxor__(self, other):
        return self ^ other
    def __ixor__(self, other):
        return self ^ other    
        
    def __and__(self, other):
        if isinstance(other, GF2):
            return GF2(value = self.value & other.value, size = self.size, mod = self.mod)
        elif isinstance(other, int):
            return self & GF2(value = other, size = self.size, mod = self.mod)        
        else:
            return NotImplemented
        
    def __rand__(self, other):
        return self & other
    def __iand__(self, other):
        return self & other    
        
    def __or__(self, other):
        if isinstance(other, GF2):
            return GF2(value = self.value | other.value, size = self.size)
        elif isinstance(other, int):
            return self | GF2(value = other, size = self.size, mod = self.mod)        
        else:
            return NotImplemented  
    
    def __ror__(self, other):
        return self | other
    def __ior__(self, other):
        return self | other     
        
    def __eq__(self, other):
        if isinstance(other, GF2):
            return self.value == other.value and self.size == other.size and self.mod == other.mod
        if isinstance(other, int):
            return self.value == other
        else:
            raise NotImplemented
        
    def __hash__(self):
        return hash((self.value, self.size, self.mod))
    
    def __add__(self, other):
        if isinstance(other, GF2):
            return self ^ other
        elif isinstance(other, int):
            return self ^ GF2(value = other, size = self.size, mod = self.mod)          
        else:
            return NotImplemented
        
    def __radd__(self, other):
        return self + other
    def __iadd__(self, other):
        return self + other    
        
    def __sub__(self, other):
        if isinstance(other, GF2):
            return self ^ other
        elif isinstance(other, int):
            return self - GF2(value = other, size = self.size, mod = self.mod)          
        else:
            return NotImplemented
        
    def __rsub__(self, other):
        return self - other
    def __isub__(self, other):
        return self - other     
    
    def __mul__(self, other):
        if isinstance(other, GF2):
            p = gf2mulmod(self.value, other.value, self.mod)
            return GF2(value = p, size = self.size, mod = self.mod)
        elif isinstance(other, int):
            return self * GF2(value = other, size = self.size, mod = self.mod)            
        else:
            return NotImplemented  
        
    def __rmul__(self, other):
        return self.__mul__(other)
    def __imul__(self, other):
        return self.__mul__(other) 
        
    def __truediv__(self, other):
        if isinstance(other, GF2):
            if self == other:
                return GF2(value = 1, size = self.size, mod = self.mod)
            elif other == GF2(value = 0, size = self.size, mod = self.mod):
                return GF2(value = 0, size = self.size, mod = self.mod)
            else:
                return self.__mul__(~other)
        elif isinstance(other, int):
            return self / GF2(value = other, size = self.size, mod = self.mod)                    
        else:
            return NotImplemented
        
    def __rtruediv__(self, other):
        return (~self).__mul__(other)
    def __itruediv__(self, other):
        return self.__truediv__(other)
    
    def __mod__(self, other):
        if isinstance(other, GF2):
            return GF2(value = gf2mod(self.value, other.value), size = self.size, mod = self.mod)
        elif isinstance(other, int):
            return GF2(value = self.value % other, size = self.size, mod = self.mod)
        else:
            return NotImplemented
        
    def __imod__(self, other):
        return self % other    
    
    def __divmod__(self, other):
        if isinstance(other, GF2):
            return (self.__truediv__(other), self.__mod__(other))
        elif isinstance(other, int):
            return divmod(self, GF2(value = self.value, size = self.size, mod = self.mod))
        else:
            raise NotImplemented()
        
    def __pow__(self, other, mod = None):
        if mod is None:
            mod = self.mod        
        GF2constructor = lambda x: GF2(value = x, size = self.size, mod = self.mod)
        
        try:
            other = int(other)
        except TypeError:
            pass
        
        if isinstance(other, int):

            return GF2constructor(gf2powmod(self.value, other, mod))
        else:
            raise NotImplemented
        
    def __ipow__(self, other):
        return self ** other
        
    def __lshift__(self, other):
        if isinstance(other, int):
            return GF2(value = self.value << other, size = self.size, mod = self.mod)
        else:
            return NotImplemented
        
    def __ilshift__(self, other):
        return self << other
        
    def __rshift__(self, other):
        if isinstance(other, int):
            return GF2(value = self.value >> other, size = self.size, mod = self.mod)
        else:
            return NotImplemented 
        
    def __irshift__(self, other):
        return self << other    
        
    def __bool__(self):
        return bool(self.value)
    
    def __lt__(self, other):
        if isinstance(other, GF2):
            return self.value < other.value
        elif isinstance(other, int):
            return self.value < other
        else:
            return NotImplemented  
        
    def __gt__(self, other):
        if isinstance(other, GF2):
            return self.value > other.value
        elif isinstance(other, int):
            return self.value > other
        else:
            return NotImplemented 
        
    def __le__(self, other):
        if isinstance(other, GF2):
            return self.value <= other.value
        elif isinstance(other, int):
            return self.value <= other
        else:
            return NotImplemented 
        
    def __ge__(self, other):
        if isinstance(other, GF2):
            return self.value >= other.value
        elif isinstance(other, int):
            return self.value >= other
        else:
            return NotImplemented
        
    def __int__(self):
        return self.value
    
    def __index__(self):
        return self.value
    
    def __float__(self):
        return float(self.value)
    
    def __abs__(self):
        return self
    
    def __neg__(self):
        return self
    
    def bit_length(self):
        return self.value.bit_length()
    
    def to_bytes(self, *args, **kwargs):
        return self.value.to_bytes(*args, **kwargs)

_irreduciblePolynomialsCache = {}
def checkIrreduciblePolynomial(poly, size):
    
    if (poly, size) in _irreduciblePolynomialsCache:
        return _irreduciblePolynomialsCache[(poly, size)]
    
    for i in range(2, 2**(size-1)):
        if gf2mod(poly, i) == 0:
            _irreduciblePolynomialsCache[(poly, size)] = False
            return False            
        
    _irreduciblePolynomialsCache[(poly, size)] = True
    return True

def checkIrreduciblePolynomialsRabin(poly, size):
    h = gf2powmod(0b10, 2**(size//2), poly)
    if gf2gcd(poly, h ^ 0b10) != 1:
        return False
    
    h = gf2powmod(0b10, 2**(size), poly)
    if gf2mod(h, poly) != 0b10:
        return False
    
    return True

def findIrreduciblePolynomials(size, *, start = None, stop = None):
    if start is None:
        start = 2**size + 1
    if stop is None:
        stop = 2**(size+1)
        
    try:
        checker = checkIrreduciblePolynomialC
    except NameError:
        checker = checkIrreduciblePolynomial
        
    for mod in range(start, stop):
        if checker(mod, size):
            yield mod

_randomIrreduciblePolynomialsCache = []
def findRandomIrreduciblePolynomial(size, random, *, cacheFile = 'randomIrreduciblePolynomialsCache.txt'):
    # 30 of size 8
    # 4080 of size 16
    # 698870 of size 24
    # 134215680 of size 32
    # 27487764474 of size 40      
    #knownCounts = {8:30, 16:4080, 24:698870}
    #knownCounts = {8:30, 16:4080}
    knownCounts = {}

    if size in knownCounts:
        polyNum = random.randrange(0, knownCounts[size])
        global _randomIrreduciblePolynomialsCache
        if len(_randomIrreduciblePolynomialsCache) > polyNum:
                return _randomIrreduciblePolynomialsCache[polyNum]
        
        import os
        from packedInt import bytesToNum, numToBytes
        mode = 'rb+'
        if not os.path.exists(cacheFile):
            mode = 'wb+'
        with open(cacheFile, mode, buffering=0) as f:
            _randomIrreduciblePolynomialsCache = []
            seek = 0
            for k,v in knownCounts.items():
                if k < size:
                    seek += k//8*v
            f.seek(seek, 0)
            for i in range(knownCounts[size]):
                data = f.read(size//8)
                if len(data) < size//8:
                    f.seek(-len(data), 1)
                    break
                
                value = bytesToNum(b'\x01' + data)
                if value == 0:
                    f.seek(-len(data), 1)
                    break                    
                _randomIrreduciblePolynomialsCache.append(value)
            
            if len(_randomIrreduciblePolynomialsCache) > polyNum:
                return _randomIrreduciblePolynomialsCache[polyNum]
            
            if len(_randomIrreduciblePolynomialsCache) == 0:
                start = None
            else:
                start = _randomIrreduciblePolynomialsCache[-1] + 1
                
            for poly in findIrreduciblePolynomials(size, start=start):
                _randomIrreduciblePolynomialsCache.append(poly)
                data = numToBytes(poly)
                f.write(data[1:])
                try:
                    return _randomIrreduciblePolynomialsCache[polyNum]
                except IndexError:
                    pass
            assert False       
            
    else:
        try:
            checker = checkIrreduciblePolynomialC
        except NameError:
            checker = checkIrreduciblePolynomial   
        
        #checker = checkIrreduciblePolynomialsRabin
        checker = checkIrreduciblePolynomialRabinC
        
        getPoly = lambda: random.randrange(2**size + 1, 2**(size+1))
        poly = getPoly()
        while not checker(poly, size):
            assert not checkIrreduciblePolynomialsRabin(poly, size)
            poly = getPoly()
        return poly

_generatorPolynomialCache = defaultdict(set)
def checkGeneratorPolynomial(poly, mod, size):
    poly = int(poly)
    if (poly, size, mod) in _generatorPolynomialCache:
        return _generatorPolynomialCache[(poly, size, mod)]
    
    count = bin(poly).count("1")
    if count <= 1:
        return False
        
    binPows = [poly]
    for i in range(1, size):
        binPows.append(gf2mulmod(binPows[-1], binPows[-1], mod))
    
    for j in _generatorPolynomialCache[(size, mod)]:
        if gf2powmod(poly, j, mod) == 1:
            _generatorPolynomialCache[(poly, size, mod)] = False
            return False 
    
    for j in range(2, 2**size-1):
        binJ = [i for i in reversed(bin(j)[2:])]
        vals = [p for i, p in enumerate(binPows) if len(binJ) > i and binJ[i] == '1']
        prod = reduce(lambda x, y: gf2mulmod(x,y,mod), vals, 1)
        assert prod == gf2powmod(poly, j, mod)
        if prod == 1:
            _generatorPolynomialCache[(poly, size, mod)] = False
            _generatorPolynomialCache[(size, mod)].add(j)
            return False
    _generatorPolynomialCache[(poly, size, mod)] = True
    return True

def findGeneratorPolynomials(size, mod, *, start = None, stop = None):
    if start is None:
        start = 2
    if stop is None:
        stop = 2**size
        
    try:
        checker = checkGeneratorPolynomialC       
    except NameError:  
        checker = checkGeneratorPolynomial
        
    for i in range(start, stop):
        if checker(i, mod, size):
            yield GF2(value=i, size=size, mod=mod)
        

def findRandomGeneratorPolynomial(size, mod, random):
    try:
        checker = checkGeneratorPolynomialC       
    except NameError:  
        checker = checkGeneratorPolynomial
    
    poly = random.randrange(2, 2**size-1)
    while not checker(poly, mod, size):
        poly = random.randrange(2, 2**size-1)
    return GF2(value=poly, size=size, mod=mod)

def cacheIrreducibleAndGeneratorPolynomials(sizes, filename):
    with open(filename, 'a', buffering=1) as f:
        for size in sizes:
            f.write('Size: %s\n' % size)
            for mod in findIrreduciblePolynomials(size):
                f.write('%s: ' % hex(mod)[2:])
                for gen in findGeneratorPolynomials(size, mod):
                    f.write('{0:>0{1}} '.format(hex(gen)[2:], math.ceil(size/4)))
                    f.flush()
                f.write('\n')
            f.write('\n')


if __name__ == '__main__':
    import random
    
    lgSize = 32
    
    for i in range(1000):
        mod = findRandomIrreduciblePolynomial(lgSize, random)
        gen = findRandomGeneratorPolynomial(lgSize, mod, random)
        
        print(i, mod, gen)