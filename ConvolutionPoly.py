import numpy as np

class ConvolutionPoly():
    def __init__(self, L, rank=0, char = 0):
        # L is given as [a_0,...,a_n]
        # set rank to 0 for a non-convolution
        # characteristic gives the ring over which we consider the coeffs
        # 0 for Z and q for Z/qZ
        
        def minimal_form(L, rank, char=0):
            if rank > 0:
                reducedL = [sum(L[i::rank]) for i in range(0,rank)]
            else:
                reducedL = L
            if char == 0:
                return reducedL
            else:
                reducedL = np.array(reducedL) % char
                return reducedL.tolist()

        self.coeffs = minimal_form(L, rank, char)

        self.rank = rank

        self.char = char

        def find_degree(L):
            i = 1
            while (L[-i] == 0) & (i < len(L)):
                i+=1
            if (i == len(L)) & (L[0] == 0):
                return -1
            else:
                return len(L) - i
    
        self.degree = find_degree(self.coeffs)
        

    def __getitem__(self, arg):
        return self.coeffs[arg]
    
    def __setitem__(self, arg, value):
        self.coeffs[arg] = value

    # my first idea for printing the polynomial in a nice way
    # TODO: handle negative coefficients
    def __str__(self):
        #create a new list of coefficients as strings with the appropriate x monomial
        newL = []
        for i in range(len(self.coeffs)):
            if not self.coeffs[i] == 0:
                if i == 0:
                    newL.append(str(self.coeffs[i]))
                if i > 0:
                    if self.coeffs[i] == -1:
                        newL.append('-x^'+str(i))
                    elif self.coeffs[i] == 1:
                        newL.append('x^'+str(i))
                    else:
                        newL.append(str(self.coeffs[i])+'x^'+str(i))
        plses = ['+'] * (len(newL) - 1)
        poly = [None]*(len(newL)+(len(newL)-1))
        poly[::2] = newL
        poly[1::2] = plses
        if newL == []:
            return "0"
        else:
            return ''.join(poly)

    # defining addition of convolution polynomials
    def __add__(self, other):
        c1 = self.coeffs
        c2 = other.coeffs
        res = [sum(t) for t in zip(c1,c2)]
        return ConvolutionPoly(res)

    # defining subtraction of convolution polynomials
    # ziplongest probably not needed anymore now that ive standardized the length of coeffs with reduced form. 
    # will want to cast an error when operating on two convolution polynomials of different rank/characteristics
    def __sub__(self, other):
        c1 = self.coeffs
        c2 = other.coeffs
        res = [t1 - t2 for t1, t2 in zip(c1,c2)]
        return ConvolutionPoly(res)
    
    # we will define multiplication using normal polynomial multiplication,
    # then reducedFrom will take care of the rest (itd be interesting to
    # consider which is faster, this way or the direct way)
    def __mul__(self, other):
        r1 = len(self.coeffs)
        r2 = len(other.coeffs)
        result = [0] * (r1+r2 -1)
        for i in range(r1):
            for j in range(r2):
                result[i+j] += self.coeffs[i]*other.coeffs[j]
        return ConvolutionPoly(result,self.rank,self.char)

    def make_normal(self):
        return ConvolutionPoly(self.coeffs,0,self.char)

    #This drops a characteristic 0 polynomial to a characteristic q polynomial
    def reduce_modulus(self, newchar):
        return ConvolutionPoly(self.coeffs,self.rank,newchar)

    # A sort of inverse of reduce_modulus. This lifts a characteristic q polynomial
    # to a characteristic 0 polynomials but with coefficients between -q/2 and q/2
    def center_lift(self):
        L = self.coeffs
        for i in range(len(L)):
            if (L[i] > self.char / 2):
                L[i] = L[i] - self.char
        return ConvolutionPoly(L, self.rank)

    #def inverse(self):         