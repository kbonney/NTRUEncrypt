# DEPRICATED CODE


# The multiplicitive inverse of a in the integers modulo p.
# Return d s.t.
# a * d == 1 mod p
# NOTE: computationally brute force, could contribute a lot of inefficiency down
# the line
def invmodp(a, p):
    for d in range(1, p):
        r = (d * a) % p
        if r == 1:
            break
    else:
        raise ValueError('%d has no inverse mod %d' % (a, p))
    return d

# returns k,r such that a = k*b + r with deg r < deg b or r = 0
# global poly_div_alg
def poly_div_alg(a,b):
    # a,b convolution polys
    # b nonzero
    p = a.char
    a = a.make_normal()
    b = b.make_normal()
    k = ConvolutionPoly([0])
    r = a

    def adjust_values(k,r,p):
        adjustlist = [0]*(r.degree - b.degree)
        adjustlist[r.degree - b.degree - 1] = r.coeffs[r.degree - 1] * \
        invmodp(b.coeffs[b.degree - 1],p) 
        adjust = ConvolutionPoly(adjustlist)
        k = k + adjust
        r = r + (adjust * b)
        return (k,r)
    
    #recursively drop the degree of r using adjust_values until deg r < deg b
    while (r.degree >= b.degree):
        k,r = adjust_values(k,r, p)
    
    return (k,r)

# returns the gcd of polynomial a and b, r
# and a history of r and s for use in extended alg
def poly_eucl_alg(a,b):
    r = a
    s = b
    pairs = [(r,s)]
    while s.degree != -1:
        r, s = poly_div_alg(r,s)
        pairs.append((r,s))
    return r, pairs