import math
import mpqs

# primes less than 212
small_primes = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
    41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211]

# pre-calced sieve of eratosthenes for n = 2, 3, 5, 7
indices = [
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
    107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157,
    163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209]

# distances between sieve values
offsets = [
    10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6,
    6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4,
    2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6,
    4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2]

# tabulated, mod 105
dindices = [
    0, 10, 2, 0, 4, 0, 0, 0, 8, 0, 0, 2, 0, 4, 0,
    0, 6, 2, 0, 4, 0, 0, 4, 6, 0, 0, 6, 0, 0, 2,
    0, 6, 2, 0, 4, 0, 0, 4, 6, 0, 0, 2, 0, 4, 2,
    0, 6, 6, 0, 0, 0, 0, 6, 6, 0, 0, 0, 0, 4, 2,
    0, 6, 2, 0, 4, 0, 0, 4, 6, 0, 0, 2, 0, 6, 2,
    0, 6, 0, 0, 4, 0, 0, 4, 6, 0, 0, 2, 0, 4, 8,
    0, 0, 2, 0, 10, 0, 0, 4, 0, 0, 0, 2, 0, 4, 2]

max_int = 2147483647


# returns the index of x in a sorted list a
# or the index of the next larger item if x is not present
# i.e. the proper insertion point for x in a
def binary_search(a, x):
    s = 0
    e = len(a)
    m = e >> 1
    while m != e:
        if a[m] < x:
            s = m
            m = (s + e + 1) >> 1
        else:
            e = m
            m = (s + e) >> 1
    return m


# divide and conquer list product
def list_prod(a):
    while len(a) > 1:
        a = [m * n for m, n in zip(a[::2], a[1::2] + [1])]
    return a[0]


# greatest common divisor of a and b
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


# extended gcd
def ext_gcd(a, m):
    a = int(a % m)
    x, u = 0, 1
    while a:
        x, u = u, x - (m // a) * u
        m, a = a, m % a
    return m, x, u


# legendre symbol (a|m)
# note: returns m-1 if a is a non-residue, instead of -1
def legendre(a, m):
    return pow(a, (m - 1) >> 1, m)


# modular inverse of a mod m
def mod_inv(a, m):
    return ext_gcd(a, m)[1]


# modular sqrt(n) mod p
# p must be prime
def mod_sqrt(n, p):
    a = n % p
    if p % 4 == 3:
        return pow(a, (p + 1) >> 2, p)
    elif p % 8 == 5:
        v = pow(a << 1, (p - 5) >> 3, p)
        i = ((a * v * v << 1) % p) - 1
        return (a * v * i) % p
    elif p % 8 == 1:
        # Shanks' method
        q = p - 1
        e = 0
        while q & 1 == 0:
            e += 1
            q >>= 1

        n = 2
        while legendre(n, p) != p - 1:
            n += 1

        w = pow(a, q, p)
        x = pow(a, (q + 1) >> 1, p)
        y = pow(n, q, p)
        r = e
        while True:
            if w == 1:
                return x

            v = w
            k = 0
            while v != 1 and k + 1 < r:
                v = (v * v) % p
                k += 1

            if k == 0:
                return x

            d = pow(y, 1 << (r - k - 1), p)
            x = (x * d) % p
            y = (d * d) % p
            w = (w * y) % p
            r = k
    else:  # p == 2
        return a


# integer sqrt of n
def isqrt(n):

    c = (n << 2) // 3
    d = c.bit_length()

    a = d >> 1
    if d & 1:
        x = 1 << a
        y = (x + (n >> a)) >> 1
    else:
        x = (3 << a) >> 2
        y = (x + (c >> a)) >> 1

    if x != y:
        x = y
        y = (x + n // x) >> 1
        while y < x:
            x = y
            y = (x + n // x) >> 1
    return x


# integer cbrt of n
def icbrt(n):
    d = n.bit_length()

    if d % 3 == 2:
        x = 3 << d // 3 - 1
    else:
        x = 1 << d // 3

    y = (2 * x + n // (x * x)) // 3
    if x != y:
        x = y
        y = (2 * x + n // (x * x)) // 3
        while y < x:
            x = y
            y = (2 * x + n // (x * x)) // 3
    return x


# strong probable prime
def is_sprp(n, b=2):
    if n < 2:
        return False
    d = n - 1
    s = 0
    while d & 1 == 0:
        s += 1
        d >>= 1

    x = pow(b, d, n)
    if x == 1 or x == n - 1:
        return True

    for r in range(1, s):
        x = (x * x) % n
        if x == 1:
            return False
        elif x == n - 1:
            return True

    return False


# lucas probable prime
# assumes D = 1 (mod 4), (D|n) = -1
def is_lucas_prp(n, D):

    Q = (1 - D) >> 2

    # n+1 = 2**r*s where s is odd
    s = n + 1
    r = 0
    while s & 1 == 0:
        r += 1
        s >>= 1

    # calculate the bit reversal of (odd) s
    # e.g. 19 (10011) <=> 25 (11001)
    t = 0
    while s:
        if s & 1:
            t += 1
            s -= 1
        else:
            t <<= 1
            s >>= 1

    # use the same bit reversal process to calculate the sth Lucas number
    # keep track of q = Q**n as we go
    U = 0
    V = 2
    q = 1
    # mod_inv(2, n)
    inv_2 = (n + 1) >> 1
    while t:
        if t & 1:
            # U, V of n+1
            U, V = ((U + V) * inv_2) % n, ((D * U + V) * inv_2) % n
            q = (q * Q) % n
            t -= 1
        else:
            # U, V of n*2
            U, V = (U * V) % n, (V * V - 2 * q) % n
            q = (q * q) % n
            t >>= 1

    # double s until we have the 2**r*sth Lucas number
    while r:
        U, V = (U * V) % n, (V * V - 2 * q) % n
        q = (q * q) % n
        r -= 1

    # primality check
    # if n is prime, n divides the n+1st Lucas number, given the assumptions
    return U == 0


# Baillie-PSW #
# this is technically a probabalistic test, but there are no known pseudoprimes
def is_bpsw(n):
    if not is_sprp(n, 2):
        return False

    # idea shamelessly stolen from Mathmatica's PrimeQ
    # if n is a 2-sprp and a 3-sprp, n is necessarily square-free
    if not is_sprp(n, 3):
        return False

    a = 5
    s = 2
    # if n is a perfect square, this will never terminate
    while legendre(a, n) != n - 1:
        s = -s
        a = s - a
    return is_lucas_prp(n, a)


# an 'almost certain' primality check
def is_prime(n):
    if n < 212:
        m = binary_search(small_primes, n)
        return n == small_primes[m]

    for p in small_primes:
        if n % p == 0:
            return False

    # if n is a 32-bit integer, perform full trial division
    if n <= max_int:
        p = 211
        while p * p < n:
            for o in offsets:
                p += o
                if n % p == 0:
                    return False
        return True

    return is_bpsw(n)


# next prime strictly larger than n
def next_prime(n):
    if n < 2:
        return 2

    n += 1
    if n < 212:
        m = binary_search(small_primes, n)
        return small_primes[m]

    # find our position in the sieve rotation via binary search
    x = int(n % 210)
    m = binary_search(indices, x)
    i = int(n + (indices[m] - x))

    # adjust offsets
    offs = offsets[m:] + offsets[:m]
    while True:
        for o in offs:
            if is_prime(i):
                return i
            i += o


# an infinite prime number generator
def primes(start=0):
    for n in small_primes[start:]:
        yield n
    pg = primes(6)
    p = next(pg)
    q = p * p
    sieve = {221: 13, 253: 11}
    n = 211
    while True:
        for o in offsets:
            n += o
            stp = sieve.pop(n, 0)
            if stp:
                nxt = n // stp
                nxt += dindices[nxt % 105]
                while nxt * stp in sieve:
                    nxt += dindices[nxt % 105]
                sieve[nxt * stp] = stp
            elif n < q:
                yield n
            else:
                sieve[q + dindices[p % 105] * p] = p
                p = next(pg)
                q = p * p


# true if n is a prime power > 0
def is_prime_power(n):
    if n > 1:
        for p in small_primes:
            if n % p == 0:
                n = n // p
                while n % p == 0:
                    n = n // p
                return n == 1

        r = isqrt(n)
        if r * r == n:
            return is_prime_power(r)

        s = icbrt(n)
        if s * s * s == n:
            return is_prime_power(s)

        p = 211
        while p * p < r:
            for o in offsets:
                p += o
                if n % p == 0:
                    n = n // p
                    while n % p == 0:
                        n = n // p
                    return n == 1

        if n <= max_int:
            while p * p < n:
                for o in offsets:
                    p += o
                    if n % p == 0:
                        return False
            return True

        return is_bpsw(n)
    return False


max_trial = 1e10
max_pollard = 1e22


def factor(n):
    if n < max_trial:
        return factor_trial(n)
    for p in small_primes:
        if n % p == 0:
            return [p] + factor(n // p)
    if is_prime(n):
        return [n]
    if n < max_pollard:
        p = pollard_rho(n)
    else:
        p = lenstra_ecf(n) or mpqs.mpqs(n)
    return factor(p) + factor(n // p)


def factor_trial(n):
    a = []
    for p in small_primes:
        while n % p == 0:
            a += [p]
            n = n // p
    i = 211
    while i * i < n:
        for o in offsets:
            i += o
            while n % i == 0:
                a += [i]
                n = n // i
    if n > 1:
        a += [n]
    return a


def pollard_rho(n):
    # Brent's variant
    y, r, q = 0, 1, 1
    c, m = 9, 40
    g = 1
    x, ys = y  # syntax highlighting -- these will be set to y anyways
    while g == 1:
        x = y
        for i in range(r):
            y = (y * y + c) % n
        k = 0
        while k < r and g == 1:
            ys = y
            for j in range(min(m, r - k)):
                y = (y * y + c) % n
                q = q * abs(x - y) % n
            g = gcd(q, n)
            k += m
        r *= 2
    if g == n:
        ys = (ys * ys + c) % n
        g = gcd(n, abs(x - ys))
        while g == 1:
            ys = (ys * ys + c) % n
            g = gcd(n, abs(x - ys))
    return g


def ec_add(tup1, tup2, tup3, n):
    x1, z1 = tup1
    x2, z2 = tup2
    x0, z0 = tup3

    t1, t2 = (x1 - z1) * (x2 + z2), (x1 + z1) * (x2 - z2)
    x, z = t1 + t2, t1 - t2
    return z0 * x * x % n, x0 * z * z % n


def ec_double(tup1, tup2, n):
    x, z = tup1
    a, b = tup2

    t1 = x + z
    t1 *= t1
    t2 = x - z
    t2 *= t2
    t3 = t1 - t2
    t4 = 4 * b * t2

    return t1 * t4 % n, t3 * (t4 + a * t3) % n


def ec_multiply(k, p, C, n):
    # Montgomery ladder algorithm
    p0 = p
    q, p = p, ec_double(p, C, n)
    b = k >> 1
    while b > (b & -b):
        b ^= b & -b
    while b:
        if k & b:
            q, p = ec_add(p, q, p0, n), ec_double(p, C, n)
        else:
            q, p = ec_double(q, C, n), ec_add(p, q, p0, n),
        b >>= 1
    return q


def lenstra_ecf(n, m=5):
    # Montgomery curves w/ Suyama parameterization.
    # Based on pseudocode found in:
    # "Implementing the Elliptic Curve Method of Factoring in Reconfigurable Hardware"
    # Gaj, Kris et. al
    # http://www.hyperelliptic.org/tanja/SHARCS/talks06/Gaj.pdf
    # Phase 2 is not implemented.
    B1, B2 = 8, 13
    for i in range(m):
        pg = primes()
        p = next(pg)
        k = 1
        while p < B1:
            k *= p ** int(math.log(B1, p))
            p = next(pg)
        for s in range(B1, B2):
            u, v = s * s - 5, 4 * s
            x = u * u * u
            z = v * v * v
            t = pow(v - u, 3, n)
            P = (x, z)
            C = (t * (3 * u + v) % n, 4 * x * v % n)
            Q = ec_multiply(k, P, C, n)
            g = gcd(Q[1], n)
            if 1 < g < n:
                return g
        B1, B2 = B2, B1 + B2


if __name__ == "__main__":
    import argparse
    import time

    parser = argparse.ArgumentParser(description='Uses a various methods to factor a composite number')
    parser.add_argument('composite', metavar='number_to_factor', type=int, help='the composite number to factor')
    args = parser.parse_args()

    time_begin = time.time()
    print(factor(args.composite))
    print(f"time elapsed: {(time.time() - time_begin):.3f} seconds")
