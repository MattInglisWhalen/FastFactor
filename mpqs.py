import math
import my_math
import time


# Multiple Polynomial Quadratic Sieve
# assumes n is composite
def mpqs(n, verbose=False):

    factor = 1
    time1 = time.time()

    # root_n = my_math.isqrt(n)
    root_2n = my_math.isqrt(n + n)

    # formula chosen by experimentation
    # seems to be close to optimal for n < 10^50
    bound = int(7.5 * math.log(n, 10) ** 2)

    prime = []
    mod_root = []
    log_p = []
    num_prime = 0

    # size of the sieve
    x_max = bound * 5
    x_max_2 = x_max + x_max

    # maximum value on the sieved range
    m_val = (x_max * root_2n) >> 1

    # find a number of small primes for which n is a quadratic residue
    pg = my_math.primes()
    p = next(pg)
    while p < bound or num_prime < 3:

        # legendre (n|p) is only defined for odd p
        if p > 2:
            leg = my_math.legendre(n, p)
        else:
            leg = n & 1

        if leg == 1:
            prime += [p]
            log_p += [math.log(p, 10)]
            r = int(my_math.mod_sqrt(n, p))
            roots = [r]
            q = p
            while q < x_max:
                # find all square roots mod p^n via Hensel Lifting
                r = int((r + (n - r * r) * my_math.mod_inv(r + r, q)) % q)
                # assert r*r%q == n%q
                roots += [r]
                q *= p
            mod_root += [roots]
            num_prime += 1
        elif leg == 0:
            if verbose:
                print('trial division found factors:')
                print(p, 'x', n // p)
            return p

        p = next(pg)

    # fudging the threshold down a bit makes it easier to find powers of partial-partial
    # relationships, but it also makes the smoothness check slower. reducing by twice the log
    # of the largest prime in the factor base results in cofactors less than that value squared
    thresh = math.log(m_val, 10) - log_p[-1] * 2

    # skip small primes. they contribute very little to the log sum
    # and add a lot of unnecessary entries to the table
    # instead, fudge the threshold down a bit, according to expected number of factors
    min_prime = int(thresh * 2)
    sp_idx = my_math.binary_search(prime, min_prime)
    sieve_primes = prime[sp_idx:]

    fudge = sum([log_p[i] // (prime[i] - 1) for i in range(sp_idx)])

    sums = [fudge] * x_max_2

    if verbose:
        print('smoothness bound:', bound)
        print('sieve size:', x_max)
        print('log threshold:', thresh)
        print('skipping primes less than:', min_prime)

    smooth = []
    used_prime = set()
    partial = {}
    num_smooth = 0

    num_partial = 0
    num_poly = 0
    root_A = my_math.isqrt(root_2n // x_max)

    if verbose:
        print('sieving for smooths...')
    while True:
        # find an integer value A such that:
        # A is =~ sqrt(2*n) / x_max
        # A is a perfect square
        # sqrt(A) is prime, and n is a quadratic residue mod sqrt(A)
        while True:
            root_A = my_math.next_prime(root_A)
            leg = my_math.legendre(n, root_A)
            if leg == 1:
                break
            elif leg == 0:
                if verbose:
                    print('dumb luck found factors:')
                    print(root_A, 'x', n // root_A)
                return root_A

        A = root_A * root_A

        # solve for an adequate B
        # B*B is a quadratic residue mod n, such that B*B-A*C = n
        # this is unsolvable if n is not a quadratic residue mod sqrt(A)
        b = my_math.mod_sqrt(n, root_A)
        B = (b + (n - b * b) * my_math.mod_inv(b + b, root_A)) % A

        # B*B-A*C = n <=> C = (B*B-n)/A
        C = (B * B - n) // A

        num_poly += 1

        # sieve for prime factors
        i = sp_idx
        for p in sieve_primes:
            logp = log_p[i]

            e = 0
            q = p
            while q < x_max:
                inv_A = my_math.mod_inv(A, q)
                # modular root of the quadratic
                a = int(((mod_root[i][e] - B) * inv_A) % q)
                b = int(((q - mod_root[i][e] - B) * inv_A) % q)

                amx = a + x_max
                bmx = b + x_max

                apx = amx - q
                bpx = bmx - q

                k = q
                while k < x_max:
                    sums[apx + k] += logp
                    sums[bpx + k] += logp
                    sums[amx - k] += logp
                    sums[bmx - k] += logp
                    k += q

                q *= p
                e += 1

            i += 1

        # check for smooths
        x = -x_max
        i = 0
        while i < x_max_2:
            v = sums[i]
            if v > thresh:
                vec = set()
                sqr = []
                # because B*B-n = A*C
                # (A*x+B)^2 - n = A*A*x*x+2*A*B*x + B*B - n
                #               = A*(A*x*x+2*B*x+C)
                # gives the congruency
                # (A*x+B)^2 = A*(A*x*x+2*B*x+C) (mod n)
                # because A is chosen to be square, it doesn't need to be sieved
                sieve_val = (A * x + B + B) * x + C

                if sieve_val < 0:
                    vec = {-1}
                    sieve_val = -sieve_val

                for p in prime:
                    while sieve_val % p == 0:
                        if p in vec:
                            # keep track of perfect square factors
                            # to avoid taking the sqrt of a gigantic number at the end
                            sqr += [p]
                        vec ^= {p}
                        sieve_val = int(sieve_val // p)

                if sieve_val == 1:
                    # smooth
                    smooth += [(vec, (sqr, (A * x + B), root_A))]
                    used_prime |= vec
                elif sieve_val in partial:
                    # combine two partials to make a (xor) smooth
                    # that is, every prime factor with an odd power is in our factor base
                    pair_vec, pair_vals = partial[sieve_val]
                    sqr += list(vec & pair_vec) + [sieve_val]
                    vec ^= pair_vec
                    smooth += [(vec, (sqr + pair_vals[0], (A * x + B) * pair_vals[1], root_A * pair_vals[2]))]
                    used_prime |= vec
                    num_partial += 1
                else:
                    # save partial for later pairing
                    partial[sieve_val] = (vec, (sqr, A * x + B, root_A))
            x += 1

            # reset the value for the next go
            sums[i] = fudge
            i += 1

        prev_num_smooth = num_smooth
        num_smooth = len(smooth)
        num_used_prime = len(used_prime)

        if verbose:
            print(100 * num_smooth // num_prime, 'percent complete\r')

        if num_smooth > num_used_prime and num_smooth > prev_num_smooth:
            if verbose:
                print(f"{num_poly} polynomials sieved ({num_poly*x_max_2} values)")
                print(f"found {num_smooth} smooths ({num_partial} from partials)"
                      f" in {time.time()-time1:.3f} seconds")
                print('solving for non-trivial congruencies...')

            # set up bit fields for gaussian elimination
            masks = []
            mask = 1
            bit_fields = [0] * num_used_prime
            for vec, vals in smooth:
                masks += [mask]
                i = 0
                for p in used_prime:
                    if p in vec:
                        bit_fields[i] |= mask
                    i += 1
                mask += mask

            # row echelon form
            col_offset = 0
            null_cols = []
            for col in range(num_smooth):
                pivot = col - col_offset == num_used_prime or bit_fields[col - col_offset] & masks[col] == 0
                for row in range(col + 1 - col_offset, num_used_prime):
                    if bit_fields[row] & masks[col]:
                        if pivot:
                            bit_fields[col - col_offset], bit_fields[row] = bit_fields[row], bit_fields[
                                col - col_offset]
                            pivot = False
                        else:
                            bit_fields[row] ^= bit_fields[col - col_offset]
                if pivot:
                    null_cols += [col]
                    col_offset += 1

            # reduced row echelon form
            for row in range(num_used_prime):
                # lowest set bit
                mask = bit_fields[row] & -bit_fields[row]
                for up_row in range(row):
                    if bit_fields[up_row] & mask:
                        bit_fields[up_row] ^= bit_fields[row]

            # check for non-trivial congruencies
            for col in null_cols:
                all_vec, (lh, rh, rA) = smooth[col]
                lhs = lh  # sieved values (left hand side)
                rhs = [rh]  # sieved values - n (right hand side)
                rAs = [rA]  # root_As (cofactor of lhs)
                i = 0
                for field in bit_fields:
                    if field & masks[col]:
                        vec, (lh, rh, rA) = smooth[i]
                        lhs += list(all_vec & vec) + lh
                        all_vec ^= vec
                        rhs += [rh]
                        rAs += [rA]
                    i += 1

                factor = my_math.gcd(my_math.list_prod(rAs) * my_math.list_prod(lhs) - my_math.list_prod(rhs), n)
                if 1 < factor < n:
                    break
            else:
                if verbose:
                    print('none found.')
                continue
            break

    if verbose:
        print('factors found:')
        print(factor, 'x', n // factor)
        print(f"time elapsed: {(time.time() - time1):.3f} seconds")
    return factor


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Uses a MPQS to factor a composite number')
    parser.add_argument('composite', metavar='number_to_factor', type=int, help='the composite number to factor')
    parser.add_argument('--verbose', dest='verbose', action='store_true', help="enable verbose output")
    args = parser.parse_args()

    if args.verbose:
        mpqs(args.composite, args.verbose)
    else:
        time_begin = time.time()
        print(mpqs(args.composite))
        print(f"time elapsed: {(time.time() - time_begin):.3f} seconds")
