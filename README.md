# FastFactor
A fast implementation of the Multiple Polynomial Quadratic Sieve. I saw 
[this post](https://codegolf.stackexchange.com/questions/8629/fastest-semiprime-factorization) and wanted to preserve it.

### Basic Concepts

In general, a quadratic sieve is based on the following observation: any odd composite $n$ may be represented as:

$n=(x+d)(x−d)=x^2−d^2⇒d^2=x^2−n$

This is not very difficult to confirm. Since $n$ is odd, the distance between any two cofactors of $n$ must be even $2d$, where $x$ is 
the mid point between them. Moreover, the same relation holds for any multiple of $n$.

$abn=(ax+ad)(bx−bd)=abx^2−abd^2⇒abd^2=abx^2−abn$

Note that if any such $x$ and $d$ can be found, it will immediately result in a (not necessarily prime) factor of $n$, since $x+d$ and $x−d$
both divide $n$ by definition. This relation can be further weakened - at the consequence of allowing potential trivial congruencies - to the following 
form:

$d^2≡x^2\mod n$

So in general, if we can find two perfect squares which are equivalent mod $n$, then it's fairly likely that we can directly produce a factor of $n$
a la $\gcd(x±d,n)$. Seems pretty simple, right?

Except it's not. If we intended to conduct an exhaustive search over all possible xxx, we would need to search the entire range from 
$\sqrt{n}$ $\sqrt{2n}$, which is marginally smaller than full trial division, but also requires an expensive `is_square` operation each iteration 
to confirm the value of $d$. Unless it is known beforehand that $n$ has factors very near $\sqrt{n}$, trial division is likely to be faster.

Perhaps we can weaken this relation even more. Suppose we chose an $x$, such that for

$y≡x^2\mod n$

a full prime factorization of $y$ is readily known. If we had enough such relations, we should be able to _construct_ an adequate $d$, if we choose 
a number of $y$ such that their product is a perfect square; that is, all prime factors are used an even number of times. 
In fact, if we have more such $y$ than the total number of unique prime factors they contain, a solution is guaranteed to exist; 
it becomes a system of linear equations. The question now becomes, how do we chose such $x$? That's where sieving comes into play.

### The Sieve

Consider the polynomial:

$y(x)=x^2−n$

Then for any prime $p$ and integer $k$, the following is true:

$y(x+kp)=(x+kp)^2−n = x^2+2xkp+(kp)^2−n = y(x)+2xkp+(kp)^2 ≡ y(x)\mod p$


This means that after solving for the roots of the polynomial mod $p$ - that is, you've found an $x$ such that $y(x)≡0 \mod p$, ergo $y$ is 
divisible by $p$ - then you have found an infinite number of such $x$. In this way, you can sieve over a range of $x$, identifying small 
prime factors of $y$, hopefully finding some for which all prime factors are small. 
Such numbers known as $k−smooth$, where $k$ is the largest prime factor used.

There's a few problems with this approach, though. Not all values of $x$ are adequate, in fact, there's only very few of them, 
centered around $\sqrt{n}$. Smaller values will become largely negative (due to the $-n$ term), and larger values will become too large, 
such that it is unlikely that their prime factorization consists only of small primes. There will be a number of such $x$, but unless the 
composite you're factoring is very small, it's highly unlikely that you'll find enough smooths to result in a factorization. And so, for 
larger $n$, it becomes necessary to sieve over _multiple_ polynomials of a given form.

### Multiple Polynomials

So we need more polynomials to sieve? How about this:

$y(x)=(Ax+B)^2−n$

That'll work. Note that $A$ and $B$ could literally be any integer value, and the math still holds. All we need to do is choose a few 
random values, solve for the root of the polynomial, and sieve the values close to zero. At this point we could just call it 
good enough: if you throw enough stones in random directions, you're bound to break a window sooner or later.

Except, there's a problem with that too. If the slope of the polynomial is large at the x-intercept, there'll only be a few suitable 
values to sieve per polynomial. It'll work, but you'll end up sieving a whole lot of polynomials before you get what you need. Can we do better?

We can do better. An observation, as a result of [Montgomery](http://en.wikipedia.org/wiki/Peter_Montgomery_(mathematician)) is as 
follows: if $A$ and $B$ are chosen such that there exists some $C$ satisfying

$B^2−n\=AC$

then the entire polynomial can be rewritten as

$y(x)\=(Ax+B)^2−n = (Ax)^2+2ABx+B2−n = A(Ax2+2Bx+C)$

Furthermore, if $A$ is chosen to be a perfect square, the leading $A$ term can be neglected while sieving, resulting in much smaller 
values, and a much flatter curve. For such a solution to exist, $n$ must be a [quadratic residue](http://en.wikipedia.org/wiki/Quadratic_residue) 
mod $\sqrt{A}$, which can be known immediately by computing the 
[Legendre symbol](http://en.wikipedia.org/wiki/Legendre_symbol): $(n|\sqrt{A})=1$. Note that in order to solve for $B$, a complete prime factorization of 
$\sqrt{A}$ needs to be known (in order to take the modular square root $\sqrt{n} \mod \sqrt{A}$), which is why $\sqrt{A}$ is typically chosen to be prime.

It can then be shown that if $A\approx\sqrt{2n}/M$, then for all values of $x \in [−M,M]$:

$|y(x)|≤M\sqrt{n/2}$

And now, finally, we have all the components necessary to implement our sieve. Or do we?

### Powers of Primes as Factors

Our sieve, as described above, has one major flaw. It can identify which values of $x$ will result in a $y$ divisible by $p$, but it cannot identify 
whether or not this $y$ is divisible by a _power_ of $p$. In order to determine that, we would need to perform trial division on the value
to be sieved, until it is no longer divisible by $p$.

> **Edit:** This is incorrect. The roots of $n \mod p^k$ can be computed directly from the roots mod $p$ through the use of 
[Hensel Lifting](https://en.wikipedia.org/wiki/Hensel%27s_lemma#Hensel_lifting). This current implementation does precisely this.

We seemed to have reached an impasse: the whole point of the sieve was so that we _didn't_ have to do that. Time to check the playbook.

$\ln(a⋅b⋅c⋅d⋅…)=\ln(a)+\ln(b)+\ln(c)+\ln(d)+\ln(…)$


That looks pretty useful. If the sum of the $\ln$ of all of the small prime factors of $y$ is close to the expected value of $ln(y)$, then it's 
almost a given that $y$ has no other factors. In addition, if we adjust the expected value down a little bit, we can also identify values as 
smooth which have several powers of primes as factors. In this way, we can use the sieve as a 'pre-screening' process, and only factor 
those values which are likely to be smooth.

This has a few other advantages as well. Note that small primes contribute very little to the $\ln$ sum, but yet they require the most sieve time. 
Sieving the value 3 requires more time than 11, 13, 17, 19, and 23 _combined_. Instead, we can just skip the first few primes, 
and adjust the threshold down accordingly, assuming a certain percentage of them would have passed.

Another result, is that a number of values will be allowed to 'slip through', which are mostly smooth, but contain a single large cofactor. 
We could just discard these values, but suppose we found another mostly smooth value, with exactly the same cofactor.
We can then use these two values to construct a usable $y$; since their product will contain this large cofactor squared, 
it no longer needs to be considered.

### Putting it all together

The last thing we need to do is to use these values of $y$ construct an adequate $x$ and $d$. Suppose we only consider the non-square factors 
of $y$, that is, the prime factors of an odd power. Then, each $y$ can be expressed in the following manner:

$y0=p00⋅p11⋅p12⋯p0n$

$y1=p10⋅p01⋅p12⋯p1n$

$y2=p00⋅p01⋅p02⋯p1n$

$y3=p10⋅p11⋅p02⋯p0n$

etc., which can be expressed in the matrix form:

|0 & 1 & 1 & \ldots & 0 |

|1 & 0 & 0 & \ldots & 1 |

|0 & 0 & 0 & \ldots & 1 |

|1 & 1 & 0 & \ldots & 0 |

| \vdots &  &  &  &  |  ```

The problem then becomes to find a vector $\vec{v}$ such that $\vec{v}M\=\vec{0} \mod 2 $, where $\vec{0}$ is the null vector. That is, to solve 
for the left null space of $M$. This can be done in a number of ways, the simplest of which is to perform Gaussian Elimination on $M^T$, replacing the 
row addition operation with a row `xor`. This will result in a number of null space basis vectors, any combination of which will produce a valid solution.

The construction of $x$ is fairly straight-forward. It is simply the product of $Ax+B$ for each of the $y$ used. The construction of $d$ is slightly 
more complicated. If we were to take the product of all $y$, we will end up with a value with 10s of thousands, if not 100s of thousands of 
digits, for which we need to find the square root. This calcuation is impractically expensive. Instead, we can keep track of the even 
powers of primes during the sieving process, and then use `and` and `xor` operations on the vectors of non-square factors to reconstruct the square root.

I seem to have reached the 30000 character limit. Ahh well, I suppose that's good enough. Saved a bunch of bytes by switching to MathJax.
