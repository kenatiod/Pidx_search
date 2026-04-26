# Products of Consecutive Integers in Primorial Intervals

## 1. Background

This note presents the heuristic and computational justification for expecting the 
OEIS sequence [A141399]( https://oeis.org/A141399), the prime-complete products of the form

$$ N_2(n)=n(n+1), $$

to be finite, with termination at $n = 633{,}555.$ The sequence comprises:

$$n = 1, 2, 3, 5, 8, 9, 14, 15, 20, 24, 35, 80, 125, 224, 384, 440, 539, 714, 1715, 
2079, 2400, 3024, 4374, 9800, 12375, 123200, 194480, 633555 $$

Thus $N_2$ ranges from $1 \times 2 = 2 \text{ to } 633{,}555 \times 633{,}556 = 401{,}392{,}571{,}580$ with
each value of $N_2$ having a complete set of prime divisors from 2 to its greatest prime divisor. The last
value of $N_2$ has the divisor set \{ $2, 3, 5, 7, 11, 13, 17, 19$ \} 
giving it a distinct-prime-divisor count $\omega(N_2(n))$
of 8, which matches the prime number index $\pi(gpf(N_2(n)))$ of its greatest prime divisor, 19.
In the text below, we will use $Pidx(N_2(n))$ to refer to the prime number index $\pi(gpf(N_2(n)))$.
For any integer greater than 1, $\omega \le Pidx$, and are equal for prime-complete numbers.

The products of consecutive integers have been extensively studied in the field of Number Theory [3]. As these
products increase, two properties become prominent: increasing number of prime divisors and increasing prime
number index of the greatest prime factor of the pair $(n, n+1)$, our $Pidx(N_2(n))$. The $\omega(N_2)$ value is
limited on the high side by the $N_2$ position in the progression of primorial numbers, and $Pidx(N_2)$ is
limited on the low side by what mathematicians call the "smoothness" of the integer pairs.

At the OEIS A141399 post, Giovanni Resta states, "a(29) > 2.29*10^25, if it exists." If true, this break in finding
successive terms is very large, and puts consecutive integer pairs well past the point where any heuristic
argument from probability of smoothness would expect to find pairs where Pidx would be low enough to match
$\omega$. However, the central point is not that the functions $ω(N_2(n))$ and $Pidx(N_2(n))$ are monotone. 
They are not. Rather, the point is that primorial-scale intervals give a rigorous upper 
barrier for $ω$, while computation and smooth-number heuristics suggest that the 
interval minimum of $Pidx$ eventually rises past that barrier.

This is a heuristic and computational note, not the formal proof. The formal non-heuristic proof 
will come from a finite Diophantine reduction, such as Størmer/Lehmer/Pell enumeration 
for fixed $ω$, together with a tail-closing argument.


## 2. Definitions

Let $N_2(n)=n(n+1).$

Let $ω(m)$ denote the number of distinct prime divisors of $m$.

Let $gpf(m)$ denote the greatest prime factor of $m$, and define 

$$Pidx(m)=π(gpf(m)),$$ 

where $π(x)$ is the prime-counting function.

A positive integer $m$ is called **prime-complete** if its set of distinct 
prime divisors is exactly \{ $2,3,5,…,p_k$ \} for some $k$. Equivalently,

$$ m \text{ is prime-complete} ⟺ ω(m)=Pidx(m).$$

For $N_2(n)=n(n+1)$, this means that all primes up to the greatest prime factor 
occur in the product, and no missing primes appear in its distinct prime divisor set.

Since $n$ and $n+1$ are coprime, $ω(n(n+1))=ω(n)+ω(n+1)$, and $Pidx(n(n+1))=max{(Pidx(n),Pidx(n+1))}.$


## 3. The primorial interval barrier

Let $P_r=2⋅3⋅5⋯p_r$ be the $r$-th primorial number. The basic deterministic observation is:

**Lemma.** If $m<P_{r+1}$, then $ω(m)≤r.$

**Proof.** If $ω(m)≥r+1$, then $m$ has at least $r+1$ distinct prime factors. 
The least possible product of $r+1$ distinct primes is $2⋅3⋅5⋯p_{r+1}=P_{r+1}.$

Therefore $m≥P_{r+1}$, is a contradiction. ◻

Now choose intervals of $n$ such that $n(n+1)$ lies between successive primorial thresholds. 
Informally, the $r-th$ interval is the range of $n$ for which $P_r≤n(n+1)<P_{r+1}.$ 
On such an interval, the lemma gives the rigorous upper bound $ω(N_2(n))≤r.$ 
Therefore, in that interval, a prime-complete value must satisfy

$$Pidx(N_2(n))=ω(N_2(n))≤r.$$ 

This gives the key exclusion test:

$$
\boxed{\min_{n\in I_r}Pidx(N_2)>r\quad\Longrightarrow\quad
\text{no prime-complete } n(n+1) \text{ values occur in } I_r.}
$$

The important feature of this test is that it does not require an exact computation of 
the interval maximum of $ω(N_2(n)).$ The value $r$ is already a rigorous upper 
bound for ω throughout the interval.


## 4. What the computation measures

The program primorial_maxO_minP.c searches intervals determined by primorial-root 
boundaries and records, among other diagnostics, $minPidx(N_2(n))$ within each interval.

The column maxOmega is useful diagnostic information, but it is not the main mathematical certificate. 
If factorization is deliberately limited, then the measured maxOmega may be incomplete. The primorial 
argument avoids this difficulty: the decisive comparison is not $minPidx>maxω,$ but rather $minPidx>r.$ 
That comparison uses the exact primorial barrier $ω(N_2(n))≤r$ for the interval.

The computation is reliable for detecting any possible counterexample with $Pidx≤r,$ provided the trial 
factorization table includes at least the first $r$ primes. If a prime-complete value occurred in the 
interval, then its greatest prime factor would be at most $p_r$, so both $n$ and $n+1$ would factor completely 
over the first $r$ primes and would be detected by such a search.



## 5. Why $Pidx$ is expected to outrun $r$

The primorial barrier gives the structural ceiling: $ω(N_2)≤r$ when $N_2<P_{r+1}.$ 

Using the prime number theorem in the Chebyshev-function form,

$$
\log(P_r)=\vartheta(p_r)\sim p_r,
$$

and the usual estimate

$$
p_r\sim r\log r,
$$

we obtain, on the scale $n(n+1)\asymp P_r$,

$$
2\log n \asymp \log(P_r)\sim p_r\sim r\log r.
$$

Thus

$$
r\asymp \frac{2\log n}{\log\log n}.
$$

This is the largest possible scale for $ω(N_2)$ in that primorial interval.

By contrast, the heuristic least attainable greatest prime factor in a large 
interval is governed by smooth-number rarity. 
A value m near size $x$ is $y-smooth$ with probability heuristically controlled 
by the Dickman-de Bruijn function $\rho(u) \text{, with }u=\frac{\log(⁡x)}{\log⁡(y)}$.

The smaller $y$ is relative to $x$, the rarer such values become. For products $n(n+1)$ near size $n^2,$ 
the computation suggests that the minPidx in primorial-scale intervals behaves roughly on a scale like

$$ 
\pi((\log(n))^2) \sim \frac{(\log(n))^2}{2\log(\log(n))}.
$$

This is much larger than the primorial barrier scale $r≍2log⁡(n)log⁡(log(⁡n)).$ Their ratio is heuristically 

$$
\frac{
(\log(n))^2/(2\log(\log(n)))
}{
2\log(n)/\log(\log(n))
}
\sim
\frac{\log(n)}{4},
$$

which tends to infinity. Thus the expected behavior is eventual permanent separation: 

$$ \min_{n\in I_r}Pidx(N_2)>r.$$

Once that inequality holds for all later primorial intervals, prime-complete products $n(n+1)$ can no longer occur.


## 6. Why this is still heuristic

The primorial barrier

$$ \omega(N_2)\le r $$

inside the interval is rigorous.

The computational exclusion of a specific interval is also rigorous if the program has actually 
verified that no $n$ in the interval satisfies

$$ Pidx(N_2)\le r. $$

However, the assertion that this exclusion persists for all sufficiently large $r$ is not 
proved merely by the smooth-number heuristic. Dickman-type estimates describe the 
density or probability of smooth numbers in large ranges; by themselves they do not exclude 
all exceptional values in every later interval. Therefore this note should be read as strong 
evidence for finiteness, not as the formal proof of finiteness.

A formal proof must supply one of the following:

- a rigorous lower bound proving that $\min_{n\in I_r}Pidx(N_2(n))>r$ for all sufficiently large $r$; or

- an independent structural argument eliminating all $\text{large-}\omega$ cases, such
  as a Størmer/Lehmer/Pell finite reduction for fixed $\omega$,
  together with a CRT/LCM obstruction closing the tail.

## 7. Relation to Størmer’s theorem

For a fixed finite set of primes $P$, Størmer’s theorem implies that there are only 
finitely many pairs of consecutive P-smooth integers. Since a prime-complete value with

$$ \omega(N_2)=r $$

requires both $n$ and $n+1$ to be smooth over

$$ \{2,3,5,\ldots,p_r\}, $$

Størmer’s theorem supplies the natural finite-reduction framework for the formal proof. In computational 
terms, fixed-r prime-complete candidates may be searched through Pell-equation families. 

The heuristic in this note explains why those families should eventually become empty; 
the formal proof must show that they do.

## 8. Summary

The primorial Pidx barrier is:

$$
P_r \le n(n+1) \lt P_{r+1} \quad\Longrightarrow\quad \omega(n(n+1))\le r. 
$$

Therefore a prime-complete value in that interval must have

$$ Pidx(n(n+1))\le r. $$

So the decisive computational test is:

$$ \min_{n\in I_r}Pidx(n(n+1))>r. $$

This test avoids reliance on a measured maximum of $ω.$ The value r itself 
supplies the rigorous upper bound.

The reason finiteness is expected is that r grows only on the scale

$$ \frac{2\log n}{\log(\log(n))}, $$

while smooth-number heuristics and computation suggest that the least possible greatest-prime-factor index in these 
intervals grows on a much larger scale, roughly

$$ \frac{(\log(n))^2}{2\log(\log(n))}. $$

That expected separation is the heuristic reason A141399 should terminate. 
The formal proof must then rule out exceptional later values.

## 9. Computational Corroboration 
Substantial computational corroboration has gone into this argument. (See [data here](https://github.com/kenatiod/Delta_min))
The smallest differences (minDelta) between
values of $Pidx(N_2)$ and $\omega(N_2)$ have been examined for values of $n$ up to $2 \times 10^{14}$ thus $N_2$ 
up to $4 \times 10^{28}$. The last locations of $n$ for small $Pidx(N_2) - \omega(N_2)$ are:

| Pidx - $\omega$ | Last Found at $n$ |
| --------------- | ------------------------: |
| 0 | 633,555 |
| 1 | 80,061,344 |
| 2 | 1,109,496,723,125 |
| 3 |1,284,729,638,049 |
| 4 |20,628,591,204,480|
| 5 |11,259,346,386,959|
| 6 |119,089,041,053,696|

Observed maxOmega to minPidx growth over doubling intervals:
| $log_2$ Start | maxOmega | minPidx |
| ------------- | -------- | ------- |
| 1 | 2 | 2 |
| 2 | 3 | 2 |
| 3 | 4 | 2 |
| 4 | 4 | 3 |
| 5 | 4 | 4 |
| 6 | 5 | 3 |
| 7 | 6 | 4 |
| 8 | 6 | 5 |
| 9 | 7 | 5 |
| 10 | 7 | 6 |
| 11 | 7 | 4 |
| 12 | 8 | 4 |
| 13 | 8 | 5 |
| 14 | 8 | 7 |
| 15 | 9 | 7 |
| 16 | 9 | 6 |
| 17 | 10 | 7 |
| 18 | 10 | 7 |
| 19 | 11 | 8 |
| 20 | 11 | 9 |
| 21 | 11 | 9 |
| 22 | 11 | 8 |
| 23 | 12 | 8 |
| 24 | 12 | 10 |
| 25 | 13 | 12 |
| 26 | 13 | 10 |
| 27 | 13 | 10 |
| 28 | 13 | 11 |
| 29 | 14 | 13 |
| 30 | 14 | 11 |
| 31 | 14 | 12 |
| 32 | 15 | 16 |
| 33 | 15 | 14 |
| 34 | 15 | 14 |
| 35 | 16 | 13 |
| 36 | 16 | 16 |
| 37 | 16 | 15 |
| 38 | 16 | 14 |
| 39 | 17 | 16 |
| 40 | 17 | 15 |
| 41 | 18 | 18 |
| 42 | 18 | 18 |
| 43 | 18 | 17 |
| 44 | 18 | 17 |
| 45 | 19 | 20 |
| 46 | 19 | 20 |

If we let $x = \log_2(n)$, then theory-guided fit gives maxOmega growing roughly like $\frac{x}{\log(x)}$, while minPidx is better fit
by a function with an $\frac{x^2}{\log(x)}$ component. In the observed data the fitted curves cross near
$x=44$ ($n = 2^{44}$ or about 18 trillion), matching the point where minPidx begins to exceed maxOmega. Here here are the growth curve fit
equations:

$$
maxOmega(x)  
\approx  
-0.30+1.45\frac{x}{\ln(x+1)}+0.50\ln(x+1),  
$$

$$
minPidx(x)  
\approx  
0.66+0.83\frac{x}{\ln(x+1)}+0.0159\frac{x^2}{\ln(x+1)}.  
$$

Replacing $x$ by $\log_2(n)$, the observed maximum omega grows like a logarithmic-over-logarithmic
function of $n$, while the observed minimum Pidx has a quadratic-logarithmic term. This is the
numerical signature of the expected separation between available prime-factor count and greatest-prime-factor index.

## References

1. K. Dickman, “On the frequency of numbers containing prime factors of a certain relative magnitude”,
    *Arkiv för Matematik, Astronomi och Fysik*, 22A(10), 1930.

2. N. G. de Bruijn, “On the number of positive integers ≤x and free of prime factors >y”, 
   *Indagationes Mathematicae*, 13, 1951, 50–60.

3. G. Tenenbaum, *Introduction to Analytic and Probabilistic Number Theory*, Graduate Studies in 
   Mathematics, Vol. 163, American Mathematical Society, 2015. See especially the 
   treatment of smooth numbers and the Dickman-de Bruijn function.

4. G. H. Hardy and E. M. Wright, *An Introduction to the Theory of Numbers*, 
   Oxford University Press. See the prime number theorem, the Chebyshev functions, 
      and standard asymptotics for primes and primorials.

5. C. Størmer, “Quelques théorèmes sur l’équation de Pell x2−Dy2=±1 
   et leurs applications”, *Skrifter Videnskabsselskabet i Kristiania*, 1897.
