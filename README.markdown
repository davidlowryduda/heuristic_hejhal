# Heuristic Hejhal #

This is an implementation of Hejhal's algorithm to compute data associated to
Maass forms on $\Gamma_0(N)$ for squarefree $N$.

The core of this began as a demonstration implementation in python and sagemath
and was written in 2020. To make it a little bit faster, computationally
intensive parts have been rewritten in cython.


## Brief overview of algorithm ##

The idea behind this algorithm is simple.
Begin with a guess $r_0$ of the spectral parameter for a Maass form on
$\mathrm{SL}(2, \mathbb{Z})$ (for simplicity, we begin with the full modular
group). A Maass form $f$ with the eigenvalue $\frac{1}{4} + r_0^2$ would have
an expansion of the form
```math
f(z) = \sum_{n \neq 0} a(n) \sqrt{y} K_{ir_0}(2 \pi |n| y) e^{2 \pi i n x},
```
where $K_{\nu}(y)$ is the modified $K$-Bessel function.
The $K$-Bessel function rapidly decays, so truncating introduces small error
```math
\widehat{f}(z) \approx \sum_{0 < |n| < M} a(n) \sqrt{y} K_{ir_0}(2 \pi |n| y) e^{2 \pi i n x}.
```

Compute these approximations along a horocycle $\{ z_m = x_m + iY \}$ for
$2Q$ equally spaced points $x_m$.
Evaluating $\widehat{f}$ at these points amounts to taking a discrete Fourier
transform.
Inverting the transform shows that
```math
a(n) \sqrt{Y} K_{ir_0}(2 \pi |n| Y)
=
\frac{1}{2Q}
\sum_{1-Q = m}^Q
\widehat{f}(z_m) e(-nx_m).
```

For fixed $r_0$ and $Y$, this is a linear system in the coefficients.
But the system is very poorly behaved.
To improve the system, we use that $f(\gamma z) = f(z)$ for all $\gamma$ in the
full modular group.

For each point $z_m$ in the horocycle, let $z_m^\dagger$ denote the corresponding
point in the fundamental domain (i.e. there should be some $\gamma_m$ such that
$\gamma z_m = z_m^\dagger$ and $z_m^\dagger$ is in the fundamental domain).
Then $f(z_m) = f(z_m^\dagger)$, the the action of $\gamma_m$ is generically highly
nonlinear.
By choosing the $Y$ in the horocycle small enough, no point $z_m$ will be in
the fundamental domain and every point will map nontrivially to another point
$z_m^\dagger$.

Set up the linear system based on
```math
    a(n) \sqrt{Y} K_{ir_0}(2 \pi |n| Y)
    =
    \frac{1}{2Q}
    \sum_{1-Q = m}^Q
    \widehat{f}(z_m^\dagger) e(-nx_m),
```
where we use $z_m^\dagger$ here instead. Expanding each of the expansions of
$\widehat{f}$ on the right, we now have a (homogenous, overdetermined,
approximate) linear system in the coefficients.

Solving for the coefficients gives a vector $\vec{a}_{r_0}$.

Of course, we don't know the actual eigenvalue.
So it won't be true that $f(z_m) = f(z_m^\dagger)$ with the claimed coefficients.
Instead, there will be some error.
We might study a **cost function** or **error function** in the coefficients,
say
```math
C(r) = a_{r}(2) \cdot a_r(3) - a_r(6).
```
If $r$ were a true eigenvalue, then $C(r) = 0$ as we would have a Hecke
eigenform.
Thus finding an $r$ such that $C(r) = 0$ is heuristically a good way to find a
Maass form.
(Better cost functions exist, or one could also check other Hecke relations,
etc).

To find a Maass form, we treat the entire linear system construction as a way
of evaluating the cost function $C(r)$.
Then we can treat this as a root-finding problem in $C(r)$.
For example, we might use the secand method and produce an iterative sequence
of guesses.
Given two guesses $r_0$ and $r_1$, we can produce the next estimate $r_2$ via
```math
r_2 = \frac{r_0 C(r_1) - r_1 C(r_0)}{C(r_1) - C(r_0)}.
```
If $r_0$ and $r_1$ are close enough to a true eigenvalue, we expect this to
converge extremely rapidly.


## Adjustments for higher level ##

When the level isn't $1$, it's necessary to simultaneously solve all the
expansions around every cusp at once.
When $N$ is squarefree, it turns out that every cusp will have either a sine or
a cosine expansion (as observed by Fredrik Stromberg) and it suffices to work
in real numbers.
Further, each expansion at each cusp is essentially the same, up to a single
sign called the Atkin-Lehner sign at that cusp.

If we don't know the signs, we check every possibility for each guess $r_0$.

See https://davidlowryduda.com/talk-computing-and-verifying-maass-forms/ or
other talks I've given for more details.


# Running this code #

This code is relatively straightforward to read but is a little fragile.
But typing `make cython` then `make python` then `make tests` is a good start.


# LICENSE #

Unless otherwise specified, the code here is available under GLPv3+. See the
included license file for more.
