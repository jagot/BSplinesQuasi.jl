# Theory

The underlying equations for the definition and efficient evaluation
of B-splines are introduced. Equations labelled (dB.A.ij) refer to
equation (ij) in chapter A of

- Carl de Boor (2001). _A Practical Guide to Splines_. New York:
  Springer. ISBN: 978-0-387-95366-3.

## Definitions

The _order_ of the polynomial is designated $k$, equal to the degree + 1, i.e. a
parabolic polynomial is of order 3.

The B-splines can be defined through the Cox–de Boor recursion relation:
$$\begin{equation}
\tag{dB.IX.14}
\B{j}{k} \defd
\omega_{jk}\B{j}{k-1} +
(1-\omega_{j+1,k})\B{j+1}{k-1},
\end{equation}$$
where
$$\begin{equation}
\tag{dB.IX.15}
\omega_{jk}(x) \defd
\frac{x-t_j}{t_{j+k-1}-t_j},
\end{equation}$$
and
$$\begin{equation}
\tag{dB.IX.11}
\B{r}{1}(x) = \begin{cases}
1, & x\in [t_r,t_{r+1}),\\
0, & \textrm{else},\\
\end{cases}\quad
r\in[1,n_t-1].
\end{equation}$$

Given a knot vector $\vec{t}$ of length $n_t$, there are $n_{tk}\defd
n_t-k$ functions of order $k$. This implies that there is a highest
order $k$ a given knot set vector can support, i.e. $k_{\textrm{max}}
= n_t - 1$.

## Knot sets

As seen above, the B-splines are completely defined by the knot set
vector $\vec{t}$.

The knot set $\vec{t}=\bmat{1&2&3&4&5&6}$ gives rise to the following
cardinal splines, of orders $k=1..5$:
![Cardinal splines](figures/cardinal-splines.svg)

By increasing the _multiplicity_ of some knots, the continuity of the
splines can be controlled. E.g. the knot set
$\vec{t}=\bmat{0&1&1&3&4&6&6&6}$, will yield the following splines:
![Discontinuous splines](figures/discontinuous-splines.svg)

Lastly, it is very common to pad the knot set such that the first and
last knot have multiplicity $k$; this simplifies the implementation of
boundary conditions when B-splines are used for solving differential
equations:
![Full multiplicity splines](figures/full-multiplicity-splines.svg)

## de Boor's algorithm

An efficient way of evaluating the B-splines is given by [de Boor's
algorithm](https://en.wikipedia.org/wiki/De_Boor%27s_algorithm). The
algorithm described in the Wikipedia article [^1] assumes full
multiplicity at the endpoints of the knot set, i.e. that the first and
last points are repeated $k$ times. In the description of the same
algorithm, de Boor additionally assumes a strictly increasing knot
set, i.e. $t_{i+1}>t_i$, do avoid any divisions by zero. Since one of
the goals for this package is maximum flexibility in choosing the knot
set, a derivation of an only slightly more general version of de
Boor's algorithm follows.

Given a knot set $\vec{t}$ of length $n_t$, the spline $f(x)$
is given by
$$\begin{equation}
\tag{dB.X.23⅓ \& ⅔}
\begin{aligned}
f(x) &= \sum_{r=1}^{n_t-k} \alpha_r \B{r}{k}(x)\\
&=
\sum_{r=1}^{n_t-k} \alpha_r\frac{x-t_r}{t_{r+k-1}-t_r}\B{r}{k-1}(x)+
\sum_{r=1}^{n_t-k} \alpha_r\frac{t_{r+k}-x}{t_{r+k}-t_{r+1}}\B{r+1}{k-1}(x)\\
&=
\sum_{r=1}^{n_t-k} \alpha_r\frac{x-t_r}{t_{r+k-1}-t_r}\B{r}{k-1}(x)+
\sum_{r=2}^{n_t-k+1} \alpha_{r-1}\frac{t_{r+k-1}-x}{t_{r+k-1}-t_{r}}\B{r}{k-1}(x)\\
&=
\frac{x-t_1}{t_k-t_1}\alpha_1\B{1}{k-1}(x) +
\left[
\sum_{r=2}^{n_t-k}
\frac{(x-t_r)\alpha_r+(t_{r+k-1}-x)\alpha_{r-1}}{t_{r+k-1}-t_r}\B{r}{k-1}(x)
\right] +
\frac{t_{n_t}-x}{t_{n_t}-t_{n_t-k+1}}\alpha_{n_t-k}\B{n_t-k+1}{k-1}(x)\\
&=
\sum_{i=1}^{n_t-k+1}
\alpha_i^{[2]}(x)
\B{i}{k-1}(x),
\end{aligned}
\end{equation}$$
where
$$\begin{equation}
\tag{dB.X.24}
\alpha_r^{[2]}(x) \defd
\begin{cases}
\displaystyle
\frac{x-t_1}{t_k-t_1}\alpha_1, & r=1,\\[2ex]
\displaystyle
\frac{(x-t_r)\alpha_r+(t_{r+k-1}-x)\alpha_{r-1}}{t_{r+k-1}-t_r}, & r\in[2,n_t-k],\\[2ex]
\displaystyle
\frac{t_{n_t}-x}{t_{n_t}-t_{n_t-k+1}}\alpha_{n_t-k}, & r = n_t - k + 1.
\end{cases}
\end{equation}$$
We have thus reexpressed the spline function $f(x)$ of order $k$ as a
linear combination of B-splines of order $k-1$. We can generalize
this, to reexpress $f(x)$ as a linear combination of B-splines of
order $k-j$, with expansion coefficients
$$\begin{equation}
\tag{dB.X.26}
\alpha_r^{[j+1]}(x) \defd
\begin{cases}
\displaystyle
\frac{x-t_1}{t_{1+k-j}-t_1}\alpha_1^{[j]}(x), & r=1,\\[2ex]
\displaystyle
\frac{(x-t_r)\alpha_r^{[j]}(x)+(t_{r+k-j}-x)\alpha_{r-1}^{[j]}(x)}{t_{r+k-j}-t_r}, & r\in[2,n_t-k+j-1],\\[2ex]
\displaystyle
\frac{t_{n_t}-x}{t_{n_t}-t_{n_t-k+j}}\alpha_{n_t-k+j-1}^{[j]}(x), & r = n_t - k + j.
\end{cases}
\end{equation}$$

The difference between this derivation and those in the Wikipedia
article and de Boor (2001), is that we here explicitly consider the
limits of the sum imposed by the length of the knot set and the order
$k$; this introduces the special cases for $r=1,n_t-k+j$.

An important philosophical difference between the Cox–de Boor
recursion relation and de Boor's algorithm, is that whereas the former
is a linear combination of basis functions evaluated at certain
position $x$, the latter is linear combination of _intervals_ (since
the first-order functions $\B{j}{1}$ are non-zero within one interval
only, and they are mutual orthogonal), with _polynomial_ expansion
coefficients $\alpha_{i}^{[k]}(x)$. To evaluate the spline function
$f(x)$, we first find the interval $i$ which contains $x$. Even if the
knot set is only non-decreasing, i.e. not strictly increasing, the
interval containing $x$ is uniquely defined, since there is only one
for which $t_i \leq x < t_{i+1}$; if the knot $t_i$ has a multiplicity
higher than unity, the additional intervals cannot contain $x$, since
they are empty: $t_{i-1}\leq x < t_i = \varnothing$ if $t_{i-1} =
t_i$. By finding the _last_ $i$, for which $t_i \leq x$, we thus
guarantee that no divisions by zero will occur.

[^1]: NB that the Wikipedia article uses 0-based indexing, whereas de Boor and BSplinesQuasi.jl use 1-based indexing.
