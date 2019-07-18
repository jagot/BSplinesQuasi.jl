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
algorithm described in the Wikipedia article assumes full multiplicity 