# Audited Formula Sheet

Let `d=L-a`, `q1=x-L/2`, and `b=x(L-x)`.

## General overlap

C_ij(a;L) = integral_0^(L-a) x^i (x+a)^j dx
            = sum_{r=0}^j binom(j,r) a^(j-r) d^(i+r+1)/(i+r+1).

K_ij = C_ij + C_ji.

## Degree 1

K_00 = 2d
K_01 = L d
K_11 = d^2(2L+a)/3
K_q1q1 = d(L^2 - 2La - 2a^2)/6

K_q1q1 < 0 when a/L > (sqrt(3)-1)/2.

## Bubble even block

K_0b = d^2(L+2a)/3
K_bb = d^3(L^2+3La+a^2)/15

det K_{1,b} = d^4(L^2-2La-14a^2)/45.

## Scalar curvature

Inside a prime-power cell, with r=e^L,

W00''(L) = 2(r^3-r-1)/(sqrt(r)(r^2-1)).

At L=log(q), q=p^k,

W00'(right)-W00'(left) = -2 Lambda(q)/sqrt(q) < 0.

## Fourier transforms

H0(t;L) = integral_0^L exp(itx) dx.

With z=Lt/2 and A(z)=sin(z)/z,
H0 = L exp(iz) A(z).

For b=x(L-x), B(z)=(sin z-z cos z)/z^3,
Hb = L^3/2 exp(iz) B(z).

Hb = -2i/t Hq1.

## Degree-2 parity factorization

Let O1,T=G_T[q1,q1], E2,T=det G_T[{1,b}], D2,T=det G_T[{1,x^2}]. Then

D2,T = E2,T + L^2 G00,T O1,T,
det G_T[{1,x,x^2}] = O1,T E2,T.

These identities are finite-cutoff identities as well as cutoff-free identities when all terms are defined consistently.
