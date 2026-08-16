# RH Research Notebook V2 — Atlas Integration Checkpoint

## Accepted research trajectory

1. Scalar f0 verifier: prime-power cell splitting, strict convexity above log(plastic constant), downward derivative jumps at breakpoints.
2. f1 audit: midpoint-odd basis q1=x-L/2; exact overlap formula; active prime shifts strengthen the odd channel under G=G0-Gp+Ginf.
3. Even degree-2 block: bubble basis b=x(L-x); exact prime kernels; midpoint parity factorization.
4. Structural bounds: endpoint-jet filtration gives omitted-tail scaling O(log(T)/T^(2r+1)) when r endpoint jets vanish.
5. Fourier cross-check: T=84 selected; stable entire low-frequency transforms and support-length Taylor jets developed.

## Latest imported notebook state

The notebook reported a hybrid uniform T=84 degree-2 certificate and later independently reproduced direct-Fourier point certificates at log(3), an interior bottleneck near L=1.10595, and log(4). The latest direct-Fourier work reduced the missing uniform independent check to interval sign coverage of E2,84'' and E2,84'.

**Atlas policy:** none of those numerical claims is promoted by this integration commit. They are research inputs until regenerated in-repo.

## Direct-Fourier next target

Certify with outward-rounded interval arithmetic:

- E2,84''(L) > 0 on [log(3), 1.20]
- E2,84'(L) > 0 on [1.20, log(4)]
- one interval point ball near L=1.1059498113

Then infer uniform E2,84(L)>0.

## Future degree-3 basis

q1 = x-L/2
b3 = x(L-x)(x-L/2)

Do not start degree 3 until the independent degree-2 Fourier certificate is reproducibly closed.
