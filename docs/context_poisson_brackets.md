# CONTEXT — Poisson Brackets

© 2026 Computational Classnotes · CC BY-NC-SA 4.0
Last updated: 23 November 2026

File: `PoissonBrackets.wl` (Wolfram Language, fresh kernel, ends with a table of 50 named checks; all True on 23 Nov 2026). Notebook: `PoissonBrackets.nb` loads it and runs the functions with plots. Companion: `LaplacianBrackets.wl` (see context_laplacian_brackets.md).

## Sections
1. Definition: canonical {f,g}; general Poisson tensor Π^{ij}; Hamilton's equations ḟ = {f,H}.
2. Intuition: Jacobian/area; rate of change along a flow (exact oscillator flow, area preserved); Moyal star product: [q,p]⋆ = iħ, Moyal bracket → Poisson as ħ→0, exact for quadratic H, ⋆ associative.
3. Algebra: bilinear, antisymmetric, Jacobi, Leibniz, canonical relations, chain rule; [X_f,X_g] = −X_{f,g} with X_f h = {h,f}; Jacobiator test (Π = {{0,y,0},{−y,0,x},{0,−x,0}} fails Jacobi); every 2D bivector is Poisson.
4. Lie–Poisson on so(3)*: Jacobi, {m_i,m_j} = ε m_k, Casimir |m|², rigid body Euler equations (minus bracket), conservation of H and |m|².
5. Uses: Kepler (L, Runge–Lenz A conserved; so(4) relations {A_i,A_j} = −2H ε L_k); L3 generates rotations; leapfrog symplectic vs explicit Euler (det = 1 + h² cos q); KdV via Gardner bracket, with conservation laws checked by the Euler operator.
6. Geometry/manifolds: Liouville (div X_H = 0); Archimedes–Darboux {φ, m3} = 1; {θ,φ} = 1/(R sin θ); KKS form = area/R, total 4πR; action–angle {θ,I} = 1 with θ = ArcTan[p,q]; anisotropic oscillator in involution; geodesic flow from H = ½gⁱʲpᵢpⱼ on S², H², an arbitrary conformal surface; Clairaut; link to Laplacian brackets: Γ + i{ , } = 4∂_z f ∂_z̄ g.

## Lessons
- Never name a variable π in Wolfram code: it is Pi, and `π_` patterns break. Use Π (\[CapitalPi]).
- VariationalMethods` shadowing in a reused session: use a hand-written Euler operator (Σ(−∂s)^k ∂/∂u^(k)) instead.
- Both 3D "test" bivectors {{0,1,0},{−1,0,x},{0,−x,0}} and {{0,z,0},{−z,0,x},{0,−x,0}} turned out to be Poisson; the failing example is {{0,y,0},{−y,0,x},{0,−x,0}}.

## Possible next steps
Dirac brackets for constraints; Poisson cohomology / modular vector field; Koszul–Brylinski bracket and the Poisson Laplacian (bridge to LaplacianBrackets.wl); Kontsevich star product; Lax pairs and r-matrices for integrable systems.
