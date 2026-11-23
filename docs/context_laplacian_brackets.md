# CONTEXT — Laplacian Brackets

© 2026 Computational Classnotes · CC BY-NC-SA 4.0
Last updated: 23 November 2026

File: `LaplacianBrackets.wl` (Wolfram Language, runs in a fresh kernel, ends with a table of 42 named checks; all True on 23 Nov 2026).

## Working definition
Laplacian bracket of a linear operator L: Γ_L(f,g) = ½[L(fg) − fLg − gLf] (carré du champ). Iterate: Γ₂(f,g) = ½[LΓ(f,g) − Γ(f,Lg) − Γ(Lf,g)]. For L = Δ on Rⁿ, Γ = ∇f·∇g.

## Sections and checks
1. Definition: Γ_Δ = ∇f·∇g; for a general 2nd-order L, Γ is the symbol (drift drops out); vector fields give 0.
2. Intuition: covariance on a sphere, mean[(f−f(x))(g−g(x))] = r²/n Γ(f,g) + O(r⁴) (from Laplacian Lab §1).
3. Algebra: symmetric, bilinear, biderivation, chain rule, diffusion property, positivity, integration by parts; not Lie (Jacobi fails); Koszul K₃ = 0 for Δ, ≠ 0 for Δ².
4. Poisson twin: {f,g} = ∇f·J∇g (antisymmetric, Jacobi, Leibniz, canonical); Γ + i{ , } = 4∂_z f ∂_z̄ g; Δ = 4∂_z∂_z̄.
5. Manifolds: Laplace–Beltrami; Γ(xⁱ,xʲ) = gⁱʲ (Δ determines the metric); eikonal Γ(r,r)=1; Bochner Γ₂ = |Hess|² + Ric on S², H², arbitrary conformal surface; CD(1,2) on S²; Lichnerowicz equality on S² (cos θ).
6. Graphs: Γ = half weighted covariance; summation by parts; not a derivation.
7. BV: superfunctions as Association(odd indices → coefficient); Δ_BV = Σ∂x∂ξ; Δ²=0; antibracket graded antisymmetric, graded Jacobi, Poisson rule, Δ a derivation of { , } (random-test over all 8 parity patterns).

## Lessons
- The Δ² non-derivation check needs degree-2 inputs (x, y, xy are too low degree: K₃ = 0 there).
- Sign conventions for BV verified: (−1)^|a|{a,b} = Δ(ab) − Δ(a)b − (−1)^|a| aΔ(b); signs (−1)^((|a|+1)(|b|+1)).

## Possible next steps
Koszul–Brylinski bracket on forms of a Poisson manifold; Bakry–Émery curvature on graphs (cycle, complete graph); weighted Laplacian Δ − ∇V·∇ and Bakry–Émery Ricci; adding a Laplacian-brackets section to Laplacian Lab.
