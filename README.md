# Poisson Brackets

Wolfram Language code for the Poisson bracket and its symmetric twin, the Laplacian bracket (carré du champ): definitions, intuition, algebraic properties, applications in mechanics and field theory, and geometry on manifolds. Every claim is a named check that evaluates to `True`.

© 2026 Computational Classnotes · CC BY-NC-SA 4.0

## Files

| File | What it is |
|---|---|
| `PoissonBrackets.wl` | Poisson brackets: definition, intuition (area, flows, Moyal limit), Poisson algebra, Lie–Poisson on so(3)* and the rigid body, Kepler and so(4), symplectic integrators, KdV, symplectic geometry and geodesic flow. 50 checks. |
| `PoissonBrackets.nb` | Notebook that loads `PoissonBrackets.wl` and runs its functions, with plots (pendulum phase portrait, rigid-body orbits, Kepler orbit, leapfrog vs Euler, geodesic on S²). |
| `LaplacianBrackets.wl` | Laplacian brackets Γ and Γ₂: sphere-covariance intuition, algebra and Koszul hierarchy, Hermitian link Γ + i{ , }, Bochner formula and curvature-dimension on manifolds, graph Laplacians, BV antibracket. 42 checks. |
| `docs/` | Context notes for each file: what is covered, conventions, lessons, next steps. |
| `tools/make_nb.py` | Script that generates `PoissonBrackets.nb`. |

## Running

Keep `PoissonBrackets.nb` and `PoissonBrackets.wl` in the same folder, open the notebook and choose **Evaluation ▸ Evaluate Notebook**.

Or, in any Wolfram kernel:

```wolfram
Get["PoissonBrackets.wl"]      (* shows a Dataset of 50 checks *)
Get["LaplacianBrackets.wl"]    (* shows a Dataset of 42 checks *)
```

Each file starts with `ClearAll["Global`*"]`, so load one at a time in a fresh kernel.

## License

Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International — https://creativecommons.org/licenses/by-nc-sa/4.0/
