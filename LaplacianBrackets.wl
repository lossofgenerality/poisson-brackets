(* ::Package:: *)

(* ::Title:: *)
(*Laplacian Brackets*)
(*  The bracket a Laplacian leaves behind when it fails the product rule:*)
(*  carre du champ \[CapitalGamma], its iterate \[CapitalGamma]2, its antisymmetric twin (Poisson), its odd twin (BV).*)
(*  Companion to Laplacian Lab (Computational Classnotes notebook, DOI 10.13140/RG.2.2.33742.01600).*)
(*  \[Copyright] 2026 Computational Classnotes \[CenterDot] CC BY-NC-SA 4.0*)

(* ::Text:: *)
(* Every section ends in named checks. Run the whole file; the last line prints*)
(* the table `checks`, every entry should be True. Works in a stateless kernel.*)

ClearAll["Global`*"];
checks = <||>;
check[name_String, test_] := (checks[name] = TrueQ[test]);
$Assumptions = 0 < \[Theta] < Pi && y > 0;


(* ::Section:: *)
(*1. Definition*)

(* ::Text:: *)
(* For any linear operator L on functions, the Laplacian bracket of L is*)
(*      \[CapitalGamma]_L(f,g) = 1/2 ( L(f g) - f L(g) - g L(f) ).*)
(* It measures exactly how far L is from being a derivation.  A first-order*)
(* operator (a vector field) gives \[CapitalGamma] = 0; a second-order operator gives*)
(* its symbol; for L = \[CapitalDelta] on R^n it gives \[CapitalGamma](f,g) = \[Del]f\[CenterDot]\[Del]g.*)
(* Iterating once more (Bakry\[Dash]\[CapitalEAcute]mery):*)
(*      \[CapitalGamma]2(f,g) = 1/2 ( L \[CapitalGamma](f,g) - \[CapitalGamma](f, L g) - \[CapitalGamma](L f, g) ).*)

grad[h_, v_List] := D[h, {v}];
lap[v_List] := Function[h, Total[D[h, {#, 2}] & /@ v]];

\[CapitalGamma]Op[L_][f_, g_] := Expand[(L[f g] - f L[g] - g L[f])/2];
\[CapitalGamma]2Op[L_][f_, g_] :=
  Expand[(L[\[CapitalGamma]Op[L][f, g]] - \[CapitalGamma]Op[L][f, L[g]] - \[CapitalGamma]Op[L][L[f], g])/2];

v3 = {x, y, z};
\[CapitalDelta]3 = lap[v3];
F = f @@ v3; G = g @@ v3; H = h @@ v3;

check["1a  \[CapitalGamma]_\[CapitalDelta](f,g) = \[Del]f\[CenterDot]\[Del]g in R^3",
  Simplify[\[CapitalGamma]Op[\[CapitalDelta]3][F, G] - grad[F, v3] . grad[G, v3]] === 0];

(* general second-order operator  L = a^{ij} \[PartialD]i\[PartialD]j + b^i \[PartialD]i  :  \[CapitalGamma] = a^{ij} \[PartialD]i f \[PartialD]j g, drift drops out *)
aM = Array[a, {3, 3}] /. a[i_, j_] /; i > j :> a[j, i];
aM = aM /. a[i_, j_] :> a[i, j] @@ v3;
bV = Array[b[#] @@ v3 &, 3];
Lgen = Function[u, Sum[aM[[i, j]] D[u, v3[[i]], v3[[j]]], {i, 3}, {j, 3}] + bV . grad[u, v3]];
check["1b  \[CapitalGamma]_L = symbol of L (drift b drops out)",
  Simplify[\[CapitalGamma]Op[Lgen][F, G] - grad[F, v3] . aM . grad[G, v3]] === 0];

(* first-order operator: no bracket at all *)
Vfield = Function[u, bV . grad[u, v3]];
check["1c  vector field has \[CapitalGamma] = 0", Simplify[\[CapitalGamma]Op[Vfield][F, G]] === 0];


(* ::Section:: *)
(*2. Intuition: \[CapitalGamma] is the covariance of f and g on a small sphere*)

(* ::Text:: *)
(* Laplacian Lab \[Section]1: mean of h over the sphere S(x,r) = h(x) + r^2/(2n) \[CapitalDelta]h + O(r^4).*)
(* Apply it to f g, f, g and subtract: the covariance of the fluctuations*)
(*      mean[(f - f(x)) (g - g(x))] = r^2/n \[CapitalGamma](f,g) + O(r^4).*)
(* So \[CapitalDelta] is "mean minus value" and \[CapitalGamma] is "covariance on the sphere": how much*)
(* f and g wiggle together around x.  \[CapitalGamma](f,f) = |\[Del]f|^2 is the variance.*)

sphereMean[u_, {x0_, y0_, z0_}, r_] := Module[{\[Theta]s, \[Phi]s},
  Integrate[(u /. {x -> x0 + r Sin[\[Theta]s] Cos[\[Phi]s], y -> y0 + r Sin[\[Theta]s] Sin[\[Phi]s],
        z -> z0 + r Cos[\[Theta]s]}) Sin[\[Theta]s], {\[Theta]s, 0, Pi}, {\[Phi]s, 0, 2 Pi}]/(4 Pi)];

fp = x^2 y + z^3 - x z; gp = x y z^2 + y^3 - 2 x;
pt = {p1, p2, p3}; at = Thread[v3 -> pt];
cov = Expand@sphereMean[(fp - (fp /. at)) (gp - (gp /. at)), pt, r];
check["2a  sphere covariance = r^2/3 \[CapitalGamma](f,g) + O(r^4)",
  Simplify[Coefficient[cov, r, 2] - (\[CapitalGamma]Op[\[CapitalDelta]3][fp, gp] /. at)/3] === 0 &&
   Coefficient[cov, r, 0] === 0 && Coefficient[cov, r, 1] === 0];

(* stochastic reading: for Brownian motion X_t, d<f(X),g(X)>_t = 2 \[CapitalGamma](f,g)(X_t) dt  (generator \[CapitalDelta]) *)


(* ::Section:: *)
(*3. Properties as an algebra*)

(* ::Text:: *)
(* (C, \[CapitalGamma]) for a second-order L with L(1) = 0:*)
(*  symmetric and bilinear;  a biderivation (Leibniz in each slot);*)
(*  chain rule \[CapitalGamma](\[CurlyPhi](f), g) = \[CurlyPhi]'(f) \[CapitalGamma](f,g);*)
(*  diffusion property L \[CurlyPhi](f) = \[CurlyPhi]'(f) L f + \[CurlyPhi]''(f) \[CapitalGamma](f,f);*)
(*  positivity \[CapitalGamma](f,f) \[GreaterEqual] 0 when the symbol is positive (elliptic);*)
(*  integration by parts \[Integral] \[CapitalGamma](f,g) d\[Mu] = -\[Integral] f L g d\[Mu]  (L symmetric w.r.t. \[Mu]).*)
(* It is NOT a Lie bracket (symmetric, no Jacobi).  Its order-detector is the*)
(* Koszul hierarchy: K_{k+1} \[Congruent] 0  \[DoubleLongLeftRightArrow]  L has order \[LessEqual] k.  For \[CapitalDelta]: K2 = 2\[CapitalGamma], K3 = 0.*)

\[CapitalGamma]3 = \[CapitalGamma]Op[\[CapitalDelta]3];
check["3a  symmetric", Simplify[\[CapitalGamma]3[F, G] - \[CapitalGamma]3[G, F]] === 0];
check["3b  bilinear", Simplify[\[CapitalGamma]3[\[Alpha] F + \[Beta] H, G] - \[Alpha] \[CapitalGamma]3[F, G] - \[Beta] \[CapitalGamma]3[H, G]] === 0];
check["3c  Leibniz (biderivation)", Simplify[\[CapitalGamma]3[F, G H] - G \[CapitalGamma]3[F, H] - H \[CapitalGamma]3[F, G]] === 0];
check["3d  chain rule", Simplify[\[CapitalGamma]3[\[CurlyPhi][F], G] - \[CurlyPhi]'[F] \[CapitalGamma]3[F, G]] === 0];
check["3e  diffusion property",
  Simplify[\[CapitalDelta]3[\[CurlyPhi][F]] - \[CurlyPhi]'[F] \[CapitalDelta]3[F] - \[CurlyPhi]''[F] \[CapitalGamma]3[F, F]] === 0];
check["3f  \[CapitalGamma](f,f) = |\[Del]f|^2 \[GreaterEqual] 0",
  Simplify[\[CapitalGamma]3[F, F] - Total[grad[F, v3]^2]] === 0];

(* integration by parts on the flat torus [0,2\[Pi]]^2 *)
v2 = {x, y}; \[CapitalDelta]2 = lap[v2];
ft = Sin[x] Cos[2 y] + Cos[3 x]; gt = Cos[x]^2 Sin[y] + Sin[x + y];
check["3g  \[Integral]\[CapitalGamma](f,g) = -\[Integral] f \[CapitalDelta]g (torus)",
  Integrate[\[CapitalGamma]Op[\[CapitalDelta]2][ft, gt] + ft \[CapitalDelta]2[gt], {x, 0, 2 Pi}, {y, 0, 2 Pi}] === 0];

(* not a Lie bracket: the Jacobi sum does not vanish *)
jac = \[CapitalGamma]3[x, \[CapitalGamma]3[x^2, y^2 x]] + \[CapitalGamma]3[x^2, \[CapitalGamma]3[y^2 x, x]] + \[CapitalGamma]3[y^2 x, \[CapitalGamma]3[x, x^2]];
check["3h  Jacobi fails (\[CapitalGamma] is not Lie)", Simplify[jac] =!= 0];

(* Koszul hierarchy *)
K3[L_][p_, q_, s_] := Expand[L[p q s] - L[p q] s - L[p s] q - L[q s] p + L[p] q s + L[q] p s + L[s] p q - L[1] p q s];
check["3i  K3 = 0 for \[CapitalDelta] (order 2)", Simplify[K3[\[CapitalDelta]3][F, G, H]] === 0];
\[CapitalDelta]sq = Function[u, \[CapitalDelta]3[\[CapitalDelta]3[u]]];
check["3j  K3 \[NotEqual] 0 for \[CapitalDelta]^2 (order 4): \[CapitalGamma] no longer a derivation",
  Simplify[K3[\[CapitalDelta]sq][x, x, x^2]] === 24 &&
   Simplify[\[CapitalGamma]Op[\[CapitalDelta]sq][x^2, y^2 x^2] - y^2 \[CapitalGamma]Op[\[CapitalDelta]sq][x^2, x^2] - x^2 \[CapitalGamma]Op[\[CapitalDelta]sq][x^2, y^2]] === 16 x^2];


(* ::Section:: *)
(*4. Antisymmetric twin: Poisson bracket, and the Hermitian bracket*)

(* ::Text:: *)
(* \[CapitalGamma](f,g) = \[Del]f\[CenterDot]g^{-1}\[CenterDot]\[Del]g uses the metric (symmetric);  {f,g} = \[Del]f\[CenterDot]J\[CenterDot]\[Del]g uses the*)
(* symplectic form (antisymmetric).  Both are biderivations; only the Poisson*)
(* one satisfies Jacobi.  On C = R^2 (a K\[ADoubleDot]hler manifold) they are the real and*)
(* imaginary parts of one Hermitian bracket:*)
(*      \[CapitalGamma](f,g) + i {f,g} = 4 \[PartialD]_z f \[CenterDot] \[PartialD]_zbar g.*)

Jmat[n_] := ArrayFlatten[{{0, IdentityMatrix[n]}, {-IdentityMatrix[n], 0}}];
v4 = {q1, q2, p1, p2};
PB[f_, g_] := Expand[grad[f, v4] . Jmat[2] . grad[g, v4]];
F4 = f @@ v4; G4 = g @@ v4; H4 = h @@ v4;
check["4a  Poisson antisymmetric", Simplify[PB[F4, G4] + PB[G4, F4]] === 0];
check["4b  Poisson Jacobi", Simplify[PB[F4, PB[G4, H4]] + PB[G4, PB[H4, F4]] + PB[H4, PB[F4, G4]]] === 0];
check["4c  Poisson Leibniz", Simplify[PB[F4, G4 H4] - G4 PB[F4, H4] - H4 PB[F4, G4]] === 0];
check["4d  canonical {q_i,p_j} = \[Delta]_ij", Table[PB[v4[[i]], v4[[j + 2]]], {i, 2}, {j, 2}] === IdentityMatrix[2]];

F2 = f[x, y]; G2 = g[x, y];
dz[u_] := (D[u, x] - I D[u, y])/2; dzb[u_] := (D[u, x] + I D[u, y])/2;
PB2[f_, g_] := D[f, x] D[g, y] - D[f, y] D[g, x];
check["4e  \[CapitalGamma] + i{ , } = 4 \[PartialD]z f \[PartialD]zbar g",
  Simplify[\[CapitalGamma]Op[\[CapitalDelta]2][F2, G2] + I PB2[F2, G2] - 4 dz[F2] dzb[G2]] === 0];
(* and \[CapitalDelta] = 4 \[PartialD]z \[PartialD]zbar, so \[CapitalDelta] itself is the "trace" of the Hermitian bracket *)
check["4f  \[CapitalDelta] = 4 \[PartialD]z\[PartialD]zbar", Simplify[\[CapitalDelta]2[F2] - 4 dz[dzb[F2]]] === 0];


(* ::Section:: *)
(*5. Manifolds: Laplace\[Dash]Beltrami, metric from \[CapitalGamma], Bochner, curvature*)

(* ::Text:: *)
(* On (M,g):  \[CapitalDelta]_g h = |g|^{-1/2} \[PartialD]i( |g|^{1/2} g^{ij} \[PartialD]j h ).*)
(* (1) \[CapitalGamma](x^i, x^j) = g^{ij}: the Laplacian alone determines the metric.*)
(* (2) \[CapitalGamma](r,r) = 1 characterises distance functions (eikonal equation).*)
(* (3) Bochner:  \[CapitalGamma]2(f,f) = |Hess f|^2 + Ric(\[Del]f,\[Del]f).  Curvature appears as the*)
(*     defect of the iterated bracket; Ric \[GreaterEqual] K and dim n give the*)
(*     curvature-dimension condition CD(K,n):  \[CapitalGamma]2(f) \[GreaterEqual] (\[CapitalDelta]f)^2/n + K \[CapitalGamma](f).*)
(* (4) CD(K,n) with K>0 \[DoubleLongRightArrow] Lichnerowicz \[Lambda]1 \[GreaterEqual] nK/(n-1); sharp on spheres.*)

LB[gm_, v_] := Function[u, With[{gi = Inverse[gm], s = Simplify[Sqrt[Det[gm]]], n = Length[v]},
    Simplify[Sum[D[s gi[[i, j]] D[u, v[[j]]], v[[i]]], {i, n}, {j, n}]/s]]];

christoffel[gm_, v_] := With[{gi = Inverse[gm], n = Length[v]},
   Table[Simplify@Sum[gi[[k, l]] (D[gm[[l, i]], v[[j]]] + D[gm[[l, j]], v[[i]]] - D[gm[[i, j]], v[[l]]])/2, {l, n}],
    {k, n}, {i, n}, {j, n}]];
ricci[gm_, v_] := With[{c = christoffel[gm, v], n = Length[v]},
   Table[Simplify@Sum[D[c[[a, d, b]], v[[a]]] - D[c[[a, a, b]], v[[d]]] +
       Sum[c[[a, a, e]] c[[e, d, b]] - c[[a, d, e]] c[[e, a, b]], {e, n}], {a, n}], {b, n}, {d, n}]];
hess[gm_, v_][u_] := With[{c = christoffel[gm, v], n = Length[v]},
   Table[D[u, v[[i]], v[[j]]] - Sum[c[[k, i, j]] D[u, v[[k]]], {k, n}], {i, n}, {j, n}]];

bochnerDefect[gm_, v_][u_] := With[{gi = Inverse[gm], Hs = hess[gm, v][u], du = grad[u, v]},
   Simplify[\[CapitalGamma]2Op[LB[gm, v]][u, u] - Tr[gi . Hs . gi . Hs] - du . gi . ricci[gm, v] . gi . du]];

(* round sphere S^2 *)
vS = {\[Theta], \[Phi]}; gS = {{1, 0}, {0, Sin[\[Theta]]^2}};
\[CapitalDelta]S = LB[gS, vS]; \[CapitalGamma]S = \[CapitalGamma]Op[\[CapitalDelta]S];
check["5a  S^2: \[CapitalGamma](x^i,x^j) = g^{ij}",
  Simplify[Outer[\[CapitalGamma]S, vS, vS] - Inverse[gS]] === {{0, 0}, {0, 0}}];
check["5b  S^2: Ric = g (K = 1)", Simplify[ricci[gS, vS] - gS] === {{0, 0}, {0, 0}}];
check["5c  S^2: Bochner \[CapitalGamma]2 = |Hess|^2 + Ric(\[Del]f,\[Del]f)", bochnerDefect[gS, vS][f[\[Theta], \[Phi]]] === 0];
check["5d  S^2: \[CapitalGamma](\[Theta],\[Theta]) = 1 (\[Theta] = distance to pole)", Simplify[\[CapitalGamma]S[\[Theta], \[Theta]]] === 1];

(* CD(1,2) on S^2: \[CapitalGamma]2 - (\[CapitalDelta]f)^2/2 - \[CapitalGamma](f) = |traceless Hess|^2 \[GreaterEqual] 0 *)
With[{u = f[\[Theta], \[Phi]], gi = Inverse[gS]},
  With[{T = hess[gS, vS][u] - \[CapitalDelta]S[u] gS/2},
   check["5e  S^2: CD(1,2) defect = |Hess\:2080 f|^2 \[GreaterEqual] 0",
    Simplify[\[CapitalGamma]2Op[\[CapitalDelta]S][u, u] - \[CapitalDelta]S[u]^2/2 - \[CapitalGamma]S[u, u] - Tr[gi . T . gi . T]] === 0]]];

(* Lichnerowicz: \[Lambda]1 = 2 = nK/(n-1); the eigenfunction cos \[Theta] has zero traceless Hessian (equality case) *)
check["5f  S^2: \[CapitalDelta] cos\[Theta] = -2 cos\[Theta]  (\[Lambda]1 = nK/(n-1) = 2)", Simplify[\[CapitalDelta]S[Cos[\[Theta]]] + 2 Cos[\[Theta]]] === 0];
check["5g  S^2: Hess cos\[Theta] = -cos\[Theta] g  (equality in Lichnerowicz)",
  Simplify[hess[gS, vS][Cos[\[Theta]]] + Cos[\[Theta]] gS] === {{0, 0}, {0, 0}}];

(* hyperbolic plane (upper half-plane) *)
vH = {x, y}; gH = {{1/y^2, 0}, {0, 1/y^2}};
check["5h  H^2: Ric = -g (K = -1)", Simplify[ricci[gH, vH] + gH] === {{0, 0}, {0, 0}}];
check["5i  H^2: Bochner", bochnerDefect[gH, vH][f[x, y]] === 0];
check["5j  H^2: \[CapitalGamma](log y, log y) = 1 (vertical geodesic distance)",
  Simplify[\[CapitalGamma]Op[LB[gH, vH]][Log[y], Log[y]]] === 1];

(* any conformal surface e^{2u}(dx^2+dy^2): \[CapitalDelta]_g = e^{-2u} \[CapitalDelta], Gaussian curvature K = -e^{-2u} \[CapitalDelta]u *)
gC = Exp[2 u[x, y]] IdentityMatrix[2];
check["5k  conformal: \[CapitalDelta]_g = e^{-2u}\[CapitalDelta]", Simplify[LB[gC, v2][F2] - Exp[-2 u[x, y]] \[CapitalDelta]2[F2]] === 0];
check["5l  conformal: Ric = K g, K = -e^{-2u}\[CapitalDelta]u",
  Simplify[ricci[gC, v2] + Exp[-2 u[x, y]] \[CapitalDelta]2[u[x, y]] gC] === {{0, 0}, {0, 0}}];
check["5m  conformal: Bochner for arbitrary u, f", bochnerDefect[gC, v2][F2] === 0];


(* ::Section:: *)
(*6. Graphs (Laplacian Lab \[Section]7): \[CapitalGamma] survives, the derivation property does not*)

(* ::Text:: *)
(* \[CapitalDelta]_G f(x) = \[Sum]_y w_xy (f(y) - f(x)).   Then*)
(*      \[CapitalGamma](f,g)(x) = 1/2 \[Sum]_y w_xy (f(y)-f(x)) (g(y)-g(x))   (discrete covariance),*)
(* summation by parts holds, but \[CapitalGamma] is no longer a biderivation and the chain rule*)
(* fails: this is why Bakry\[Dash]\[CapitalEAcute]mery curvature on graphs is subtle.*)

nG = 5; W = Table[If[i == j, 0, w[Min[i, j], Max[i, j]]], {i, nG}, {j, nG}];
\[CapitalDelta]G[fv_] := W . fv - Total[W, {2}] fv;
\[CapitalGamma]G[fv_, gv_] := (\[CapitalDelta]G[fv gv] - fv \[CapitalDelta]G[gv] - gv \[CapitalDelta]G[fv])/2;
fv = Array[\[Xi], nG]; gv = Array[\[Eta], nG]; hv = Array[\[Zeta], nG];
check["6a  graph \[CapitalGamma] = half weighted covariance",
  Simplify[\[CapitalGamma]G[fv, gv] - Table[Sum[W[[i, j]] (fv[[j]] - fv[[i]]) (gv[[j]] - gv[[i]]), {j, nG}]/2, {i, nG}]] === ConstantArray[0, nG]];
check["6b  graph summation by parts", Simplify[Total[\[CapitalGamma]G[fv, gv]] + fv . \[CapitalDelta]G[gv]] === 0];
check["6c  graph \[CapitalGamma] is not a derivation", Simplify[\[CapitalGamma]G[fv, gv hv] - gv \[CapitalGamma]G[fv, hv] - hv \[CapitalGamma]G[fv, gv]] =!= ConstantArray[0, nG]];


(* ::Section:: *)
(*7. Odd twin: the BV Laplacian and its antibracket*)

(* ::Text:: *)
(* Add odd (Grassmann) partners \[Xi]_i of the coordinates x_i and put*)
(*      \[CapitalDelta]_BV = \[Sum]_i \[PartialD]/\[PartialD]x_i \[PartialD]/\[PartialD]\[Xi]_i ,   \[CapitalDelta]_BV^2 = 0.*)
(* The same "failure of Leibniz" now produces a genuine Lie bracket of degree 1:*)
(*      (-1)^|a| {a,b} = \[CapitalDelta](ab) - \[CapitalDelta](a) b - (-1)^|a| a \[CapitalDelta](b).*)
(* Result: a Batalin\[Dash]Vilkovisky algebra (Gerstenhaber bracket generated by an*)
(* order-2 square-zero operator).  With \[Xi]_i read as \[PartialD]/\[PartialD]x_i, superfunctions are*)
(* polyvector fields and {,} is the Schouten bracket; the master equation*)
(* \[CapitalDelta]S + (1/2\[HBar]){S,S} = 0 is BV quantisation of gauge theories.*)
(* Representation: Association  sortedOddIndices -> even coefficient c(x).*)

nBV = 3; xs = Array[Symbol["x" <> ToString[#]] &, nBV];
sClean[a_Association] := Select[Expand /@ a, # =!= 0 &];
sPlus[a___Association] := sClean@Merge[{a}, Total];
sTimes[a_Association, b_Association] := sClean@Merge[Flatten@KeyValueMap[
      Function[{s1, c1}, KeyValueMap[Function[{s2, c2},
         If[IntersectingQ[s1, s2], Nothing, Sort[Join[s1, s2]] -> Signature[Join[s1, s2]] c1 c2]], b]], a], Total];
sScale[c_, a_Association] := sClean[c a];
dx[i_][a_Association] := sClean[D[#, xs[[i]]] & /@ a];
d\[Xi][i_][a_Association] := sClean@Merge[KeyValueMap[Function[{s, c},
      If[MemberQ[s, i], DeleteCases[s, i] -> (-1)^(First@FirstPosition[s, i] - 1) c, Nothing]], a], Total];
\[CapitalDelta]BV[a_Association] := sPlus @@ Table[dx[i][d\[Xi][i][a]], {i, nBV}];
par[a_Association] := If[a === <||>, 0, Mod[Length[First@Keys[a]], 2]];
anti[a_Association, b_Association] := With[{e = (-1)^par[a]},
   sScale[e, sPlus[\[CapitalDelta]BV[sTimes[a, b]], sScale[-1, sTimes[\[CapitalDelta]BV[a], b]], sScale[-e, sTimes[a, \[CapitalDelta]BV[b]]]]]];
sZeroQ[a_Association] := sClean[a] === <||>;

X[i_] := <|{} -> xs[[i]]|>; \[CapitalXi][i_] := <|{i} -> 1|>;
check["7a  {x_i, \[Xi]_j} = \[Delta]_ij  (canonical antibracket)",
  Table[anti[X[i], \[CapitalXi][j]], {i, nBV}, {j, nBV}] === Table[If[i == j, <|{} -> 1|>, <||>], {i, nBV}, {j, nBV}]];

SeedRandom[2026];
rpoly[] := Sum[RandomInteger[{-3, 3}] Times @@ (xs^RandomInteger[{0, 2}, nBV]), {4}];
rsuper[p_] := sClean@Association@Table[s -> rpoly[], {s, Select[Subsets[Range[nBV]], Mod[Length[#], 2] == p &]}];
trials = Flatten[Table[{rsuper[p1], rsuper[p2], rsuper[p3]}, {p1, 0, 1}, {p2, 0, 1}, {p3, 0, 1}], 2];
sg[a_, b_] := (-1)^((par[a] + 1) (par[b] + 1));

check["7b  \[CapitalDelta]_BV^2 = 0", AllTrue[trials, sZeroQ[\[CapitalDelta]BV[\[CapitalDelta]BV[#[[1]]]]] &]];
check["7c  graded antisymmetry {a,b} = -(-1)^((|a|+1)(|b|+1)) {b,a}",
  AllTrue[trials, sZeroQ@sPlus[anti[#[[1]], #[[2]]], sScale[sg[#[[1]], #[[2]]], anti[#[[2]], #[[1]]]]] &]];
check["7d  graded Jacobi",
  AllTrue[trials, Function[t, With[{a = t[[1]], b = t[[2]], c = t[[3]]},
     sZeroQ@sPlus[anti[a, anti[b, c]], sScale[-1, anti[anti[a, b], c]], sScale[-sg[a, b], anti[b, anti[a, c]]]]]]]];
check["7e  Poisson (Leibniz) rule {a,bc} = {a,b}c + (-1)^((|a|+1)|b|) b{a,c}",
  AllTrue[trials, Function[t, With[{a = t[[1]], b = t[[2]], c = t[[3]]},
     sZeroQ@sPlus[anti[a, sTimes[b, c]], sScale[-1, sTimes[anti[a, b], c]],
       sScale[-(-1)^((par[a] + 1) par[b]), sTimes[b, anti[a, c]]]]]]]];
check["7f  \[CapitalDelta] is a derivation of {,}: \[CapitalDelta]{a,b} = {\[CapitalDelta]a,b} + (-1)^(|a|+1){a,\[CapitalDelta]b}",
  AllTrue[trials, Function[t, With[{a = t[[1]], b = t[[2]]},
     sZeroQ@sPlus[\[CapitalDelta]BV[anti[a, b]], sScale[-1, anti[\[CapitalDelta]BV[a], b]],
       sScale[-(-1)^(par[a] + 1), anti[a, \[CapitalDelta]BV[b]]]]]]]];


(* ::Section:: *)
(*Report*)

Dataset[checks]
