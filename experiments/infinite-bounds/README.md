# Landau-continued Cauchy transforms by sinc/cot product quadrature

A numerical-math note, self-contained: the problem, the method, why it is
optimal, its limits, and a Base-only reference implementation (`cauchy.jl`).
Application notes for this repository live in `application-vmd.md`.

## 1. The problem

Given `g: ℝ → ℂ` smooth, decaying, analytic in a strip `|Im u| < d`, evaluate
its **Cauchy transform** with the **Landau prescription**:

    C⁺[g](ζ) = ∫_ℝ g(u)/(u − ζ) du   continued analytically from Im ζ > 0
             = integral + 2πi·g(ζ)   for Im ζ < 0
             = PV + iπ·g(ζ)          for Im ζ = 0.

The Hilbert transform is the on-axis PV case. The plasma dispersion function
Z is the g = Gaussian instance.

Structure the method exploits:
- **g is the expensive factor**; the kernel is known in closed form.
- **Ladder evaluation**: applications often need C⁺[g](ζ_n) at many points for
  the SAME g — amortizing samples across points matters as much as
  single-point cost.
- ζ can be anywhere: upper half-plane, lower (continuation active), exactly
  real, or far away (|g(ζ)| huge or underflowed).

## 2. The method

Map `u = ψ(t) = u_c + S·sinh(t)` (u_c, S ≈ center and rms width of |g|), take
the plain trapezoid on uniform nodes `t_j = j·h`, and correct the pole in
closed form:

    C⁺[g](ζ) = h·Σ_j G(t_j)/(ψ(t_j) − ζ)  +  π·g(ζ)·( cot(π·t*/h) + i ),
    G = g(ψ)·ψ′,   t* = asinh((ζ − u_c)/S).

Derivation: residue calculus on `(π/h)·cot(πw/h)·G(w)/(ψ(w)−ζ)` over the
analyticity strip. Key facts:
- The residue at t* is **exactly g(ζ)** — the map's Jacobian cancels.
- The `+iπ` constant merges Im ζ ≷ 0 AND the Landau 2πi into ONE analytic
  formula: no crossed bookkeeping, no branch logic, exactly-real ζ included.
- Error ~ `e^(−2πd_t/h)`, d_t = strip width of G in t. h=0.2 gives ~1e-9 for
  thermal-like g; halve h for ~1e-13.
- **Nodes must be exact multiples of h.** An offset grid (e.g. `(-T):h:T` with
  T/h ∉ ℤ) shifts the cot phase → O(1) errors; non-monotone h-convergence is
  the symptom.
- Saturate `cot+i` past |Im| ≈ 20 (`_cot_i` in cauchy.jl): raw cot overflows
  to NaN near |Im w| ≈ 700; the exact limits are 0 (above) / 2i (below).
- `g` must not overflow at large real |u| (scaled special functions, not raw).
- The sinh map compresses tails logarithmically: algebraic decay only
  lengthens the node window by a log factor, never changes the convergence
  rate (the rate is set by the strip, i.e. by g's complex singularities).

Validation of `cauchy.jl` against the exact oracle `C⁺[e^(−u²)] = √π·Z(ζ)`,
`Z(ζ) = i√π·erfcx(−iζ)`: ≤7e-11 at h=0.2 across damped / growing / exactly
real / strongly damped / far-pole / marginal (Im ζ = 1e-12); algebraic-tail g
`(1+u²/3)^(−2.6)`: h vs h/2 agree to 2e-15 on an 11-point ladder.

## 3. Why the trapezoid is fundamental (not incidental)

The integrand factorizes: [g: strip-analytic, decaying, expensive] × [Cauchy
kernel: singular but known]. Rank any method by how much of that it uses:

1. **Product quadrature**: interpolate only g, integrate the interpolant
   against the exact kernel → nodes resolve g alone; node count independent of
   pole count/placement; one sample set serves the whole ladder.
2. **Equispacing ⇔ closed-form pole sum**: `h·Σ_j 1/(jh−t*)` IS the
   partial-fraction expansion of `−π·cot(πt*/h)` — the quadrature's own error
   at a pole is the correction, in one term. Any other node family loses this
   (Chebyshev → Q_k recurrences per pole; Gauss → secondary functions).
   Equivalent view: trapezoid = exact integration of the sinc interpolant,
   whose Cauchy transform is closed-form.
3. **Trapezoid optimality**: for strip-analytic decaying functions on ℝ the
   equispaced trapezoid converges at `e^(−2πd/h)` and is asymptotically
   unbeatable (Trefethen–Weideman). Compactifying to (−1,1) discards this and
   turns algebraic tails into endpoint branch points (algebraic convergence);
   under sinh the same tails only shorten the node window.
4. **Meromorphy in ζ = free Landau continuation**: cot is one analytic
   function across the axis. Log-based corrections have their branch cut ON
   the integration path → explicit 2πi bookkeeping (a recurring bug source,
   §6).
5. **No engineered cancellation**: subtraction methods rely on `g(u)−g(ζ)`
   cancelling and go ill-conditioned for large |g(ζ)|; here the correction is
   O(g(ζ))·(cot+i) with (cot+i) → 0 exponentially as the pole leaves the
   strip.

General principle: match the interpolant class to g's smoothness, integrate it
exactly against the Cauchy kernel. Gauss-weighted polynomials → Z functions;
piecewise polynomials → per-cell logs; strip-analytic → sinc/cot. When g is
only piecewise-smooth (splines, tabulated data) the trapezoid drops to
algebraic order — use the piecewise closed form instead.

## 4. Known limitations

- **No a-posteriori error control**: h is fixed a priori; multi-scale g (e.g.
  narrow core + wide halo, or a secondary population far from u_c) degrades
  silently. Cheap fix: h vs h/2 check — nested grids share every node, so
  validation costs 2× once, not per call.
- **Single global scale**: u_c/S come from first/second moments; features
  |u − u_c| ≫ S sit where node spacing `S·h·cosh(t)` is coarse.
- **Node collision**: Re t* within ~1e-3·h of a node with Im t* → 0 loses
  relative accuracy to pair cancellation (measure-zero; nudge or half-shift
  the grid if it ever matters).
- Truncation at `u_c ± S·sinh(T)` (T=7 → ~550·S).

## 5. Reference implementation

`cauchy.jl`, Base-only, composable:

    cauchy_landau(g, ζ)                        # scalar, default alg
    cauchy_landau(g, ζs, alg)                  # ladder: g-samples shared across all ζs
    alg = SincCot(map = SinhMap(uc, S), h = 0.2, window = 7.0)

The map is a composable axis because the cot correction is MAP-INDEPENDENT:
for any analytic bijection ψ: ℝ → ℝ the residue of G(t)/(ψ(t)−ζ) at
t* = ψ⁻¹(ζ) is exactly g(ζ). A map (`mapto`/`mapjac`/`mapinv`) only chooses
where nodes land and the strip geometry that sets the rate:

- `SinhMap(uc, S)` (default): logarithmic tail compression — algebraic decay
  costs only a log-wider window.
- `LinearMap(uc, S)`: plain trapezoid — optimal for entire g with fast decay
  (Gaussian test: 7e-15 at h=0.25, window=6), but the window must cover the
  support in units of S, so heavy tails want SinhMap.

Validated: SinhMap vs Z oracle ≤7e-11 (§2); drifted g via `uc` against the
shifted oracle 1e-12; kappa ladder h vs h/2 to 2e-15.

## 6. Superseded approaches (each fixed the previous obstacle)

The naive route — adaptive quadrature with constant Plemelj subtraction
`∫(g−g(ζ))/(u−ζ) + g(ζ)·log((U−ζ)/(L−ζ))` — needs FINITE endpoints: the
subtracted constant leaves a log-divergent tail on ℝ and the log term is
undefined at ∞. Two intermediate repairs, both correct, both beaten:

1. **Decaying-kernel subtraction**: subtract `g(ζ)·k(u)` with
   k(u) = W²/(W²+(u−Re ζ)²); then `C(ζ) = ∫k/(u−ζ)` closes:
   iπW/(W+Im ζ) above the axis, 2πi − iπW/(W−Im ζ) below. Converges on ℝ
   (≤6e-8) but the kernel's own curvature costs ~1.8× in nodes and W is a
   free knob. A Gaussian kernel (erfcx closed form) uses identical node
   counts and is slower (an exp per pole per node). Trap: on a finite
   interval the Landau term is the FULL 2πi·g(ζ) once — scaling it by k(ζ)
   leaves a (Im ζ/W)²·g(ζ) error.
2. **Compactify + constant subtraction**: map u = S·t/(1−t²) to (−1,1); the
   residue in t is exactly g(ζ), so constant subtraction works again on the
   finite t-interval. Machine-ε accurate, ~2× the adaptive-box cost.
   Structural flaws the sinc/cot method removes: the map's second preimage
   (ghost pole at −1/t*) and endpoint branch points for fractional-power
   tails.

Adaptive quadrature retains one real advantage over all fixed-node schemes:
it measures its own error and localizes into multi-scale structure (§4).
