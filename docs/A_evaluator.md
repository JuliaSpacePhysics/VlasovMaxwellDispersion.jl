## Residue extraction

`T(a,z)` is meromorphic in `p∥`; its `1/sin πa` poles are the resonances `a=n`, i.e.
`p∥=ζ_n`, linear at fixed outer coordinate (`ζ_n(γ)=(ωγ−nΩ₀)/k∥`). Peel the in-range poles,

    ∫ I dp∥ = ∫ [ I − Σ_n ρ_n/(p∥−ζ_n) ] dp∥  +  Σ_n ρ_n·[ log((b−ζ_n)/(a−ζ_n)) (+2πi) ],

with residue `ρ_n = 𝓣_n /(∂a/∂p∥)`. The bracket is smooth (peaks removed ⇒ coarse 2-D
cubature); the pole term carries the Landau `+2πi·ρ_n` for `Im ζ_n<0`, so A handles damped modes. 
The pole count `~ k∥·(support width)/Ω₀` is **independent of `k⊥`**, so A stays flat in `k⊥ρ`; the `Σ_{m∉range} 𝓣_m/(a−m)` tail is what the closed form sums for free.

## Regularizing kernel at `z→0` (`k⊥→0`)

`T(a,z)` carries explicit `1/z, 1/z²`, and `Jₐ(z)=(z/2)^a·[even entire series]` diverges at `z=0` for Re a < 0. Yet `T(a,z)` is entire in z — the singularity is removable. 

`T(a,z)` needs only the **regularized quartet** below, each entire in z:

    σ0 = π J_{−a}Jₐ/sin πa = Σ_n Jₙ²/(a−n)
    σ1 = (a σ0 − 1)/z² = (1/z²) Σ_n n Jₙ²/(a−n)
    σD = (z/2) π P'/sin πa / z² = (1/z²) Σ_n z Jₙ Jₙ'/(a−n)
    σJ = π J'_{−a}J'ₐ/sin πa + a/z² = Σ_n (Jₙ')²/(a−n)

In closed form it is the hypergeometric `σ0 = (1/a)·₁F₂(½; 1+a, 1−a; −z²)`.

The kernel 𝓣 is then a polynomial in `(z,u,w)·σ` (`w=p⊥, u=p∥`):

    𝓣(a,z) = [ a σ1 w²    i a σD w²    σ1 z w u ;
              −i a σD w²   σJ w²      −i σD z w u ;
               σ1 z w u    i σD z w u   u² σ0     ].

It is manifestly finite at `z=0`, where it equals the `n=0,±1` harmonic sum (`Jₙ(0)=δₙ₀` kills `|n|≥2`); the only surviving poles are `a=0,±1`.

**Evaluation.** `σ` come from one `|z|`-split (the small-/large-argument switch every Bessel
routine already has):

- `|z| ≥ 1`: the closed Bessel form (`P, P'` from `Jₐ, J'ₐ`).
- `|z| < 1`: the entire z²-series (`₁F₂`). With `qₖ ≡ (π/sin πa)·pₖ`,
  `pₖ=(−1)^k (2k)!/(k!²\,Γ(a+k+1)Γ(1−a+k))`, the recurrence
  `qₖ = qₖ₋₁·[−(2k)(2k−1)/(k²(k²−a²))]`, `q₀=1/a`;
  then `σ0=Σ_k qₖ xᵏ`, `σ1=(a/4)Σ_{k≥1} qₖ x^{k−1}`, `σD=(1/4)Σ_{k≥1} k qₖ x^{k−1}`,
  `σJ = ½Σ_{k≥1} k² qₖ x^{k−1} + σ0 − a σ1`, with `x=(z/2)²`.

The series is ~10× faster per call at small `z` (no complex-order Bessel) and more accurate near `z=0`.
The stop uses rigorous tail bound (Johansson Thm 1): gate `k>|a|` (past the `k≈|a|` spike of `1/(k²−a²)`), then bound the slowest (k²-weighted σJ) tail by a negligible term.

A maintained uniform `pFq` (e.g. Slevinsky's sequence-transformation method) evaluates it without the `|z|`-split, but is several× slower per call and needs the three companion `₁F₂`'s for `σ1,σD,σJ`.

See [qin_sigmas_compare](../benchmark/qin_sigmas_compare.jl) for comparison of different implementations.

**Pitfalls / drawbacks.**

- *Series spike near integer `a`.* `qₖ` carries `1/(k²−a²)`, spiking at `k≈|a|`. The stop
  gates on `k>|a|` (never stops before the spike) and bounds the slow k²-weighted tail (Thm 1
  above); the 100-term cap covers `|a|` well beyond any physical case. Near-integer `a` is
  still conditioning-limited at f64 (`~1e-13`) — shared by the state of the art (Arb's only
  recourse there is more precision).
  — the principled fix is ε-power-series arithmetic (Johansson §2.3, §7.2, used by mpmath/Arb). Arb evaluates the removable singularity g(a)/sin(πa) at integer a by computing in ℂ[[ε]]/⟨ε²⟩ (auto-diff in the parameter), formally cancelling the zero. That cures exactly-integer a.
- *Poles are physical, not removed.* `σ` still diverge at `a=0,±1` (and the closed form at
  every integer `a`) — these are the cyclotron resonances, peeled by the residue machinery
  (§6). The regularization removes only the **spurious** `z=0` singularity, never a resonance.

Literature verdict — no uniform cancellation-free formula exists. Johansson, Arb, and mpmath all use the same approach: series for small |z|, asymptotic/acceleration for large |z| as the cancellation at large |z| is fundamental to ₁F₂.

Sources:
- F. Johansson, [Computing Hypergeometric Functions Rigorously](https://dl.acm.org/doi/10.1145/3328732) (series + asymptotic; parameter integer/fractional separation)
  - > "it would be a mistake to use (8) for all z … it would lead to exponentially large cancellation when z→+∞."
- [HypergeometricFunctions.jl](https://github.com/JuliaMath/HypergeometricFunctions.jl) 
  - R. M. Slevinsky, [Fast and stable rational approximation of generalized hypergeometric functions](https://doi.org/10.1007/s11075-024-01808-w)
- [mpmath Hypergeometric functions](https://mpmath.org/doc/current/functions/hypergeometric.html)