# Derivation

The susceptibility of a hot magnetized plasma is **one helical-orbit phase-space integral**.
This document shows how this integral collapses, step by step, under specific assumptions.

## Conventions

Normalize $v → v/c$, $ω → ω/Ω_{ref}$, $k → k c/Ω_{ref}$.

Let $Ω=qB/(γmc)$ be the relativistic gyrofrequency, $Ω₀=γΩ$ the signed constant rest
gyrofrequency, and $Π ≡ ω_p/Ω_{ref}$.

## 1. Linear Vlasov → orbit integral

Linearize relativistic Vlasov–Maxwell about a gyrotropic `f₀(p⊥,p∥)`.
With `θ=0` the orbit phase is

    β = −z[sin(φ + Ω τ) − sin φ] + (ω − k∥ v∥)τ.

Set `s=Ωτ`, and define the three dimensionless orbit parameters

    a ≡ (ω − k∥ v∥)/Ω = (γω − k∥p∥)/Ω₀,   z ≡ k⊥ v⊥/Ω = k⊥ p⊥/Ω₀,

with `φ` the azimuthal angle of `p`. Every field component then reduces to the single
orbit integral

    g(φ,z) = Ω⁻¹ ∫₀^∞ exp[−i z sin(φ+s) + i a s] ds.

## 2. Gyrophase average → susceptibility tensor

Assembling the current moment from `g(φ,z)` by a gyrophase average gives

    χ = (Π²)/(ω Ω₀) ∫ d³p · S

with

$$
S = -i\Omega
\begin{pmatrix}
p_⊥ U G_{33} & p_⊥ U G_{32} &
  p_⊥ ∂_{p_∥}f_0 G_{31} - p_⊥ V G_{33} \\
p_⊥ U G_{23} & p_⊥ U G_{22} &
  p_⊥ ∂_{p_∥}f_0 G_{21} - p_⊥ V G_{23} \\
p_∥ U G_{13} & p_∥ U G_{12} &
  p_∥ ∂_{p_∥}f_0 G_{11} - p_∥ V G_{13}
\end{pmatrix},
$$

the operators

$$
U = ∂_⊥f₀ + \frac{k_∥}{ωγ}(p_⊥ ∂_∥f₀ − p_∥ ∂_⊥f₀),
\qquad
V = \frac{k_⊥}{ωγ}(p_⊥ ∂_∥f₀ − p_∥ ∂_⊥f₀)
$$

and the gyrophase matrix `G[g]`

$$
G_{i j} \equiv \frac{1}{2 \pi} \int_0^{2 \pi} d φ\, e^{i z \sin φ }\left(\begin{array}{ccc}
g & i \frac{∂ g}{∂ z} & \frac{i}{z} \frac{∂ g}{∂ φ} \\
g \sin φ & i \frac{∂ g}{∂ z} \sin φ & \frac{i}{z} \frac{∂ g}{∂ φ} \sin φ \\
g \cos φ & i \frac{∂ g}{∂ z} \cos φ & \frac{i}{z} \frac{∂ g}{∂ φ} \cos φ
\end{array}\right).
$$

`S` has identical form in physical or normalized variables; the per-species `m_s, q_s`
sit only in the `Ω₀, Π` prefactor.

The `V` dependence (last column of `S`) factors into a compact closed form, the `e∥e∥`-nonresonant term.

Two evaluators diverge: **A** uses the periodicity of `g`; **B** expands it in harmonics.

## 3. Closing the orbit integral

### (A) Finite-period reduction (Qin's complex-order form)

Gyrophase symmetry `g(φ,z)=g(φ+2π,z)` splits the semi-infinite `s` integral into gyroperiods,

    g(φ,z) = c₀/(2π) ∫₀^{2π} exp[-i z sin(φ+η) + i a η] dη,
    c₀ = −e^{-iπa} π / (i Ω sin πa),

retaining every resonance in the `1/sin(πa)` factor (zeros at `ωγ−k∥p∥=nΩ₀`). The gyrophase matrix, for example `G₁₁`, can be reduced to

$$
G_{11} = -\frac{\pi J_{−a}(z)J_a(z)}{iΩ \sin πa},
$$

And S = e∥e∥ (Ω/ω)(p∥/p⊥)(p⊥ ∂∥f₀ − p∥ ∂⊥f₀) + p⊥ U T where the resonant 3×3 `T(a,z)` closes in `Q, Q', J'_{−a}J'_a, a, z` and `r=p∥/p⊥`:

    T = [
      a/z² (Q−1)              i/(2z) Q'                 (Q−1) r/z
     −i/(2z) Q'               π J'_{−a}J'_a/sin πa + a/z²   −i r Q'/(2a)
      (Q−1) r/z               i r Q'/(2a)                (Q/a) r²
    ]

where $σ_0 ≡ π J_{−a}(z)J_a(z)/\sin π a,   Q ≡ a σ_0$.

Note: Qin Eq. 35's printed `T₁₂,T₂₁,T₂₃,T₃₂` have typos.

The susceptibility splits into a resonant and a non-resonant (Bernstein) term:

$$
\chi = \chi_{T} + \chi_{B}, \qquad
\chi_{T} ≡ \frac{\Pi^2}{\omega \Omega_0} \int d^3p\; p_\perp U T(a,z),\quad
\chi_B ≡ \frac{\Pi^2}{\omega \Omega_0} \int d^3p\;\frac{\Omega}{\omega}\left(\frac{1}{p_\|}\frac{\partial f_0}{\partial p_\|} - \frac{1}{p_\perp}\frac{\partial f_0}{\partial p_\perp}\right)p_\|^2 .
$$

Price: `a` is complex ⇒ complex-order Bessel; and relativistically `sin πa=0` traces
curves in the momentum plane, leaving a 2-D integral plus residues (§3.1) — versus the
clean 1-D reduction B reaches in the separable case.

See [A_evaluator.md](A_evaluator.md) for more details.

### (B) Harmonic sum

Insert $e^{iz \sinφ}=Σ_m Jₘ(z)e^{imφ}$; the orbit integral becomes
harmonic sums like `Σₙ Jₙ(z)²/(a−n)`, and

$$
S = \mathbf{e}_{\|} \mathbf{e}_{\|}\, \frac{Ω}{ω}\left(\frac{1}{p_{\|}} \frac{∂ f_0}{∂ p_{\|}}-\frac{1}{p_{\perp}} \frac{∂ f_0}{∂ p_{\perp}}\right) p_{\|}^2  +\sum_{n=-\infty}^{\infty} \frac{Ω\, p_{\perp} U}{ω-k_{\|} v_{\|}-n Ω}\, \mathbf{T}_n,
$$

with the per-harmonic tensor

$$
\mathbf{T}_n \equiv\begin{pmatrix}
\frac{n^2 J_n^2}{z^2} & \frac{i n J_n J_n^{\prime}}{z} & \frac{n J_n^2 p_{\|}}{z p_{\perp}} \\
-\frac{i n J_n J_n^{\prime}}{z} & \left(J_n^{\prime}\right)^2 & -\frac{i J_n J_n^{\prime} p_{\|}}{p_{\perp}} \\
\frac{n J_n^2 p_{\|}}{z p_{\perp}} & \frac{i J_n J_n^{\prime} p_{\|}}{p_{\perp}} & \frac{J_n^2 p_{\|}^2}{p_{\perp}^2}
\end{pmatrix}.
$$

Using the ring kernel $R_n \equiv (n/z)J_n=(J_{n-1}+J_{n+1})/2$, $𝓣_n \equiv p_\perp^2 T_n$ reads

$$
𝓣_n=
\begin{pmatrix}
p_\perp^2 R_n^2 & i\,p_\perp^2 R_n J_n' & p_\parallel p_\perp R_n J_n\\
-i\,p_\perp^2 R_n J_n' & p_\perp^2 J_n'^2 & -i\,p_\parallel p_\perp J_nJ_n'\\
p_\parallel p_\perp R_n J_n & i\,p_\parallel p_\perp J_nJ_n' & p_\parallel^2 J_n^2
\end{pmatrix}.
$$

Finite for $z\to 0$ and $p_\perp\to0$.
Note: For _real_ $(p_\|,p_\perp)$, $𝓣_n$ is the Hermitian
outer product $\mathbf{v}_n\mathbf{v}_n^\dagger$ of $\mathbf{v}_n=(p_\perp R_n,-i p_\perp J_n',p_\|J_n)$.

The susceptibility splits into harmonic contributions plus one non-resonant term:

$$
\chi = \frac{\Pi^2}{\omega^2}\sum_n X_n + \mathbf{e}_\| \mathbf{e}_\|\, \chi_B,\qquad
X_n ≡ \int d^3p\;\frac{\omega p_\perp U /  \gamma}{\omega - k_\parallel v_\parallel - n\Omega}\;\mathbf{T}_n
 = \int d^3p\;\frac{𝒰 /  \gamma^2}{\omega - k_\parallel v_\parallel - n\Omega} 𝓣_n.
$$

where 𝒰 ≡ (ωγ/p⊥) U = k∥ ∂∥f₀ + (ωγ − k∥p∥)/p⊥ · ∂⊥f₀ = ω ∂f₀/∂γ + k∥ ∂f₀/∂p∥.

Insert $1=\sum_n J_n^2(z)$ into the $\chi_B$ integrand, with per-harmonic non-resonant term defined as

$$
X_n^{B} \equiv \int d^3p\;\frac{J_n^2(z)}{\gamma}\,p_\|^2\left(\frac{1}{p_\|}\frac{\partial f_0}{\partial p_\|}-\frac{1}{p_\perp}\frac{\partial f_0}{\partial p_\perp}\right),
$$

$\chi_B$ can also be expressed as a sum: $\chi_B = (\Pi^2/\omega^2) \sum_n X_n^{B}$. Defining $\chi_n \equiv \frac{\Pi^2}{\omega^2} (X_n + \mathbf{e}_\| \mathbf{e}_\|\, X_n^{B})$, then

$$
\chi=\sum_n\chi_n.
$$

Keep `|n| ≤ nmax ≈ k⊥ρ`; the sum converges slowly for large `z`. The payoff: each raw `X_n` **factorizes when `f₀` is separable** (§5).

(`σ0 = Σ_n Jₙ(z)²/(a−n) = π J_{−a}J_a/sin πa` is the Lerche–Newberger Bessel identity linking A and B.)

### The momentum integral

Either closure leaves the same shape,

    χ = (Π²/ω Ω₀) ∫ d³p [ e∥e∥-Bernstein + (p⊥ U) · K ],

with kernel $K$ = closed $T(a,z)$ (A) or harmonic sum $Σ_n T_n/(a−n)$ (B).

**Maxwellian / Maxwell–Jüttner / cold.** Gaussian `f` closes both primitives in _closed form_: `𝒞→Z` (plasma dispersion function), `P⊥→Γ_n=Iₙ(λ)e^{−λ}`. No quadrature
— fastest. (Still evaluator B: the `n`-sum is truncated; only per-`n` moments are closed.)

- See [Maxwellian.md](Maxwellian.md) for a worked example.

## 5. The outer-coordinate density `I`

The cyclotron resonance is a pole at $p_\parallel=\zeta_n$, at fixed outer
coordinate. Integrating inner $p_\parallel$ first — continued past the pole (Plemelj boundary value + growing-sheet Landau residue, §5.3) — leaves a regular outer-coordinate density $I$:

$$\chi=\frac{\Pi^2}{\omega^2}\int I(\xi)\,d\xi.$$

Closures A and B differ only in how they reach the poles: B sums one single-pole transform per
harmonic (truncating the tail), A peels the in-range poles from the closed kernel and sums the rest.

The following reduction is **not** closure-general: the non-relativistic $I(p_\perp)$ factors into parallel moments
× perp Bessel weights (§5.1) only because B's $T_n$ cleanly separates $p_\parallel$-polynomials from the
$z(p_\perp)$-Bessel factors. A's complex-order kernel ties $a(p_\parallel)$ to $z(p_\perp)$, forbidding such factorization.

### 5.1 Non-relativistic

In $(p_\perp,p_\parallel)$ coordinates

$$
I = I(p_\perp)=\sum_n I_n(p_\perp)+\mathbf e_\parallel\mathbf e_\parallel\,I_B(p_\perp),
$$

$$
I_n(p_\perp)=2\pi\!\int dp_\parallel\;
\frac{p_\perp\mathcal U\; 𝓣_n}{\omega-k_\parallel p_\parallel-n\Omega},
\qquad
I_B(p_\perp)=2\pi\!\int dp_\parallel\;\Bigl(p_\perp p_\parallel\,\partial_\parallel f_0-p_\parallel^2\,\partial_\perp f_0\Bigr).
$$

With $\gamma=1$, three things happen: the pole $\zeta_n=(\omega-n\Omega)/k_\parallel$ loses its $\gamma$-dependence,
the Bessel argument $z=k_\perp p_\perp/\Omega$ depends on $p_\perp$ alone, and $𝓣_n$
is polynomial in $p_\parallel$. Now $p_\perp$ and $p_\parallel$ decouple, so

Unlike $I(\gamma)$, $I(p_\perp)$ **does** factor into parallel moments and perp Bessel weights.

**Parallel Landau moments.**
With the recurring Landau combination ($p_\perp\mathcal U=(\omega-k_\parallel p_\parallel)\,\partial_\perp f_0 +k_\parallel p_\perp\,\partial_\parallel f_0$), the moment of the full numerator $p_\perp\mathcal U 𝓣_n$ contains the following base form:

$$
D^m(p_\perp)\equiv\int\!\frac{p_\parallel^{\,m}\,p_\perp\mathcal U}{\omega-k_\parallel p_\parallel-n\Omega}\,dp_\parallel
=\omega\,M^m_F-k_\parallel M^{m+1}_F+k_\parallel p_\perp\,M^m_T,
$$

where

$$
M^m_F(p_\perp)\equiv\int\frac{p_\parallel^{\,m}\,\partial_\perp f_0}{\omega-k_\parallel p_\parallel-n\Omega}\,dp_\parallel,
\qquad
M^m_T(p_\perp)\equiv\int\frac{p_\parallel^{\,m}\,\partial_\parallel f_0}{\omega-k_\parallel p_\parallel-n\Omega}\,dp_\parallel,
$$

$$
M^m_F=-\frac1{k_\|}\mathcal C[p_\|^m\,\partial_\perp f_0](\zeta_n),\qquad
M^m_T=-\frac1{k_\|}\mathcal C[p_\|^m\,\partial_\| f_0](\zeta_n).
$$

And the Bessel factors leave the $p_\parallel$ integral untouched. Note without separability, the moments are recomputed at every $p_\perp$, otherwise it would pull $f_\perp(p_\perp)$ out of $M^m$ and kill the outer integral.

Assembled, the resonant harmonic block reads

$$
I_n(p_\perp)=2\pi
\begin{pmatrix}
p_\perp^2 R_n^2\,D_0 & i\,p_\perp^2 R_n J_n'\,D_0 & p_\perp R_n J_n\,D_1\\[2pt]
-i\,p_\perp^2 R_n J_n'\,D_0 & p_\perp^2 J_n'^2\,D_0 & -i\,p_\perp J_nJ_n'\,D_1\\[2pt]
p_\perp R_n J_n\,D_1 & i\,p_\perp J_nJ_n'\,D_1 & J_n^2\,D_2
\end{pmatrix},
$$

**The 33 entry (summed).** Per harmonic the resonant $2\pi J_n^2 D_2$ reaches the top moments
$M^3_F,M^2_T$, which carry $\zeta_n$ and diverge as $k_\parallel\!\to\!0$; the divergence cancels only after summing against $I_B$ (using $\Sigma_n J_n^2=1$). The finite `e∥e∥` density is

$$
I_{33}
=\sum_n 2\pi J_n^2 D_2+I_B
=\sum_n 2\pi J_n^2\Bigl[\,n\Omega\,M^2_F+(\omega-n\Omega)\,p_\perp\,M^1_T\,\Bigr].
$$

All six distinct entries draw on the five parallel moments
$\{M^0_F,M^1_F,M^2_F,M^0_T,M^1_T\}$ at the slice $p_\perp$ and the six symmetric
Bessel bilinears of the regular triple $\{R_n,J_n',J_n\}$ at $z=k_\perp p_\perp/\Omega$.

For relativistic f₀ see [relativistic.md](relativistic.md).

### 5.2 Non-relativistic Separable `f₀=f∥·f⊥`

Now `∂∥f₀=f∥′f⊥`, `∂⊥f₀=f∥f⊥′`, so the `p⊥`-slice pulls out of `𝒞`:

    𝒞[ p^m ∂∥f₀(·,p⊥) ](ζ_n) = f⊥(p⊥) · 𝒞[ p^m f∥′ ](ζ_n),
    𝒞[ p^m ∂⊥f₀(·,p⊥) ](ζ_n) = f⊥′(p⊥) · 𝒞[ p^m f∥ ](ζ_n).

So the `p⊥` integral peels `f⊥` out of every parallel moment, and `X_n` **fully factors** into two independent 1-D primitives — the parallel `𝒞` (§5.3) and the perpendicular `P⊥` (§5.4):

    X_n = Σ_terms ( parallel 𝒞-moment of f∥ ) × ( perp P⊥-moment of f⊥ ).

### 5.3 Parallel Cauchy transform and resonance moments

The parallel primitive is the Landau-causal Cauchy transform over `supp g` (finite per cell, infinite for a Gaussian; `Im ζ>0` physical):

$$
\mathcal C[g](\zeta) \equiv \int \frac{g(p)}{p-\zeta}\,dp,
\qquad \zeta_n=\frac{\omega-n\Omega}{k_\parallel}.
$$

The surrounding derivation uses this same `𝒞`. “Cauchy transform” is the precise
off-real-axis name; its real-axis boundary value is a Hilbert transform plus the Plemelj
residue, so the code field `hilbert` is a (retained) misnomer for the off-axis object. The actual resonance moments are the single
indexed functional

$$
\mathcal C_m[g](\zeta_n)
\equiv -\frac1{k_\parallel}\mathcal C[p^m g](\zeta_n)
=\int\frac{p^m g(p)}{\omega-k_\parallel p-n\Omega}\,dp.
$$

For a separable distribution, the two inputs are `g=f∥` and `g=f∥′`:

$$
M_F^m\equiv\mathcal C_m[f_\parallel],\qquad
M_T^m\equiv\mathcal C_m[f_\parallel'].
$$

Thus `F` and `T` are implementation-facing labels, not distinct transforms. No universal
symbol exists for a distribution-agnostic Cauchy moment; $\mathcal C_m[g]$ states its
meaning directly. Reserve $Z_m$ for Maxwellian-only moments. Use `M_F^m,M_T^m` only where
matching the assembler or code fields helps. The `1/k∥` and sign come from
`ω−k∥p−nΩ = −k∥(p−ζ)`. Both are evaluations of the _same_ $\mathcal C$, dispatched on
how `g` is represented:

| `g`                              | $\mathcal C[g]$                                  |
| -------------------------------- | ------------------------------------------------ |
| normalized Gaussian `e^{−p²}/√π` | `Z(ζ)`, plasma dispersion function               |`                                        |
| piecewise polynomial             | per cell `∫q dp + P(ζ)·log((p_{i+1}−ζ)/(p_i−ζ))` |


**Piecewise-poly cell.** Synthetic-divide `P(p)=q(p)(p−ζ)+P(ζ)`; then
`∫P/(p−ζ)dp = ∫q dp + P(ζ)·log((p_{i+1}−ζ)/(p_i−ζ))`. **Branch-cut invariant**: one
complex `log` of the _ratio_ (not a difference of logs) keeps the continuation
single-valued as `Im ζ→0`. Landau continuation to the growing sheet (`Im ζ<0`, `Re ζ` in
cell) adds `2πi·p(ζ)`.

### 5.4 Perpendicular primitive `P⊥`

The perpendicular primitive collects the `p⊥` integrals (`2π` = azimuthal factor of
`d³p`). The Bessel argument `z=k⊥p⊥/Ω₀` runs with `p⊥` (`Jn′≡∂Jn`). Index the three
Bessel bilinears by `j=0,1,2`, `Wⱼ∈{Jₙ², JₙJₙ′, Jₙ′²}`; the `p⊥`-power is fixed by `j` and
the density (`T_n` entries, §3B). The **perp moments**:

$$
P_j \equiv 2\pi\!\int W_j(z)\,f_\perp\,p_\perp^{\,j+1}\,dp_\perp,\qquad
P_j^\partial \equiv 2\pi\!\int W_j(z)\,f_\perp'\,p_\perp^{\,j}\,dp_\perp,
$$

(`∂` superscript = the `f⊥′` slice).

Dispatched on `f⊥`:

| `f⊥`                 | moments                                                         |
| -------------------- | --------------------------------------------------------------- |
| Gaussian             | `Γ_n(λ)=Iₙ(λ)e^{−λ}`, `λ=(k⊥ p_th⊥/Ω₀)²/2`, plus recurrences    |
| piecewise polynomial | per-cell Bessel-product power series (Schläfli ₂F₃) `∫pᵈ⁺¹ J J` |
| arbitrary analytic   | direct adaptive quadrature of the six integrals                 |

The Gaussian closed forms `{P₀=Γ_n, P₀^∂=−2Γ_n/p_th², P₁^∂=−(k⊥/Ω₀)Γ_n′, …}` are the
normalized `f⊥` evaluations of exactly these integrals — closed per harmonic, but the
**sum over `n`** of `Γ_n` is still truncated (Gaussian is evaluator B).

### Assembling χ from the two primitives

With both primitives in hand, the `p⊥` integral of the §5.1 block turns each entry into a single
**perp–parallel contraction**: the `f⊥′`-slice `P_j^∂` pairs with the `f⊥`-slice `P_j`,

$$
\Sigma_j^{(m)}\;\equiv\;P_j^\partial\bigl(\omega\,𝖬_F^{\,m}-k_\|𝖬_F^{\,m+1}\bigr)+P_j\,k_\|𝖬_T^{\,m},
$$

where `𝖬_F^m=−\tfrac1{k_\|}𝒞[p_\|^m f_\|](ζ_n)` and `𝖬_T^m=−\tfrac1{k_\|}𝒞[p_\|^m f_\|'](ζ_n)` are the
_pure-parallel_ moments (`Z,Z'` of the bi-Maxwellian, [Maxwellian.md]). The `f⊥`/`f⊥′` split of `D_m`
(§5.1) is what lets `P_j` and `P_j^∂` factor out of the `p⊥` integral. Using `R_n=(n/z)J_n`
(so `p⊥R_n=(n/β)J_n`, `β=k⊥/Ω₀`), the §5.1 matrix integrates to

$$
\chi=\sum_n\chi_n,\qquad
\chi_n=\frac{\Pi^2}{\omega^2}\begin{pmatrix}
\tfrac{n^2}{\beta^2}\Sigma_0^{(0)} & i\tfrac n\beta\Sigma_1^{(0)} & \tfrac n\beta\Sigma_0^{(1)}\\[3pt]
-i\tfrac n\beta\Sigma_1^{(0)} & \Sigma_2^{(0)} & -i\,\Sigma_1^{(1)}\\[3pt]
\tfrac n\beta\Sigma_0^{(1)} & i\,\Sigma_1^{(1)} & \Sigma_{zz}
\end{pmatrix},\quad
\Sigma_{zz}=n\Omega_0\,P_0^\partial𝖬_F^2+(\omega-n\Omega_0)\,P_0𝖬_T^1
$$

This is where `P⊥` enters χ — linearly, once per entry. A new VDF is then **only** a new `{P_j,P_j^∂}` recipe; the
wiring above is fixed (Gaussian→`Γ_n`, gyro-ring, ring-beam all reuse it).
