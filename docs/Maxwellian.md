## Drifting bi-Maxwellian closes both primitives

Take the normalized drifting bi-Maxwellian `f₀ = f⊥ ⋅ f∥`,

$$
f_\|(p)=\frac{e^{-(p-v_d)^2/p_{th\|}^2}}{\sqrt\pi\,p_{th\|}},\qquad
f_\perp(p)=\frac{e^{-p^2/p_{th\perp}^2}}{\pi\,p_{th\perp}^2}.
$$

Two master functions close the primitive integrals, and in each the derivative form needs _no_ new integral (since `f'=−2(p−v_d)f/p_th²`).

**Parallel `𝒞` → plasma dispersion function `Z`.** Rescale `u=(p−v_d)/p_th∥`; the pole
maps to

$$
\xi_n=\frac{\zeta_n-v_d}{p_{th\|}}=\frac{\omega-k_\| v_d-n\Omega_0}{k_\| p_{th\|}},
\qquad
\mathcal C[f_\|](\zeta_n)=\frac1{p_{th\|}}\,Z(\xi_n).
$$

The derivative collapses onto the same family via `f∥′=−2(p−v_d)f∥/p_th∥²`:

$$
\mathcal C[f_\|'](\zeta_n)=-\frac2{p_{th\|}^2}\bigl[1+\xi_n Z(\xi_n)\bigr]=\frac1{p_{th\|}^2}Z'(\xi_n),
$$

giving

$$
M^0_F=-\frac{Z(\xi_n)}{k_\| p_{th\|}},\qquad
M^0_T=-\frac{Z'(\xi_n)}{k_\| p_{th\|}^2}.
$$

**Perpendicular `P⊥` → ring sum $Γ_n$**. Weber's second exponential integral,

$$
\int_0^\infty e^{-p^2/p_{th\perp}^2}J_n(\beta p)^2\,p\,dp
=\frac{p_{th\perp}^2}{2}\,e^{-\lambda}I_n(\lambda),
\qquad \beta=\frac{k_\perp}{\Omega_0},\quad
\lambda\equiv\frac{(\beta\,p_{th\perp})^2}{2},
$$

collapses the base moment to the modified-Bessel ring sum (`z=k⊥p/Ω₀`):

$$
P_0=2\pi\!\int_0^\infty J_n(z)^2 f_\perp\,p\,dp=e^{-\lambda}I_n(\lambda)\equiv\Gamma_n(\lambda).
$$

Again the derivative needs no new integral — `f⊥′=−2p f⊥/p_th⊥²` gives

$$
P_0^\partial=2\pi\!\int_0^\infty J_n(z)^2 f_\perp'\,dp=-\frac{2}{p_{th\perp}^2}\,P_0=-\frac{2\,\Gamma_n(\lambda)}{p_{th\perp}^2},
$$

and the mixed combo from differentiating Weber w.r.t. `β`: `P₁^∂=−(k⊥/Ω₀)Γ_n′`, with
`Γ_n′(λ)=½(Γ_{n−1}+Γ_{n+1})−Γ_n`.

The perp moments $\{P_j,P_j^\partial\}$ ($j=0,1,2$) are the _only_ `f⊥`-dependent input to χ: the shared
perp–parallel contraction $\Sigma_j^{(m)}$ (`derivation.md` §5.4) pairs each with the parallel $Z$-moments
identically for every model below. So each ring variant is just a new $\{P_j,P_j^\partial\}$ recipe.

## Ring generalization: the shifted-perp Maxwellian

The clean closure survives a **perpendicular ring shift** $p_r$. Averaging of a 2-D Maxwellian _beam_ centred at $\mathbf p_r=p_r\hat{\mathbf x}$ over gyrophase $\phi$
(with $\mathbf p_\perp=p_\perp(\cos\phi,\sin\phi)$, so $|\mathbf p_\perp-\mathbf p_r|^2=p_\perp^2+p_r^2-2p_\perp p_r\cos\phi$):

$$
\begin{aligned}
f_\perp^{ring}(p_\perp) &= \Big\langle\frac{e^{-|\mathbf p_\perp-\mathbf p_r|^2/p_{th\perp}^2}}{\pi p_{th\perp}^2}\Big\rangle_\phi \\
&=\frac{e^{-(p_\perp^2+p_r^2)/p_{th\perp}^2}}{\pi p_{th\perp}^2}\;
\frac1{2\pi}\!\int_0^{2\pi}\!e^{\,2p_r p_\perp\cos\phi/p_{th\perp}^2}d\phi\\
&=\frac1{\pi p_{th\perp}^2}
\exp\!\Big(-\frac{p_\perp^2+p_r^2}{p_{th\perp}^2}\Big)\,
I_0\!\Big(\frac{2p_r p_\perp}{p_{th\perp}^2}\Big)
\end{aligned}
$$

the last term is the integral representation of the modified Bessel function of the first kind of order 0: $I_0(x)=\frac1{2\pi}\int_0^{2\pi}e^{x\cos\phi}d\phi$.

$$
\qquad \Lambda_r\equiv\beta\,p_r=\frac{k_\perp p_r}{\Omega_0}.
$$

Equivalently it is the cold ring $\delta(p_\perp-p_r)/(2\pi p_r)$ **2-D-convolved** with the Gaussian.
The two views matter: because $f_\perp^{ring}$ is a _bona-fide 2-D-plane_ Gaussian (just gyro-averaged),
its perp moments are full-plane integrals → Weber closes them. For $p_r\gg p_{th\perp}$,
$I_0(x)\sim e^x/\sqrt{2\pi x}$ ⇒ $f_\perp^{ring}\propto e^{-(p_\perp-p_r)^2/p_{th\perp}^2}$, the
shifted Gaussian — but only _gyrotropized_; the literal magnitude-Gaussian (last section) is the **$\phi=0$ slice**, and loses the 2-D-plane structure (hence no Weber closure).

**Base moment closes by a second Weber integral.** The ring's base perp moment is the
direct analogue of the Maxwellian's $P_0=\Gamma_n(\lambda)$:

$$
\Gamma_n^{ring}(\lambda,\Lambda_r)\equiv P_0^{ring}
=2\pi\!\int_0^\infty J_n(z)^2\,f_\perp^{ring}\,p_\perp\,dp_\perp .
$$

It closes the same way: sum the perp generating function
$Σ_n J_n(z)^2 e^{-in\chi}=J_0(2z\sin\tfrac\chi2)$ against $f_⊥^{ring}$ and apply Weber's
_same-order_ integral $∫_0^∞ t e^{-pt^2}I_0(at)J_0(ct)\,dt=\tfrac1{2p}e^{(a^2-c^2)/4p}J_0(ac/2p)$,
which collapses everything (the $e^{-p_r^2/\dots}$ prefactor cancels):

$$
\boxed{\;\sum_n \Gamma_n^{ring}\,e^{-in\chi}
= e^{-\lambda(1-\cos\chi)}\,J_0\!\big(2\Lambda_r\sin\tfrac\chi2\big)\;}
$$

— the Maxwellian generating function times a **cold-ring `J_0`**. Re-expand the `J_0` with the
same identity ($J_0(2\Lambda_r\sin\tfrac\chi2)=Σ_m J_m(\Lambda_r)^2 e^{im\chi}$) and read off the
$e^{-in\chi}$ coefficient:

$$
\boxed{\;\Gamma_n^{ring}(\lambda,\Lambda_r)=\sum_{m=-\infty}^{\infty}
J_m(\Lambda_r)^2\,\Gamma_{n+m}(\lambda)\;}
$$

A **discrete convolution**: the cold-ring Bessel spectrum $J_m^2(\Lambda_r)$ (peaked at `m≈\Lambda_r`,
dead for `|m|>\Lambda_r`) smears the thermal `Γ_n(λ)`. Both limits are exact:

- `Λ_r→0`: `J_m(0)^2=\delta_{m0}` ⇒ `Γ_n^{ring}=Γ_n` (Maxwellian).
- `λ→0`: `Γ_{n+m}(0)=\delta_{n+m,0}` ⇒ `Γ_n^{ring}=J_n^2(\Lambda_r)` (cold ring).

The `λ`-derivative commutes with the convolution, so `Γ_n′^{ring}=½(Γ_{n-1}^{ring}+Γ_{n+1}^{ring})-Γ_n^{ring}`
holds verbatim. With the substitution `Γ_•→Γ_•^{ring}` (explicit `n`, `λ`, `β` kept) the
**`n`-weighted and `Γ′`-type tensor entries close unchanged** — `Rn²=-n²Γ_n/λ`, `RnJn=-βnΓ_n/λ`,
`RnJn′=-nΓ_n′`, `Jn′Jn=-βΓ_n′` — because those are exact Fourier identities of `e^{-in\chi}`
(differentiation in `χ` ⇒ factor `n`; differentiation in `λ` ⇒ `Γ′`; the cold-ring `J_0`
carries no `λ`).

**The remaining entries close through one extra moment.** The substitution as-is breaks for the
entries whose Maxwellian form used `f_⊥'=-2p f_⊥/p_{thperp}^2` (the ring's `f_⊥'` carries an extra
`I_1`) and for the `F`-slice `Jₙ′²`. Geometry first collapses the tensor: with
`p_⊥R_n=(n/β)J_n`, every entry is built from **three** Bessel structures `{Jₙ², p_⊥JₙJₙ′, p_⊥²Jₙ′²}`
in each slice, and the `(n/β)` factors are exact. Two structures per slice reduce to `Γ_n^{ring}`
and its `(λ,Λ_r)` partials directly; the third needs the single `v⊥²`-moment

$$
K_n\equiv 2\pi\!\int_0^\infty p_\perp^3\,f_\perp^{ring}\,J_n(z)^2\,dp_\perp .
$$

`K_n` is **also** closed. Its `v⊥³`-Weber is `−∂_p` of the _same-order_ `I_0 J_0` integral, and the
generating function reduces — using $−2\psi J_1(\psi)=2\Lambda_r\,\partial_{\Lambda_r}J_0(\psi)$,
$\psi=2\Lambda_r\sin\tfrac\chi2$ — to a clean combination of `Γ_n^{ring}` and its partials:

$$
\boxed{\;K_n=(2\sigma^2+p_r^2)\,\Gamma_n^{ring}+2\sigma^2\lambda\,\partial_\lambda\Gamma_n^{ring}
+2\sigma^2\Lambda_r\,\partial_{\Lambda_r}\Gamma_n^{ring}\;},\qquad \sigma^2=\tfrac12 p_{th\perp}^2,
$$

with the clean convolutions $\partial_\lambda\Gamma_n^{ring}=Σ_m J_m^2\,\Gamma_{n+m}'$ and
$\partial_{\Lambda_r}\Gamma_n^{ring}=Σ_m 2J_mJ_m\,\Gamma_{n+m}$. So **all 12 perp moments are
closed**, sharing the one cold-ring spectrum `J_m(\Lambda_r)^2`. (The lone `n=0` value of the
`f_⊥'`-base `Jₙ²` entry is never needed — it enters `χ_{zz}` only through `n\Omega`, which vanishes.)

This is the cold⊕Gaussian principle: momentum `\mathbf p=\mathbf p_c+\mathbf p_g` (cold-ring vector
⊕ independent Gaussian); scalar `J_n^2` moments convolve with `J_m^2`, and the `v⊥²` insertion adds
the `p_r^2` (cold) and `\partial_{\Lambda_r}` (cross) pieces above.

> Implemented as `Maxwellian(; vr=…)` (`vr=0` keeps the bit-identical fast path); validated against
> `SeparableVDF` to ~1e-10, ~30–45× faster than the quadrature path.

## Literal shifted-Gaussian ring-beam (the `exp[−(p_⊥−p_r)²/p_{th⊥}²]` form)

A widely-used ring-beam writes the perp factor as a Gaussian in the **magnitude** `p_⊥`,
$f_⊥(p_⊥)=e^{−(p_⊥−p_r)²/p_{th⊥}²}/N_⊥$,
$N_⊥=2π\!\int_0^\infty p_⊥ e^{−(p_⊥−p_r)²/p_{th⊥}²}dp_⊥=\pi p_{th⊥}^2 A_e$,
$A_e=e^{−(p_r/p_{th⊥})^2}+\sqrt\pi(p_r/p_{th⊥})\,\mathrm{erfc}(−p_r/p_{th⊥})$. This is **not** the
gyro-averaged ring above:
$e^{−(p_⊥−p_r)²/p_{th⊥}²}=e^{−(p_⊥^2+p_r^2)/p_{th⊥}^2}\sum_k I_k(2p_r p_⊥/p_{th⊥}^2)$, an infinite tower
of `I_k`-rings of which the closed $\Gamma_n^{ring}$ is only `k=0`. The `k≠0` terms pair `I_k` with `J_0`
(different order), so **no finite Bessel closure exists** — magnitude-positivity forces an
`erfc`/parabolic-cylinder transcendental (already visible in `A_e`).

### Route A — exact-shift parabolic-cylinder derivation

This is a new `{P_j,P_j^∂}` recipe for the §5.4 perpendicular primitive (the §5.4 assembly into χ is
unchanged; the parallel `Z`-moments are the drifting-Gaussian ones above). With `σ²≡½p_{th⊥}²`,
`β=k_⊥/Ω₀`, `z=βp_⊥`, the moments are (§5.4)

$$
P_j=2\pi\!\int_0^\infty W_j(z)\,f_⊥\,p_⊥^{\,j+1}\,dp_⊥,\qquad
P_j^\partial=2\pi\!\int_0^\infty W_j(z)\,f_⊥'\,p_⊥^{\,j}\,dp_⊥,\qquad
f_⊥'=-\tfrac{p_⊥-p_r}{\sigma^2}f_⊥,
$$

`W_j∈{J_n²,J_nJ_n',J_n'²}`.

**(i) Bessel-product series.** Reduce $J_n'=\tfrac12(J_{n-1}-J_{n+1})$, so every $W_j$ is a sum of
products $J_\mu J_\nu$, $\mu,\nu\in\{n-1,n,n+1\}$. Each product is an entire power series (the Schläfli
$_2F_3$ already used on the piecewise-poly path, §5.4):

$$
J_\mu(z)J_\nu(z)=\sum_{l\ge0}e^{\mu\nu}_l\,z^{\mu+\nu+2l},\qquad
e^{\mu\nu}_0=\frac1{\mu!\,\nu!\,2^{\mu+\nu}},\qquad
\frac{e^{\mu\nu}_{l+1}}{e^{\mu\nu}_l}=\frac{-(\mu+\nu+2l+1)(\mu+\nu+2l+2)}{4(l+1)(\mu+\nu+l+1)(\mu+l+1)(\nu+l+1)}
$$

(coded by this ratio — no factorials; negative orders fold by $J_{-m}=(-1)^mJ_m$).

**(ii) Term-by-term integration → shifted moments.** With $z=\beta p_⊥$, define the master integral

$$
\mathcal S_q(\mu,\nu)\equiv\int_0^\infty p_⊥^{\,q}\,e^{-(p_⊥-p_r)^2/p_{th⊥}^2}J_\mu(\beta p_⊥)J_\nu(\beta p_⊥)\,dp_⊥
=\sum_{l\ge0}e^{\mu\nu}_l\,\beta^{\mu+\nu+2l}\,\mathcal E_{q+\mu+\nu+2l},
$$

where the only velocity integral left is the **shifted (parabolic-cylinder) moment**
$\mathcal E_k\equiv\int_0^\infty p_⊥^{\,k}e^{-(p_⊥-p_r)^2/p_{th⊥}^2}dp_⊥$.

**(iii) $\mathcal E_k$ recurrence (erfc seed).** Integrate $\partial_{p}\,e^{-(p-p_r)^2/p_{th⊥}^2}=-\tfrac{p-p_r}{\sigma^2}e^{-(p-p_r)^2/p_{th⊥}^2}$
by parts against $p^{k-1}$ (boundary vanishes for $k\ge2$):
$-\tfrac1{\sigma^2}(\mathcal E_k-p_r\mathcal E_{k-1})=-(k-1)\mathcal E_{k-2}$, i.e.

$$
\boxed{\;\mathcal E_k=p_r\,\mathcal E_{k-1}+(k-1)\sigma^2\,\mathcal E_{k-2}\;},\qquad
\mathcal E_0=\tfrac{\sqrt\pi}{2}\,p_{th⊥}\,\mathrm{erfc}(-p_r/p_{th⊥}),\quad
\mathcal E_1=p_r\mathcal E_0+\sigma^2 e^{-(p_r/p_{th⊥})^2}=\sigma^2 A_e .
$$

So $\mathcal E_1=\sigma^2 A_e$ and $N_⊥=2\pi\mathcal E_1=\pi p_{th⊥}^2 A_e$.

**(iv) Assemble the moments.** Insert (i) into `P_j`; the `p_⊥^{\,j+1}` measure raises `q`:

$$
P_0=\frac{\mathcal S_1(n,n)}{\mathcal E_1},\quad
P_1=\frac{\mathcal S_2(n,n{-}1)-\mathcal S_2(n,n{+}1)}{2\,\mathcal E_1},\quad
P_2=\frac{\mathcal S_3(n{-}1,n{-}1)-2\mathcal S_3(n{-}1,n{+}1)+\mathcal S_3(n{+}1,n{+}1)}{4\,\mathcal E_1}.
$$

The $f_⊥'$ slice replaces the weight $p_⊥^{\,q}f_⊥\to p_⊥^{\,q}f_⊥'$, i.e.
$\mathcal S_q\to\tfrac1{\sigma^2}\big(p_r\,\mathcal S_q-\mathcal S_{q+1}\big)$ at the matching $q$:

$$
P_0^\partial=\frac{p_r\,\mathcal S_0(n,n)-\mathcal S_1(n,n)}{\sigma^2\mathcal E_1},\quad
P_1^\partial=\frac{p_r\,\Delta_1-\Delta_2}{2\sigma^2\mathcal E_1},\ \dots,\quad
\Delta_q\equiv\mathcal S_q(n,n{-}1)-\mathcal S_q(n,n{+}1).
$$

`p_r=0` collapses $\mathcal E_k=\tfrac12\Gamma(\tfrac{k+1}2)p_{th⊥}^{k+1}$ and $P_0=\Gamma_n(\lambda)$,
$P_0^\partial=-\Gamma_n/\sigma^2$ — the bi-Maxwellian (§ above). Everything is a single `l`-series of the
erfc-seeded $\mathcal E_k$; $~\Lambda_r+\sqrt\lambda$ terms, **independent of $p_r/p_{th⊥}$**
($\Lambda_r=\beta p_r$, $\lambda=\sigma^2\beta^2$).

### The cancellation wall in evaluating $P_j$

Evaluating $P_j$ means summing the series $\mathcal S_q=\sum_l e^{\mu\nu}_l\,\beta^{\mu+\nu+2l}\mathcal E_{q+\mu+\nu+2l}$
of step (ii). The moments $\mathcal E_k$ peak where the ring sits ($p_⊥\approx a$, so $\mathcal E_k\sim a^k$),
hence the term scales as $e^{\mu\nu}_l(\beta a)^{\mu+\nu+2l}$ — the sum **is** $J_\mu(\Lambda_r)J_\nu(\Lambda_r)$
up to scale, $\Lambda_r=\beta a$. Its terms **alternate** (the $(-1)^l$ in $e^{\mu\nu}_l$) and swell to the
all-positive _modified_-Bessel envelope — exactly the $I_k$-tower the literal ring sums to — before
collapsing to the true, small $P_j$:

$$
\max_l\big|\text{term}_l\big|\ \sim\ I_\mu(\Lambda_r)I_\nu(\Lambda_r)\ \sim\ \frac{e^{2\Lambda_r}}{2\pi\Lambda_r},
\qquad P_j=O(1/\Lambda_r).
$$

So computing $P_j$ loses $\approx\log_{10}e^{2\Lambda_r}\approx0.8\,\Lambda_r$ digits to cancellation (verified
on $J_0^2$, $z\equiv\Lambda_r$):

| $z$ | $\max\text{term}/ J_0^2 $ | digits lost | $\varepsilon\cdot$ ratio |
| --- | ------------------------- | ----------- | ------------------------ |
| 5   | $5.9\times10^3$           | 3.8         | $1\times10^{-12}$        |
| 10  | $2.3\times10^7$           | 7.4         | $5\times10^{-9}$         |
| 14  | $8.6\times10^{10}$        | 10.9        | $2\times10^{-5}$         |
| 20  | $5.5\times10^{15}$        | 15.7        | $\sim1$ (total loss)     |

(The Gaussian tail samples $p_⊥>a$, pushing the controlling $z$ a little above $\Lambda_r$.) So Route A
holds ≤1e-9 for $\Lambda_r\lesssim8$, ~5 digits at $\Lambda_r\approx12$, total loss by $\Lambda_r\approx18$.
**This is the same $I_k$ tower that blocks closure**: the $I_n$ envelope those rings sum to is precisely the
$e^{2\Lambda_r}$ magnitude the alternating $\mathcal S_q$ must annihilate — no-closure and finite precision are
one obstruction.

**Beyond the wall (`Λr≳10`).** Two routes were tried; neither is a clean win, so `SeparableVDF` is the
fallback:

- **Large-argument asymptotic (closed form) — _fails_ for χ.** Using
  $J_μ(z)J_ν(z)\sim\tfrac1{\pi z}[\cos\tfrac{(\mu-\nu)\pi}2+\cos(2z-\tfrac{(\mu+\nu)\pi}2-\tfrac\pi2)]$,
  each moment becomes a Gaussian moment of `1` and `\cos2\beta p_⊥` (closed, `e^{-2\lambda}`-suppressed
  oscillation = the Bernstein structure). The _low_ harmonics close well (`P_0` to `1e-3` at `Λr=12`), but
  the **dominant harmonics `n≈Λr` sit at the Bessel turning point `z≈n`**, where this asymptotic is invalid
  (it needs the uniform Airy form). The per-moment error grows with `n` and the harmonic sum amplifies it:
  full-`χ` error is `~14–90%` for `Λr=8–40` — not usable.
- **Stable `besselj` at fixed nodes.** A fixed quadrature on the ring weight with library `besselj` (its own
  large-`z` asymptotics ⇒ no cancellation) _is_ accurate, but the integrand oscillates `~(p_r/p_{th⊥})^{-1}Λr`
  times across the ring, so it needs `O(Λr)` nodes — no better than `SeparableVDF`'s adaptive quadrature.

There is no cheap closed form at large `Λr`: the integrand genuinely oscillates `~Λr` times and the turning-
point harmonics dominate. Route A owns the `Λr≲8` fast-exact niche; for `Λr≳8` use `SeparableVDF`.

### Route B (DGH) — for comparison

Expanding the _shift_ instead, $e^{2p_r p_⊥/p_{th⊥}^2}=\sum_j (2p_r/p_{th⊥}^2)^j p_⊥^{\,j}/j!$, gives
$M_n=e^{-(p_r/p_{th⊥})^2}\sum_j \tfrac{(2p_r/p_{th⊥}^2)^j}{j!}\sum_l c_{n,l}\beta^{2n+2l}Q_{2n+2l+1+j}$,
$Q_k=\tfrac12\Gamma(\tfrac{k+1}2)p_{th⊥}^{k+1}$ the unshifted (`erf`-class) moments. A **double** sum;
the `j`-series needs $~2(p_r/p_{th⊥})^2$ terms — slower than quadrature for sharp rings, no accuracy gain.
Both routes are forced through the `erf`/`erfc` family, the computational signature of the no-Bessel-closure proof.

### Accuracy & speed

`Γ_n^{beam}` vs adaptive quadrature (single-moment study):

| `p_r/p_{th⊥}`, `λ` | `\Lambda_r` | A error / `l` | B error / `j` | A         | B     | quad  |
| ------------------ | ----------- | ------------- | ------------- | --------- | ----- | ----- |
| 0.5, 0.5           | 0.5         | 5e-16 / 20    | 3e-15 / 22    | **2.8µs** | 40µs  | 23µs  |
| 2.0, 0.5           | 2.0         | 1e-15 / 24    | 4e-14 / 53    | **2.2µs** | 117µs | 32µs  |
| 5.0, 0.5           | 5.0         | 3e-13 / 32    | 2e-13 / 139   | **2.2µs** | 376µs | 55µs  |
| 5.0, 2.0           | 10          | 2e-7 / 58     | 6e-9 / 140    | **2.3µs** | 657µs | 106µs |

Full-tensor `χ` (`GaussianRing` vs `SeparableVDF`): ~1e-9 and **~400–470× faster** for `\Lambda_r\lesssim8`,
breaking (≈30%) by `\Lambda_r\approx12`.

> Implemented as `GaussianRing(; vth_para, vth_perp, vd, vr)` (Route A, full tensor); `vr=0` reduces to the
> bi-Maxwellian. Past `\Lambda_r\approx10` use `SeparableVDF`.
