#import "@preview/cetz:0.4.2"
#set page(paper: "a4", margin: (x: 20mm, y: 20mm), numbering: "1")
#set text(size: 10pt)
#set par(justify: true)
#set heading(numbering: "1.1")

#text(16pt)[*Relativistic response and continuation*]

= Overview

Dispersion-root tracing evaluates the response at damped frequencies,
$"Im" omega < 0$. This note derives the continuation used for relativistic, oblique, damped (superluminal) modes.
The requested sheet is reached from the upper half-plane through
the subluminal interval $|"Re" omega| < |k_parallel|$.
Its correction is a sum of residues integrated over transported open paths.
Convergence of the harmonic sum imposes a separate accuracy limit.

= Continuation of moving pole integrals <crossing>

For a holomorphic numerator $g$,

$
  phi(p) = integral_(-w)^(w) (g(x)) / (x - p) dif x
$

is holomorphic off the segment. Continuing the upper branch downward through the segment gives

$
  phi_"cont" (p) = phi(p) + 2 pi i space g(p), quad "Im" p < 0.
$

Now let the pole and integration interval depend on $y$:

$
  I(omega) = integral dif y integral_(-w(y))^(w(y))
  (g(x, y)) / (x - p(y, omega)) dif x.
$

On the real axis, define

$
  Gamma(omega) := {y : p(y, omega)^2 < w(y)^2}.
$

Slice-wise boundary-value gives:
$I(omega plus.minus i 0) = "p.v." I(omega) plus.minus i pi integral_(Gamma(omega)) g(p(y, omega), y) dif y.$
Hence the analytic continuation is

$
  I_"cont" (omega) = I(omega)
  + 2 pi i integral_(Gamma(omega)) g(p(y, omega), y) dif y
$

Off-axis, $Gamma(omega)$ is transported holomorphically: its finite endpoints follow the roots $p^2 = w^2$, while unbounded ends remain in transported sectors.
If $g(x, y)$ is entire and decays at infinity, Cauchy's theorem means that only  endpoints and the connected sector $S(omega)$ matter.
Singular $g$ also requires continuation to record the winding.


#figure(
  {
    include "assets/relativistic-continuation-contour.typ"
  },
  caption: [
    After $omega$ moves off the axis: the endpoint continues holomorphically
    to complex $y_0 (omega)$, and $Gamma(omega)$ is any path from
    it whose tail stays in the transported decay sector $S(omega)$.
  ],
) <fig-transport>

= General relativistic response <problem>

Set $a = k_perp \/ Omega_0$ (subscript 0 omitted for simplicity when clear from context), $u = p_parallel$ and $q = p_perp$,
and write the distribution as $f(gamma, u)$.
With $q^2 = gamma^2 - 1 - u^2$ and $w(gamma) = sqrt(gamma^2 - 1)$,
the harmonic response is

$
  F(omega) = sum_(n = -oo)^(oo) I_n (omega),
  quad
  I_n = integral_1^oo dif gamma integral_(-w(gamma))^(w(gamma)) dif u
  space (cal(U) 𝒯_n) / D_n,
  quad
  D_n := gamma omega - k_parallel u - n Omega_0,
$

where $cal(U)(gamma, u; omega) = omega partial_gamma f + k_parallel partial_u f$ and $𝒯_n (a q, u, q)$ is
the harmonic Bessel tensor, bilinears in $q J_(n plus.minus 1)(a q)$ and $u J_n (a q)$.

On resonance $D_n (u_n) = 0$,

$
  u_n (gamma) = (gamma omega - n Omega) / k_parallel,
  quad
  q_n^2 (gamma) = gamma^2 - 1 - u_n^2(gamma).
$

The pole crosses downward for $k_parallel > 0$ and upward for
$k_parallel < 0$. Combining the crossing orientation with the residue gives

#block(inset: 5pt, fill: rgb("f0f6ff"))[
  $
    Delta I_n (omega) = -(2 pi i) / (|k_parallel|) integral_(Gamma_n(omega))
    cal(U)(gamma, u_n; omega)
    𝒯_n (a q_n, u_n, q_n) dif gamma,
    quad
    F_"germ" = F + sum_n Delta I_n.
  $
]

== Endpoint and contour <endpoint>

The contour starts at the entering root of $q_n^2 = 0$:

$
  gamma_(0n) = (n^2 Omega^2 + k_parallel^2) /
  (n Omega omega + k_parallel Delta_n)
  = (n Omega omega - k_parallel Delta_n) /
  (omega^2 - k_parallel^2),
  quad
  Delta_n = sqrt(n^2 Omega^2 + k_parallel^2 - omega^2).
$
Let $omega = omega_r + i nu$. For $nu < 0$, the radicand $n^2 Omega^2 + k_parallel^2 - omega^2$ avoids the principal-square-root branch cut $(-oo, 0]$,  so this is the branch continued from the real subluminal interval. At finite damping, the endpoint remains finite even while passing around the light-line singularity.

Parameterize contour tail of $Gamma_n$ by $gamma = gamma_(0n) + t e^(i theta)$.
Asymptotically $q_n tilde c gamma$ with $c = sqrt(1 - omega^2 \/ k_parallel^2)$,
$cal(T)_n$ (Bessel bilinears) grows as $e^(2 |"Im"(a q_n)|)$ against the distribution decay $e^(-mu "Re" gamma)$, giving

$
  |"integrand"| tilde e^(-t r(theta)), quad "where" r(theta) := mu cos theta - 2|a||c||sin(theta + arg c)|.
$

Note that $S(omega)$ always contains $theta^* = -arg c$, which eliminates Bessel growth and leaves $mu cos(arg c) > 0$.

== Summing continuations: Harmonic convergence <convergence>

Each transported harmonic is finite; their sum need not converge.
For $omega_r > k_parallel > 0$ and $n Omega < 0$, $gamma_(0n) approx |n Omega| \/ (k_parallel - omega)$ has
$"Re" gamma_(0n) < 0$, with magnitude linear in $|n|$. Large-order Bessel decay
offsets the distribution growth only when

$
  2 eta(x) > rho(omega) quad arrow.l.r.double quad "sum converges as"
  |Delta I_n| tilde exp(|n| (rho(omega) - 2 eta(x)))
$

$
  rho(omega) = (mu |Omega| (omega_r - k_parallel)) /
  (|k_parallel - omega|^2),
  quad
  eta(x) = ln((1 + sqrt(1 - x^2)) \/ x) - sqrt(1 - x^2)
  quad
  x = k_perp |c| \/ (|k_parallel - omega|)
$

for $x<1$, the Bessel sub-turning-point exponent $eta = 0$ for $x >= 1$.
Inside this domain, harmonic shells decay to tolerance; outside, they decrease then grow. Termwise harmonic continuation may therefore diverge even when a continuation of the summed response exists (Paley–Wiener theorem).

= Parametrize by $p_perp$

Although $D_n$ is not rational in $u$, its conjugate factor gives

$
  D_n tilde(D_n) = A u^2 + B u + C, quad
  tilde(D_n) := omega gamma + k_parallel u + n Omega,
$

$
  A = omega^2 - k_parallel^2, quad
  B = -2 k_parallel n Omega, quad
  C = omega^2 (1 + q^2) - n^2 Omega^2.
$

Thus $1\/D_n=tilde(D_n)\/[A(u-u_+)(u-u_-)]$. One root solves $D_n=0$; the other solves $tilde(D_n)=0$ and has zero residue. Keeping both avoids branch-dependent root classification. The residue at either root is evaluated stably as

$
  r_s = (g(u_s) tilde(D_n)(u_s)) / (2 A u_s + B)
$

Its two rationalized roots coalesce at where its discriminant $B^2 - 4 A C = 4 omega^2 (n^2 Omega^2 + (k_parallel^2 - omega^2) (1 + q^2))$ vanishes,

$
  q_"apex"^2 = (n^2 Omega^2) / (omega^2 - k_parallel^2) - 1,
$

the branch point lies below the path for $nu > 0$, touches it at $nu = 0$, and moves above it for $nu < 0$.

This continuation is exact for $omega_r^2 < k_parallel^2$. There the discriminant does not vanish on the outer path, and the physical pole crosses uniformly. Finite support endpoints sit where $f_0 approx 0$, so pole-endpoint collisions carry negligible coefficients.

For $omega_r^2 > k_parallel^2$, inside the resonant band $k_parallel^2 < omega_r^2 < k_parallel^2 + n^2 Omega^2$, this apex reaches the real $q$ path, and the straight real-sliced integral omits the cut discontinuity.

Comment: this formulation is more computationally efficient. The Bessel argument $z = k_perp p_perp\/Omega_0$ is independent of the inner $p_parallel$, so the Bessel tensor can be computed once per slice.

// Branch points occur when boundary roots collide (discriminant zeros) or escape
// through infinity (leading-coefficient zeros). Different paths around these
// points add periods of the residue integrand.

== Why energy is the transport coordinate

In energy coordinates, $u_n$ and $q_n^2$ are polynomial in $gamma$, and $𝒯_n$ is even in $q_n$. For entire $f(gamma, u)$, the integrand is entire in $gamma$ and $q_n^2 = 0$ is a regular endpoint (the artificial singularity at the apex is just a regular interior point of the integrand).

Using $q$ as transport coordinate instead fuses physical $gamma>0$ and mirror $gamma<0$ roots at the apex. It also invites the wrong branch $sqrt(1+q_n^2+u_n^2)=sqrt(gamma^2) = plus.minus gamma$ when evaluating at complex momenta. The principal root fails once the transported endpoint enters $Re(gamma)<0$.

= Maxwell-Jüttner fast paths <mj-fast>

For $f = e^(-mu gamma)$ the Swanson gyrophase-time
form $chi prop integral_0^oo e^(-sqrt(R)) K_nu (sqrt(R))"-terms" dif xi$,

$
  R = ((mu Omega - i xi omega)^2 + 2 k_perp^2 (1 - cos xi) + (k_parallel
    xi)^2) / Omega^2
$

- *Parallel ($k_perp = 0$):* $R_0 = A B \/ Omega^2$
  factorizes; per-factor principal square roots are the germ branch along any
  path avoiding $s_(*minus.plus) = -i mu |Omega| \/ (omega minus.plus
    k_parallel)$ (radial cuts). Integrate each folded harmonic on a ray in the
  wedge between the cut angles; the cyclotron-resonant harmonic may need a
  dogleg with sign flips at cut crossings. Certified to $10^(-8)"–"10^(-10)$
  against the corrected López closed form.
- *Oblique $k_perp <= 0.3|Omega|$:* expansion about $R_0$ in $k_perp^2$
  (trig-polynomial assembler over parallel primitives $J(nu, j, q)$, order $M
  <= 3$, per-$J$ `rtol` scaling). Truncation $approx (C k_perp^2 \/
    Omega^2)^(M+1)$, $C approx 5"–"10$.

#block(inset: 8pt, fill: rgb("fff4f0"))[
  *Cliff-weighted validity.* The expansion's
  effective constant multiplies the *cliff strength* of its $|j| >= 2$
  primitives. At strong-cliff points (shallow damping deep in-band, e.g. $mu
  = 2$, $omega = 0.7 - 0.1i$, $k_parallel = 0.5$) those primitives reach
  $e^(28)$ and the series diverges from its first oblique term ($|chi_"exp"|
  = 10^9$ at $k_perp = 0.02$ vs the parallel value $8 times 10^4$ — it does
  not limit to its own $k_perp = 0$ value). At the certified fixtures (A/IC
  regime) the cliff is weak and the $0.3|Omega|$ radius is real.
]

= The same bound in every representation

The domain boundary is a property of the sheet, not of the harmonic representation. The same obstruction appears as:

1. *Time domain.* $chi = integral_0^oo e^(i omega t) C(t) dif t$: damping
  demands a contour tilted upward forever, while the gyration phase confines
  $C$ to a strip $|"Im" t| lt.tilde log(mu Omega \/ k_perp) \/ Omega$. No
  contour does both.
2. *$k_perp^2$ Taylor series* about the parallel problem: coefficients are
  finite, but the $|j| >= 2$ time-harmonic primitives carry monodromy growing
  exponentially at shallow damping — an asymptotic series with
// *cliff-weighted* radius (§#ref(<mj-fast>, supplement: none)).

At $k_perp != 0$ the monodromy content of this sheet is an
infinite reorganizing ladder. Every finite evaluator is exact where its
representation converges, and an optimally-truncated surrogate with a
quantifiable floor beyond.

= Scope and implementation <scope>

Path $Gamma(omega)$ now requires an analytic complex energy-form gradient and an entire integrand. Exponential families satisfy these assumptions; kappa families need unimplemented singularity and winding tracking and extend the formulation to general energy forms.

This note describes the subluminal-germ sheet. It can be exponentially far from the physical boundary value near marginal in-band frequencies.

//  The Maxwell–Jüttner evaluator therefore keeps $|\operatorname{Im}\omega|\leq10^{-4}|\Omega|$ on the straight integral; locate near-real superluminal roots on the real boundary.
