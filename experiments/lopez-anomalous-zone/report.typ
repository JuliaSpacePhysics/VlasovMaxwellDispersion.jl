#set document(title: [The "anomalous zone" of López et al. (2014) is an analytic-continuation artifact])
#set page(margin: (x: 2.2cm, y: 2.2cm), numbering: "1")
#set text(size: 10.5pt)
#set heading(numbering: "1.1")
#show link: set text(fill: blue.darken(30%))
#set table(stroke: 0.4pt + gray)

#context if target() != "html" { text(16pt, weight: "bold", document.title) }

#v(3mm)
#block(inset: (x: 8mm), stroke: (left: 2pt + gray), outset: (y: 2mm))[
  *Claim.* The finite-$omega_r$ descending segment of the Alfvén/ion-cyclotron
  branch ("anomalous zone", $partial omega_r slash partial k < 0$) published for
  relativistic Maxwell–Jüttner pair plasmas by López et al. 2014
  (#link("https://doi.org/10.1063/1.4894679")[PoP 21, 092107], Fig.~1) and
  largely reproduced by ALPS (Verscharen et al. 2018,
  #link("https://doi.org/10.1017/S0022377818000739")[JPP 84, 905840403],
  Fig.~5) is not a plasma mode. In López's calculation it is a zero of a
  *non-holomorphic* function:
  López's closed-form continuation is continuous across the real-$omega$ axis
  but violates the Cauchy–Riemann equations below it. The unique analytic
  continuation instead carries
  the branch structure: the propagating root rises toward the light line
  while a separate purely imaginary ($omega_r = 0$, aperiodic) family carries
  the strong damping.
]

*Setup*: electron–positron pair plasma, both species
Maxwell–Jüttner with $mu = m c^2 slash T$; $omega_(p e)^2 slash Omega_c^2 = 1$
per species; López normalization $x = omega slash Omega_c$,
$y = c k slash Omega_c$, $z = x slash y$, $t = 1 slash y$.
VlasovMaxwellDispersion.jl (VMD) evaluates the factorized parallel L-mode
(`mode = :L` at exactly $k_perp = 0$) throughout. At the tested $mu = 2$ roots, VMD's two
independent relativistic evaluators (Swanson time-integral closed form;
general momentum-grid `CoupledVDF`) agree to
$|Delta omega| lt.tilde 6 times 10^(-4)$, including deep-damped points.

= López's formula (`lopez.jl`)

Write

$
  w_(e,p)(gamma; x,y) = (y^2-x^2) gamma^2 minus.plus 2x gamma - (1+y^2).
$

Then López's L-mode dispersion function (their Eq. 26, in the normalization
above) is

$
  Lambda_(L)(x,y) = 1 - y^2/x^2 - mu/y^2
  + mu^2/(4 K_2(mu) x y^3)
  sum_(j=e,p) integral_1^oo e^(-mu gamma) J_(j)(gamma,z,t)
  w_(j)(gamma;x,y) dif gamma.
$

Each Lorentz factor $gamma$ contributes the resonance integral

$
  J_j (gamma, z, t) = integral_(-S_(1 j))^(S_(2 j)) (dif xi) / (xi - z), quad
  S_(1 e) = sqrt(1 - gamma^(-2)) + t / gamma, quad
  S_(2 e) = sqrt(1 - gamma^(-2)) - t / gamma
$

(positron swaps $S_1 arrow.l.r S_2$). Their Appendix A gives the
"whole-complex-plane" closed form

$
  J_j = 1/2 ln [((R - S_(2 j))^2 + I^2) / ((R + S_(1 j))^2 + I^2)]
  + i [arctan ((S_(2 j) - R) / I) + arctan ((S_(1 j) + R) / I)
    + theta_j (gamma; R, "sign" I)],
$

with $R = "Re" z$ and $I = "Im" z$. In the subluminal region $|R|<1$, the continuation term contains the anomalous zone:

$
  theta_j = pi sigma(I) H_(j)(gamma;R), quad
  H_j = Theta(gamma-gamma_(1j)(R)), quad
  sigma(I) = cases(2 "if" I<0, 1 "if" I=0, 0 "if" I>0),
$

where the real support endpoints are

$
  gamma_(1e)(R) = (R t + sqrt(t^2+1-R^2))/(1-R^2), quad
  gamma_(1p)(R) = (-R t + sqrt(t^2+1-R^2))/(1-R^2).
$

The closed-form $J$ matches direct quadrature for $"Im" z > 0$, and weakly damped roots.

Tracing $Lambda_L = 0$ reproduces the published curves (script `02`): the
$mu = 10$ hump peaking at $omega_r = 0.443$ at $k = 1.7$ (paper: $0.444$ at
$1.65$), descending to $omega_r arrow 0$ by $k approx 2.4$; and the $mu = 2$
hump peaking $omega_r approx 0.125$ near $k approx 0.47$, descending to $0$ by
$k approx 0.6$–$0.65$. Thus the implementation matches their formula; the
disagreement with VMD is an artifact in their formula.

= Non-holomorphy of the continued $Lambda_L$: analytic diagnosis and measurement

Fix $"Im" z < 0$, so $sigma = 2$, set $u = "Re" z$, and let
$H_j (gamma; u)$ denote López's piecewise Heaviside support. With $y$ fixed, the
$theta$ contribution to Eq. 26 is

$
  C_theta (z, macron(z)) =
  (i pi sigma mu^2) / (4 K_2(mu) y^4)
  sum_j integral_1^oo e^(-mu gamma) a_j (gamma, z) H_j (gamma; u) dif gamma,
$

where

$
  a_(e,p)(gamma,z)
  = ((y^2-y^2 z^2) gamma^2 minus.plus 2 y z gamma - (1+y^2))/z.
$

This includes both the $1/x$ prefactor and the $x$, $x^2$ terms in $w_j$;
$a_j$ is holomorphic away from $z=0$. Only $H_j$ depends on $macron(z)$, through
$u=(z+macron(z))/2$. Away from changes in support topology,

$
  (partial C_theta) / (partial macron(z)) =
  (i pi sigma mu^2) / (8 K_2(mu) y^4)
  sum_j integral_1^oo e^(-mu gamma) a_j (gamma,z)
  (partial H_j (gamma;u)) / (partial u) dif gamma.
$

The derivative of $H_j$ is a moving-endpoint boundary term. Below the
axis, $partial Lambda_L slash partial macron(z)$ is nonzero, while both the
$theta$-less expression and VMD satisfy Cauchy–Riemann to finite-difference
accuracy.

López imposes only the *continuity* condition
$lim_(Gamma arrow 0^+) Lambda^+ = lim_(Gamma arrow 0^-) Lambda^-$ (their
Eqs. 19–20). Continuity of boundary values is necessary but not sufficient for
analytic continuation; the missing requirement is holomorphy in $z$. Moving
support through $("Re" z, "sign Im" z)$ introduces the nonzero boundary
terms.

#figure(
  table(
    columns: 4,
    align: (left, center, center, center),
    table.header([point $omega$], [López $Lambda_L$ (with $theta$)], [López without $theta$], [VMD L-mode]),
    [$0.124 + 0.30 i$ (UHP)], [$1.4 times 10^(-7)$], [$1.4 times 10^(-7)$], [$5.7 times 10^(-9)$],
    [$0.124 - 0.06 i$], [*$8.9 times 10^(-2)$*], [$1.0 times 10^(-6)$], [$1.2 times 10^(-7)$],
    [$0.124 - 0.136 i$ (their descent root)], [*$4.1 times 10^(-1)$*], [$5.3 times 10^(-7)$], [$1.4 times 10^(-7)$],
    [$0.124 - 0.25 i$], [*$3.1$*], [$2.0 times 10^(-7)$], [$1.0 times 10^(-7)$],
  ),
  caption: [Holomorphy defect: numerical measured ratio $|partial f \/ partial macron(z)| \/ |partial f \/ partial z|$
    ($approx 0$ for analytic functions), at $k = 0.5$, $mu = 2$. The $theta$-less base is
    locally holomorphic (scripts `04`–`05`).],
)

#figure(
  image("fig_holomorphy.png", width: 88%),
  caption: [Holomorphy-defect maps at $k = 0.5$, $mu = 2$ (script `08`).
    Left: López's continued $Lambda_L$ — the defect reaches $O(1)$ (yellow)
    across the lower-half-plane; their descent root (red #text(fill: red)[✗],
    $0.12415 - 0.13610 i$, a genuine zero of this non-holomorphic function,
    lies inside the defective region.
    Right: the corrected continuation of §4 has measured defect
    $lt.tilde 10^(-5)$ on the same grid; the VMD root (#text(fill: orange)[●],
    $0.16297 - 0.11048 i$) is shared with VMD.],
) <fig-holo>

= Upper-half-plane continuation cross-check

On $"Im" omega > 0$, no continuation is involved and López and VMD agree.
Scripts `03`–`04` fit $1/f$ there with AAA and extrapolate its poles (zeros of
$f$) into the lower half-plane. At $mu=2$, fits to either source select the VMD
root rather than the López descent:

#figure(
  table(
    columns: 3,
    align: (left, left, center),
    table.header([case], [fit source], [extrapolated zero]),
    [$k = 0.5$], [López UHP], [$0.1630 - 0.1105 i$],
    [$k = 0.5$], [VMD UHP], [$0.1630 - 0.1105 i$],
    [$k = 0.6$], [López UHP], [$0.205 - 0.185 i$],
    [$k = 0.6$], [VMD UHP], [$0.207 - 0.184 i$],
  ),
  caption: [AAA cross-check. Direct VMD roots are
    $0.16297 - 0.11048 i$ at $k=0.5$ and $0.20540 - 0.18513 i$ at $k=0.6$.],
)

#figure(
  image("fig_traces.png", width: 100%),
  caption: [Root of each function continued in $k$ ($mu = 2$; script `08`).
    Left: the López-formula root (red) departs from the shared branch near
    $k approx 0.3$ and executes the published descent to $omega_r = 0$; the
    root of VMD's determinant — and of the corrected López formula, which are
    identical to $10^(-5)$ — rises (blue). Right: the same trajectories in the
    complex $omega$ plane; the AAA-extrapolated zero (◆, fit to López's own
    UHP data at $k = 0.5$) coincides with the blue curve at shown precision.],
)

= The fix: an analytic continuation term (script `06`)

The correct continuation keeps the resonant-support endpoint *analytic in*
$z$ instead of freezing it at $"Re" z$. For the subluminal window
$|"Re" z| < 1$:

$
  Lambda_L^"corr" (z) = Lambda_L^(theta"-less") (z)
  + (i pi sigma mu^2)/(4 K_2(mu) x y^3)
  sum_j integral_(gamma_(1 j)(z))^oo
  e^(-mu gamma) w_j(gamma;x,y) dif gamma, quad
  gamma_(1 e)(z) = (z t + sqrt(t^2 + 1 - z^2)) / (1 - z^2),
$

$gamma_(1p)(z)$ replaces $z t$ by $-z t$. The tail integral runs along any
complex path from the endpoint to the positive real ray; its integrand is
entire. Here $sigma=2$ below the axis, while $sigma=1$ denotes the boundary
value on the axis. Results:

- on the axis it equals López to $lt.tilde 3 times 10^(-10)$, preserving continuity condition;
- holomorphy defect $10^(-6)$–$10^(-7)$ in the lower half-plane
  (vs their $0.09$–$3.1$; @fig-holo right panel);

Root polishing from the López descent seeds now converges to the VMD roots:
$0.16297 - 0.11048 i$ at $k=0.5$ and $0.20540 - 0.18512 i$ at $k=0.6$,
with $|Lambda| tilde 10^(-15)$. No corrected root remains at the descent
locations.

== Superluminal extension: following the branch across the light line (script `09`)

The corrected term is not confined to the subluminal window. On the sheet
reached by continuous continuation through the lower half-plane, its branch
functions are analytic along paths avoiding
$z = plus.minus 1$ (the principal $sqrt(t^2+1-z^2)$ cut needs real
$z^2 > 1 + t^2$; the endpoint denominators $1 - z^2$ pole only at
$z = plus.minus 1$), so the same formula followed continuously through the gap
$1 < "Re" z < sqrt(1+t^2)$ — i.e. $k c < omega_r < sqrt(k^2 c^2 + Omega^2)$ —
*is* the unique continuation of the subluminal germ.

The physics it reveals ($mu = 2$): the A/IC branch crosses the light line at
$k approx 1.9$ and stays *slightly superluminal* ($omega_r slash k -> 1.04$ by
$k = 3$) with slowly recovering damping ($gamma = -0.484$ at $k = 1.6$,
$-0.451$ at $3.0$, $-0.348$ at $6.0$). Relativistic cyclotron resonance
survives $v_"ph" > c$: the resonant ellipse in momentum space persists while
$omega_r^2 - k^2 c^2 lt.tilde Omega^2$, fading only as the resonant Lorentz
factors climb into the $e^(-mu gamma)$ tail.

= What the physical dispersion looks like

The anomaly-relevant fact is the separation of two roots. The propagating A/IC
root reaches $gamma approx -0.49$, then rises through the light line; the
aperiodic root remains on $omega_r = 0$ and carries the published damping dive.
They never coalesce. At $mu = 10$ the aperiodic root lies below the published
frame until $k approx 3.05$.

#figure(
  image("fig_fig5_replica.png", width: 84%),
  caption: [VMD reproduction of Verscharen et al. Fig. 5 / López et al.
    Fig. 1 (`docs/src/case-studies/relativistic_pair.jl`): top $mu=2$, bottom
    $mu=10$; × = digitized published curves. Solid curves are propagating
    A/IC-like roots (dashed past the light line), dash-dot curves are O-modes,
    and dotted curves are aperiodic ($omega_r=0$). The solid roots rise where
    the published curves descend; the aperiodic roots remain separate.],
) <fig-replica>

= Why ALPS partially reproduces the same artifact

ALPS does not evaluate López's closed form, yet its Fig. 5 chases the same
descent. The reason is that its damped-side machinery makes the same *class*
of continuation through real-part criteria.

In their scheme: a pole contributes only when its *real* part is deemed inside the
integration domain via $-P_(max, parallel j) lt.eq "Re"(macron(p)_"pole") lt.eq +P_(max, parallel j)$ (their Eq. 3.18), with $macron(p)_"pole" = Gamma omega slash (k_parallel c) - n Omega_(0 j) slash (k_parallel c)$ (their Eq. 3.15);

Their lower-half-plane term selects resonant support using $"Re" z$. As in §2,
their Wirtinger derivative contains moving-boundary terms; absent exact cancellation,
this construction is not holomorphic below the axis.

= Conclusions

The analysis establishes the defect in López's closed form.

1. The anomalous-zone $omega_r$ descent of López et al. 2014 — and its
  Fig.-5 reproduction in Verscharen et al. 2018, whose fitted continuation
  partially chases the same artifact — is a zero of a non-holomorphic function,
  not a plasma mode.
2. Continuing the support endpoint as a complex function of $z$ restores
  Cauchy–Riemann behavior and yields the VMD roots. AAA extrapolation from
  upper-half-plane data corroborates the same roots.
3. In the computed $mu=2$ and $mu=10$ cases, the propagating root rises instead
  of descending to the imaginary axis. Separate aperiodic roots carry the
  strong damping seen in the published curves.

#pagebreak()

= Appendix: The $k -> 0$ ladder <sec-ladder>

On the imaginary axis below the origin ($omega = -i|gamma|$, $sigma = 2$) the
corrected continuation of §4 reads

$
  Lambda_L = Lambda_L^(theta"-less")
  + i pi sigma dot (mu^2) / (4 K_2(mu) x y^3)
  sum_(j = e, p) integral_(gamma_(1 j)(z))^oo e^(-mu gamma) w_j (gamma) dif gamma,
$

with weights $w_(e,p) = (y^2 - x^2) gamma^2 minus.plus 2 x gamma - (1 + y^2)$
and the *analytic* support endpoints $gamma_(1 e, 1 p)(z) =
(plus.minus z t + sqrt(t^2 + 1 - z^2)) slash (1 - z^2)$. Take $y -> 0$ at
fixed $x$ (so $z = x slash y -> oo$). Everything simplifies exactly:

+ *Endpoints.* $gamma_(1 e)(z) -> -1 slash x$ and $gamma_(1 p)(z) -> +1 slash x$
  (with an $O(y)$ real correction). These are precisely the analytic
  continuations of the $k = 0$ resonant Lorentz factor: a particle resonates
  with the wave iff $gamma omega = minus.plus Omega$, i.e.
  $gamma_"res" = minus.plus Omega slash omega$.
+ *Weights.* $w_(e,p) -> -(x gamma plus.minus 1)^2$.
+ *Tails.* The tail integrals become elementary,
  $integral_(minus.plus 1 slash x)^oo e^(-mu gamma) (x gamma plus.minus 1)^2
  dif gamma = (2 x^2 slash mu^3) thin e^(plus.minus mu slash x)$, so the whole
  continuation term collapses to
  $
    C = -i (2 pi x) / (mu K_2(mu) y^3) cosh(mu / x)
    quad arrow.r.double quad
    C|_(x = -i|gamma|) = -(2 pi |gamma|) / (mu K_2(mu) y^3) cos((mu Omega) / (|gamma|)),
  $
  manifestly real on the axis.
+ *Background.* The $theta$-less part stays $O(1)$: expanding its
  finite-support $xi$-integral in $v slash w$, the $O(t^(-1))$ moment cancels
  the $-mu slash y^2$ cold term exactly (via
  $integral_1^oo e^(-mu gamma) gamma sqrt(gamma^2 - 1) dif gamma = K_2 slash mu$),
  leaving $D(|gamma|) = 1 - (mu^2 slash 3 K_2) integral_1^oo e^(-mu gamma)
  gamma (gamma^2 - 1)^(3 slash 2) slash (|gamma|^2 gamma^2 + 1) dif gamma$.

Numerically, $Lambda_L (-i |gamma|, y) = D(|gamma|) + C$ holds with relative
error $tilde y$ (checked at $y = 0.02 -> 0.005$, $mu = 2$ and $10$; script
`09`). The $y^(-3)$ term therefore dominates as $k -> 0$ and pins every axis
zero to

$
  cos((mu Omega) / (|gamma|)) = 0
  quad arrow.r.double quad
  gamma_n = -(2 mu Omega) / (pi (2 n - 1)). qed
$

The physics is transparent: continuing the resonance to an aperiodic
$omega = -i |gamma|$ makes the resonant energy *imaginary*,
$gamma_"res" = i Omega slash |gamma|$, so the Jüttner factor
$e^(-mu gamma_"res")$ becomes a pure phase $e^(i mu Omega slash |gamma|)$.
Electron and positron contribute conjugate phases, and their interference —
the cosine — vanishes at odd quarter-periods. The aperiodic ladder is an
interference pattern of the analytically continued cyclotron resonance in the
relativistic Boltzmann factor, with the odd harmonics fixed by the pair
symmetry.

*At finite $k$ the quantization survives as a stack of in-band quasimodes.*
The ladder members leave the imaginary axis pairwise as $k$ grows and become
propagating damped modes hugging the light line from below — the "second
family" is simply the least-damped member. At $k = 1.5$, $mu = 2$ carries one
($1.494 - 0.19 i$), while $mu = 10$ carries at least *five*
($1.482 - 0.24 i$, $1.475 - 0.29 i$, $1.434 - 0.45 i$, $1.361 - 0.63 i$,
$1.051 - 1.05 i$; each polished to $|Lambda| < 10^(-9)$): the member count at
fixed $k$ scales with $mu$ like the $k -> 0$ ladder, because a larger $mu$
means more accumulated resonance phase. Each member keeps nearly
$k$-independent damping ($mu = 10$ third member: $gamma = -0.45 plus.minus
0.01$ over $k = 0.6$–$1.4$, descending to its axis pair below $k approx 0.4$;
least-damped member: $-0.24 -> -0.27$ over $k = 1.5$–$4.5$; $mu = 2$:
$approx -0.19$). These are damped quasimodes of the relativistic cyclotron
continuum — discrete interference zeros tied to the band where
$gamma_L (omega - k_parallel v_parallel) = plus.minus Omega$ is solvable —
with no nonrelativistic counterpart. Their large-$k$ fate depends on where
the EM branch sits relative to the band edge: at $mu = 2$,
$2 Pi^2 K_1 slash K_2 = 1.102 approx Omega^2$, so the arriving quasimode
finds a marginal real mode at the edge and lands on it (the $k approx 6.1$
shutoff of §5); at $mu = 10$ the EM branch ($1.716$) is well separated, and
the least-damped quasimode crosses the real-$omega$ band edge at finite depth
($gamma approx -0.28$ near $k approx 5.8$, the edge broadened by its own
damping) and persists as a distinct damped mode (traced to $k = 9$,
$omega_r - k -> 0.16$, $gamma -> -0.31$).

#figure(
  image("fig_quasimodes.png", width: 92%),
  caption: [Representative finite-$k$ members omitted from @fig-replica for
    clarity. Blue: the $mu = 2$ member. Red: two members of the deeper
    $mu = 10$ stack. Dashed line is the light line.],
) <fig-quasimodes>

*Corollary: the non-holomorphic continuation misses the ladder entirely.*
López's $theta$ term freezes the support at $gamma_1 ("Re" z = 0) approx
sqrt(t^2 + 1) approx 1 slash y$, so *their* continuation term is
$O(e^(-mu slash y))$ — invisible. Checked at $k = 0.01$, $mu = 2$: their
$Lambda_L$ is smooth and positive across $s = 4 slash pi$ ($1.47$) while the
correct function swings through zero from $+3 times 10^6$ to $-3 times 10^6$.
The entire aperiodic sector — the roots that carry the published damping dive
— does not exist in their closed form at small $k$; it is recovered only by
the holomorphic continuation.

#v(4mm)
#line(length: 100%, stroke: 0.4pt + gray)
*Reproduction.* From `experiments/lopez-anomalous-zone/`:
`julia --project=../.. 0N_….jl` for `N = 1..7, 9`, and
`julia --project=../../docs 08_make_figures.jl` for report figures.
