#set document(title: [Transverse O#super[+] heating does not explain sub-threshold auroral ion-cyclotron waves — it suppresses them])
#set page(margin: (x: 2.2cm, y: 2.2cm), numbering: "1")
#set text(size: 10.5pt)
#set heading(numbering: "1.1")
#show link: set text(fill: blue.darken(30%))
#set table(stroke: 0.4pt + gray)

#context if target() != "html" { text(16pt, weight: "bold", document.title) }

#v(3mm)
#block(inset: (x: 8mm), stroke: (left: 2pt + gray), outset: (y: 2mm))[
  *Claim.* Electrostatic ion-cyclotron (EIC) waves are seen in the topside
  auroral ionosphere at field-aligned currents below the Kindel & Kennel (1971)
  threshold. Transversely heated O#super[+] is co-observed with them and carries
  perpendicular free energy the current-based threshold ignores, so it is a
  natural candidate to close the gap. It does not. A *minority* heated
  population — the observed situation in the topside — raises the current
  threshold, by 30% at a ring fraction $delta = 0.3$ and by 58% when the whole
  O#super[+] is a ring too broad to be unstable on its own. The currentless
  channel exists but is demanding: it needs a sharp ring
  ($v_r slash w gt.tilde 1.44$, $E_perp gt.tilde 0.63$~eV) that also *dominates*
  the O#super[+] population ($delta gt.tilde 0.6$–$0.9$). The two free-energy
  sources compete rather than add. Separately, and robustly: at identical
  $(n, T_perp, T_parallel)$ a ring and a bi-Maxwellian sit at opposite extremes —
  no current required versus four times the isotropic threshold — so this
  threshold is not a function of the moments that observational studies report.
]

= Why this question

Auroral field-aligned currents flow through a topside ionosphere that is cold
($T tilde 0.3$~eV), dense, and dominated by O#super[+]. Kindel & Kennel
(#link("https://doi.org/10.1029/JA076i013p03055")[JGR 76, 3055], 1971) showed
that such a current destabilises electrostatic ion-cyclotron waves before
anything else, and that result has framed auroral wave physics since: EIC waves
are read as a current diagnostic.

The reading has a long-standing problem. EIC waves are observed when the measured
field-aligned drift sits *below* the critical value. The accepted resolution
(Gavrishchaka et al. 2000; Ganguli et al. 2002; Peñano & Ganguli,
#link("https://doi.org/10.1029/2001JA000279")[JGR 107, 1315], 2002) adds a second
free-energy source: shear in the ion flow changes the sign of the cyclotron
damping term, so EIC waves grow with no current at all.

But the same altitudes host *transverse ion heating*. Ion conics and
transversely accelerated ions (TAI) are among the best-documented features of the
auroral topside (Retterer et al.,
#link("https://doi.org/10.1029/93JA03570")[JGR 99, 13189], 1994). A transversely
heated population is not a displaced Maxwellian: it is a ring or shell in
$v_perp$ with $partial f_0 slash partial v_perp > 0$ on its inner flank — free energy
of exactly the kind cyclotron harmonics feed on. That this works in principle is
established elsewhere in geospace: oxygen cyclotron harmonic waves seen by
Akebono, Van Allen Probes and MMS are attributed to the O#super[+] ion Bernstein
instability driven by ring-like O#super[+] (Min et al.,
#link("https://doi.org/10.1002/2017JA023979")[JGR 122, 2017]; Liu et al.,
#link("https://doi.org/10.1029/2020GL090575")[GRL 47, 2020] and
#link("https://doi.org/10.1029/2022JA030828")[JGR 127, 2022]), and pickup-ion
rings drive the same family at Mars.

So the question is narrow: *does the transversely heated O#super[+] already
present in the auroral topside destabilise the ion-cyclotron band by itself, and
does it lower the current the band would otherwise need?* If it does, the
sub-threshold observations need no shear.

The reason this has not been settled is decisive: the free energy lives in the
shape of $f_0$, and the dispersion solvers in routine use accept only
bi-Maxwellian or kappa distributions — a ring cannot be written in their inputs.
Answering it needs a hot magnetized kinetic solver taking an arbitrary gyrotropic
$f_0$; here VlasovMaxwellDispersion.jl (VMD).

= Setup

*Reference plasma.* Topside auroral, $B_0 = 2.8 times 10^(-5)$~T,
$n = 3 times 10^3$~cm#super[-3], O#super[+] with $T_i = 0.3$~eV and
$T_e slash T_i = 1$ unless stated. Then $Omega_(O^+) slash 2 pi = 26.7$~Hz,
$rho_(O^+) = 11.3$~m, $w_i = 1.9$~km/s,
$(omega_(p i) slash Omega_i)^2 = 1.16 times 10^4$. Ion-neutral collisions are
negligible above $tilde 250$~km ($nu_(i n) slash Omega_i lt.tilde 10^(-2)$).

*Model.* An O#super[+] core (fraction $1 - delta$, isotropic Maxwellian), an
O#super[+] ring (fraction $delta$, gyro-ring $I_0$ Maxwellian at ring speed $v_r$
and the same thermal width $w$), and Maxwellian electrons drifting at $u$ along
$bold(B)$. The knobs are the *ring sharpness* $v_r slash w$ — equivalently
$E_perp = T_i (v_r slash w)^2$ — the *ring fraction* $delta$, and the drift $u$,
reported in electron thermal speeds and as $j = n e u v_e$.

*Solve.* Electrostatic reduction (`mode = :P`), which is the Kindel & Kennel
dispersion relation and is accurate at these $beta tilde 10^(-4)$. Each parameter
point runs seedless box surveys over $k_perp rho in [0.5, 3.5]$,
$k_parallel slash k_perp in [10^(-4), 0.12]$, taking the peak growth; thresholds
are bisected to $gamma = 10^(-3) Omega_(O^+)$. Every growing root in §4 was
re-solved with a seeded Muller polish and returned to the same place.

= Gate: the current-driven threshold

The survey finds the unstable root just above the first gyroharmonic at
$k_perp rho_i approx 1.5$ — textbook EIC — with a threshold falling as
$T_e slash T_i$ rises and ion Landau damping weakens:

#figure(
  table(
    columns: 7,
    align: (right,) * 7,
    [$T_e slash T_i$], [$u_c slash v_e$], [$u_c slash w_i$], [$j_c$ [µA/m#super[2]]], [$k_perp rho$], [$k_parallel slash k_perp$], [$omega_r slash Omega$],
    [0.5], [0.208], [25.2], [23.0], [1.50], [0.035], [1.111],
    [1.0], [0.106], [18.1], [16.5], [1.50], [0.050], [1.186],
    [2.0], [0.059], [14.2], [12.9], [1.25], [0.080], [1.293],
    [4.0], [0.034], [11.7], [10.7], [1.25], [0.120], [1.433],
  ),
  caption: [Current-driven EIC threshold, pure Maxwellian O#super[+].
    $u_c slash w_i approx 12$–$25$ reproduces Kindel & Kennel.],
)

The threshold drift is a property of the temperature ratio, not the density:
$u_c slash v_e = 0.1055$ to four digits across a hundredfold change in $n$, since
$k lambda_(D e) tilde 10^(-2)$ leaves the Debye term negligible. Density enters
the *current* only through $j = n e u_c v_e$.

This is the quantitative form of the puzzle: the classical threshold is
#box[$16.5$~µA/m#super[2]] here, while auroral field-aligned currents are
typically #box[$1$–$10$~µA/m#super[2]]. The waves are observed anyway.

= The ring on its own does work — if it is sharp and it dominates

Set $u = 0$: no current anywhere. Sharpening the ring destabilises the band at

#align(center)[$v_r slash w = 1.44$, i.e. $E_perp = 0.63$~eV, at $omega_r = 0.99 Omega_(O^+)$]

and growth then rises steeply — $gamma = 0.04 Omega$ by $E_perp = 1.2$~eV and
$0.19 Omega$ by $1.9$~eV, e-folding within a few gyroperiods. The unstable
frequency hops between gyroharmonics as the ring sharpens
($omega_r slash Omega = 0.99 → 2.00 → 0.51 → 4.06$), the multi-band signature
characteristic of ring-driven Bernstein modes.

That threshold energy is undemanding — well below observed conic energies. The
binding constraint is elsewhere. Requiring instability at fixed $v_r$ and solving
instead for the *fraction* of O#super[+] that must be ring-like:

#figure(
  table(
    columns: 3,
    align: (right,) * 3,
    [$v_r slash w$], [$E_perp$ [eV]], [$delta_c$],
    [2.0], [1.2], [0.91],
    [2.5], [1.9], [0.82],
    [3.0], [2.7], [0.59],
  ),
  caption: [Ring fraction needed for currentless instability. The heated
    population must be most of the O#super[+], not a suprathermal minority.],
)

= The mechanisms compete, they do not add

The observationally relevant configuration is a heated minority over a dense cold
core. There, transverse heating moves the current threshold the *wrong way*:

#figure(
  table(
    columns: 4,
    align: (right,) * 4,
    [$E_perp$ [eV]], [$j_c$, $delta = 0.1$], [$j_c$, $delta = 0.3$], [$j_c$, $delta = 1$],
    [0.00], [16.5], [16.5], [16.5],
    [0.075], [16.6], [16.9], [18.9],
    [0.30], [16.9], [18.9], [26.1],
    [0.67], [17.4], [20.4], [*0*],
    [1.20], [17.5], [21.5], [*0*],
    [1.88], [17.4], [20.1], [*0*],
  ),
  caption: [Threshold current [µA/m#super[2]] versus how heated the O#super[+]
    already is. A minority ring *raises* it; a dominant ring raises it further
    still — until it crosses its own threshold, where the requirement vanishes.],
)

The $delta = 1$ column is the sharpest statement. Below its own currentless
threshold a fully ring-like O#super[+] makes the current-driven mode *harder* to
excite by 58% ($16.5 → 26.1$~µA/m#super[2] at $E_perp = 0.30$~eV); just past it,
the requirement drops discontinuously to zero. There is no regime in which
transverse heating gently assists the current.

Physically this is not surprising once seen. The current-driven mode at
$omega approx 1.19 Omega$, $k_perp rho approx 1.5$ survives because the ion
perpendicular distribution is narrow enough to keep cyclotron damping small.
A hot ring component broadens $⟨v_perp^2⟩$, pushing the ring ions to larger
$k_perp rho$ and *increasing* the damping of precisely that mode, while the
ring's own unstable modes live at different $(k, omega)$ and need the ring to
dominate before they overcome the cold core's damping. The two channels never
share a mode: they compete for the same population, and below the ring's own
threshold the heating is pure damping.

#figure(
  image("fig_collapse.png", width: 100%),
  caption: [*Left:* currentless growth versus O#super[+] perpendicular energy
    ($delta = 1$); the dashed line is the $gamma = 10^(-3) Omega$ threshold.
    *Right:* the current the same band needs, normalised to the unheated value.
    Every curve starts by moving into the red — heating makes the current-driven
    mode harder — and only $delta = 1$ escapes, by crossing into its own
    currentless instability.],
)

= The moments do not carry the information

A gyro-ring at sharpness $v_r slash w$ has
$⟨v_perp^2⟩ = w^2 + v_r^2$, hence an apparent anisotropy
$T_perp slash T_parallel = 1 + (v_r slash w)^2$. Pair each ring with the
bi-Maxwellian carrying *identical* $n$, $T_perp$, $T_parallel$:

#figure(
  table(
    columns: 5,
    align: (right,) * 5,
    [$T_perp slash T_parallel$], [ring $gamma slash Omega$], [ring $u_c slash v_e$], [bi-Max $gamma slash Omega$], [bi-Max $u_c slash v_e$],
    [1.00], [0], [0.1055], [0], [0.1055],
    [1.25], [0], [0.1211], [0], [0.1201],
    [2.00], [0], [0.1670], [0], [0.1582],
    [3.25], [+0.0055], [*0*], [0], [0.2197],
    [5.00], [+0.0379], [*0*], [0], [0.2939],
    [7.25], [+0.1937], [*0*], [0], [0.4053],
  ),
  caption: [Same first two velocity moments, opposite answers. At
    $T_perp slash T_parallel = 7.25$ the ring needs no current at all; the
    bi-Maxwellian needs four times the isotropic threshold.],
)

Note the bi-Maxwellian trend: raising $T_perp$ at fixed moments-of-record makes
the current-driven mode *harder* to excite, monotonically, because perpendicular
broadening is pure cyclotron damping when there is no positive slope to go with
it. The ring reaches the opposite extreme. Anywhere between "no current needed"
and "4× the isotropic threshold" is compatible with the same reported
$(n, T_perp, T_parallel)$.

#figure(
  image("fig_moments.png", width: 100%),
  caption: [Ring versus moment-matched bi-Maxwellian. Every pair shares $n$,
    $T_perp$, $T_parallel$ exactly.],
)

= What this means

*Transverse heating is not the explanation for sub-threshold EIC.* We set out to
test whether the co-observed conic removes the need for the shear mechanism. It
does not: in the regime the topside actually presents — a suprathermal heated
minority over a cold core — it makes the discrepancy worse. The shear explanation
is not displaced by this alternative.

*Where the ring channel does operate, the population is right.* The requirement
$delta gt.tilde 0.6$ is a strong statement about where to look, and it is
consistent with where ring-driven O#super[+] cyclotron harmonics are actually
reported — the inner magnetosphere and plasma sheet boundary layer, where the
energetic O#super[+] is not sitting on top of a dense cold ionospheric core. The
topside ionosphere is the one place the mechanism is *suppressed*, and it is
suppressed for a reason that is measurable: the cold core fraction.

*The prediction that survives is a distribution measurement.* Where the
O#super[+] perpendicular distribution is a ring with positive inner slope,
$E_perp gt.tilde 0.6$~eV, *and* the heated component dominates the local
O#super[+], the ion-cyclotron band should be unstable irrespective of the
current, with power hopping between harmonics as the ring sharpens. Where the
heated component is a minority, EIC should require *more* current than the
classical threshold, not less — a signed, falsifiable prediction.

*Threshold curves drawn in moment space are not thresholds.* Any stability
criterion for ion outflow expressed in $(n, T_perp, T_parallel)$ — as
parameterisations of auroral heating and outflow generally are — is
under-determined by the range shown in §6. This is the same lesson the solar-wind
anisotropy constraint teaches, arriving in the ionosphere.

*A caution on self-consistency.* The threshold is on ring *sharpness*
$v_r slash w$, not energy: what destabilises the band is the positive
perpendicular slope, and a heating process that broadens the shell as fast as it
displaces it never reaches threshold at any energy. The energies quoted assume
the heated population keeps roughly the core thermal width. Whether real TAI
distributions do is measurable and is the sharpest way to falsify this.

*Scope.* Local, homogeneous, electrostatic, collisionless, single ion species,
linear. The companion note `experiments/shear-kinetic/theory.md` treats the
inhomogeneous (sheared) problem and shows its drive lives entirely in the
flow-advection term — so the two candidate mechanisms are physically distinct and
do not overlap.

= Reproducing

`01_gate_current.jl` → `02_ring_threshold.jl` → `03_joint_threshold.jl` →
`04_moment_control.jl` → `05_make_figures.jl`. Shared model and drivers in
`common.jl`; each script writes its own `out_*.tsv`.
