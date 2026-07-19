# Mode reduction

A dispersion solve chases zeros of a scalar built from the 3×3 dispersion
tensor $\mathcal{D}(\omega, 𝐤)$. The default is $\det\mathcal{D}$, 
exact at every $𝐤$. At the two symmetry axes ($𝐤 \parallel \mathbf{B}_0$
and $𝐤 \perp \mathbf{B}_0$), the tensor block-diagonalizes and the determinant factors;
the `mode` keyword (`:det`, `:L`, `:R`, `:P`, `:O`, `:X`) selects the tracked factor.
This page derives the factorizations, their validity conditions, and why factors matter.

## Why factor the determinant

- **Simple vs. double zeros.** Where a symmetry makes two factors degenerate,
  every determinant zero is a double zero: bracketing solvers see no sign change
  and Newton-type solvers lose their quadratic convergence. The canonical case is 
  the equal-mass pair plasma, where charge-conjugation symmetry forces $R \equiv L$:
  transverse root is double in det, while each circular factor crosses zero simply.
- **Mode selection.** A factor tracks one polarization family; the det mixes all.
- Cheaper per evaluation.

## Parallel propagation: `det = L·R·P`

At $𝐤 \parallel 𝐁_0$ the system is rotationally symmetric about $𝐁_0 = B_0\hat{z}$. 
Gyrotropic $f_0(p_\perp, p_\parallel)$
response must commute with every rotation $R_\phi$ about $\hat{z}$,
$R_\phi \mathcal{D} R_\phi^T = \mathcal{D}$, which forces
$\mathcal{D}_{11} = \mathcal{D}_{22}$, $\mathcal{D}_{12} = -\mathcal{D}_{21}$ and the $13/23/31/32$ entries to vanish.

In the circular basis $\hat{e}_\pm = (\hat{x} \pm i\hat{y})/\sqrt{2}$
the tensor is diagonal:

```math
\det\mathcal{D}
 = \underbrace{(\mathcal{D}_{11} + i\mathcal{D}_{12})}_{L}\,
   \underbrace{(\mathcal{D}_{11} - i\mathcal{D}_{12})}_{R}\,
   \underbrace{\mathcal{D}_{33}}_{P}.
```

With the $e^{-i\omega t}$ convention, the $\sigma = +1$ factor `:L` is the
polarization rotating with **positive** charges — its cyclotron resonance is
$\gamma(\omega - k_\parallel v_\parallel) \to + \Omega$. `:R` ($\sigma = -1$)
resonates with negative charges.

The `P` factor $\mathcal{D}_{33} = \hat{k}^T \mathcal{D} \hat{k}$ is the
parallel electrostatic branch: at $k_\perp = 0$ the curl term
$𝐤𝐤^T - k^2 I$ is purely transverse to $𝐤$.

## Oblique propagation

Away from $\hat{z}$ no exact scalar factor exists — the polarization
eigenvectors depend on $\omega$, and eigenvalue branches are non-analytic.
The one useful oblique reduction is the **electrostatic approximation** `:P`,
$\hat{k}^T \mathcal{D}\, \hat{k}$, which keeps only the field component along
$𝐤$; and is a good approximation where $|E_\perp| \ll |E_\parallel|$ 
(e.g. short-wavelength Bernstein-like and Langmuir-like branches).

## Perpendicular propagation: `det = O·X`

At $k_\parallel = 0$ the relevant symmetry is the **mirror** $z \to -z$,
$\Sigma = \mathrm{diag}(1, 1, -1)$. It preserves the whole system:
$\mathbf{B}_0$ is an axial vector, so $B_0\hat{z}$ is invariant;
$\mathbf{k} = k_\perp\hat{x}$ has no $z$-component; and the equilibrium maps
$f_0(p_\perp, p_\parallel) \to f_0(p_\perp, -p_\parallel)$ — invariant iff
$f_0$ is even. The response then obeys $\Sigma \mathcal{D} \Sigma^T = \mathcal{D}$,
and since $(\Sigma \mathcal{D} \Sigma^T)_{13} = -\mathcal{D}_{13}$, the $13/23/31/32$
entries vanish. (At $k_\parallel \neq 0$ the mirror flips $k_\parallel$, so it only 
relates $\mathcal{D}(k_\perp, k_\parallel)$ to $\mathcal{D}(k_\perp, -k_\parallel)$.)
The tensor splits into $\mathcal{D}_{33}$ and the transverse–longitudinal 2×2 block:

```math
\det\mathcal{D}
 = \underbrace{\mathcal{D}_{33}}_{O}\;
   \underbrace{(\mathcal{D}_{11}\mathcal{D}_{22} - \mathcal{D}_{12}\mathcal{D}_{21})}_{X}.
```

`:O` is the ordinary mode ($E \parallel B_0$); `:X` holds the extraordinary 
mode and the Bernstein family. Because the condition is a *plasma* symmetry,
a field-aligned drift breaks the mirror and recouples the blocks.

## API

```@docs; canonical = false
TensorReduction
```
