# Run with: julia --project=. test/test-cold-stix.jl
# Validate the cold (ColdVDF) dispersion against Stix textbook closed forms.
# Two-species e-p plasma, Omega_ref = |Omega_e| (so Omega_e=-1, Omega_i=+1/(mp/me)).

using Test
using VlasovMaxwellDispersion
using VlasovMaxwellDispersion: Species, Wavenumber, 𝒟
using LinearAlgebra
using Roots
const VM = VlasovMaxwellDispersion

# --- shared plasma: electron + proton, arbitrary Pi2_e (wpe/Omega_ref)^2 ---
const mp_me = 1836.15
const Omega_e = -1.0
const Omega_i = 1 / mp_me
const Pi2_e = 4.0
const Pi2_i = Pi2_e / mp_me   # quasineutral equal-density e-p: wpi^2 = wpe^2 * (me/mi)

const electrons = Species(Omega_e, Pi2_e, VM.ColdVDF())
const protons = Species(Omega_i, Pi2_i, VM.ColdVDF())
const plasma = VM.Plasma(electrons, protons)

# Stix S, D, P closed forms (two-fluid cold), summed over species (Ch.1-2).
stix_S(om) = 1 - Pi2_e / (om^2 - Omega_e^2) - Pi2_i / (om^2 - Omega_i^2)
stix_D(om) = Omega_e * Pi2_e / (om * (om^2 - Omega_e^2)) + Omega_i * Pi2_i / (om * (om^2 - Omega_i^2))
stix_P(om) = 1 - Pi2_e / om^2 - Pi2_i / om^2
stix_R(om) = stix_S(om) + stix_D(om)
stix_L(om) = stix_S(om) - stix_D(om)


@testset "dielectric tensor matches S,D,P closed forms" begin
    ω = 0.5 + 0im
    k = Wavenumber(0.2, 0.3)  # oblique; S,D,P are k-independent for cold
    ε = VM.dielectric(plasma, ω, k)
    @test real(ε[1, 1]) ≈ stix_S(ω) rtol = 1.0e-12
    @test real(im * ε[1, 2]) ≈ stix_D(ω) rtol = 1.0e-12   # ε_xy = -iD (Stix convention)
    @test real(ε[3, 3]) ≈ stix_P(ω) rtol = 1.0e-12
    @test ε[1, 1] ≈ ε[2, 2]
    @test ε[1, 2] ≈ -ε[2, 1]
end

@testset "vacuum sanity (cold limit, Pi2=0): n^2=1 light line" begin
    vac = VM.Plasma(Species(0.0, 0.0, VM.ColdVDF()))
    ω = 1.0 + 0im
    k = Wavenumber(0.0, 1.0)  # n = kc/ω = 1
    @test abs(det(𝒟(vac, ω, k))) < 1.0e-12
end

@testset "parallel propagation (k_perp=0): R, L, plasma(O) wave factorization" begin
    # Stix: det(D) for k‖B0 factors as (R-n^2)(L-n^2)*P, n^2=(kz/ω)^2.
    # omega=5.0 chosen so R,L,P>0 (propagating roots exist; see S,D,P scan above).
    ω = 5.0 + 0im
    kz = 0.3
    k = Wavenumber(0.0, kz)
    n2 = (kz / ω)^2
    R, L, P = stix_R(ω), stix_L(ω), stix_P(ω)
    expected = (R - n2) * (L - n2) * P
    @test det(𝒟(plasma, ω, k)) ≈ expected rtol = 1.0e-10

    # R-wave root: det(D) should vanish at n^2=R (set kz to match).
    kz_R = sqrt(real(R)) * real(ω)
    kR = Wavenumber(0.0, kz_R)
    @test abs(det(𝒟(plasma, ω, kR))) < 1.0e-8 * abs(R * L * P)

    # L-wave root similarly at n^2=L.
    kz_L = sqrt(real(L)) * real(ω)
    kL = Wavenumber(0.0, kz_L)
    @test abs(det(𝒟(plasma, ω, kL))) < 1.0e-8 * abs(R * L * P)
end

@testset "perpendicular propagation (k_par=0): O-mode and X-mode factorization" begin
    # Stix: det(D) for k⊥B0 factors as (P-n^2)*(R*L - S*n^2), n^2=(kperp/ω)^2.
    ω = 5.0 + 0im
    kperp = 0.3
    k = Wavenumber(kperp, 0.0)
    n2 = (kperp / ω)^2
    S, R, L, P = stix_S(ω), stix_R(ω), stix_L(ω), stix_P(ω)
    expected = (P - n2) * (R * L - S * n2)
    @test det(𝒟(plasma, ω, k)) ≈ expected rtol = 1.0e-10

    # O-mode: n^2 = P exactly (electrostatic-free, decoupled zz block).
    kperp_O = sqrt(real(P)) * real(ω)
    kO = Wavenumber(kperp_O, 0.0)
    @test abs(det(𝒟(plasma, ω, kO))) < 1.0e-8 * abs(R * L * P)

    # X-mode: n^2 = R*L/S.
    n2_X = real(R * L / S)
    n2_X > 0 || error("X-mode n^2<0 at this test point; pick another omega")
    kperp_X = sqrt(n2_X) * real(ω)
    kX = Wavenumber(kperp_X, 0.0)
    @test abs(det(𝒟(plasma, ω, kX))) < 1.0e-7 * abs(R * L * P)
end

@testset "cold two-fluid cutoffs (n^2=0, k->0)" begin
    # R-cutoff: R(ω)=0. L-cutoff: L(ω)=0. P-cutoff: ω=wp (P=0).
    k0 = Wavenumber(0.0, 1.0e-6)  # k->0 numerically (S,D,P independent of k)
    Rfun(om) = real(stix_R(om))
    Lfun(om) = real(stix_L(om))
    Pfun(om) = real(stix_P(om))

    ω_Rcut = find_zero(Rfun, 2.0)
    ω_Lcut = find_zero(Lfun, 1.5)
    ω_Pcut = find_zero(Pfun, 2.0)  # ≈ wpe for wpe >> wpi

    # det(D) must vanish at the cutoffs in the k->0 limit (n^2->0).
    @test abs(det(𝒟(plasma, complex(ω_Rcut), k0))) < 1.0e-6
    @test abs(det(𝒟(plasma, complex(ω_Lcut), k0))) < 1.0e-6
    @test abs(det(𝒟(plasma, complex(ω_Pcut), k0))) < 1.0e-6

    # P-cutoff should sit near the (total) plasma frequency wp=sqrt(Pi2_e+Pi2_i).
    wp = sqrt(Pi2_e + Pi2_i)
    @test ω_Pcut ≈ wp rtol = 1.0e-3
end

@testset "cyclotron resonances (S,D divergent poles at omega->|Omega_s|)" begin
    # R diverges (pole) as ω -> |Omega_e| from above (electron cyclotron resonance).
    ωe = abs(Omega_e)
    @test abs(stix_R(ωe * 1.001)) > 1.0e3
    @test abs(stix_L(ωe * 1.001)) < 10  # L regular there (no resonance for L at Ωe)

    # L diverges at ω -> |Omega_i| (ion cyclotron resonance, L-wave pole).
    ωi = abs(Omega_i)
    @test abs(stix_L(ωi * 1.001)) > 1.0e2
end

@testset "upper hybrid resonance (S=0, k->0 electrostatic limit)" begin
    k0 = Wavenumber(0.0, 1.0e-6)
    Sfun(om) = real(stix_S(om))
    ω_UHR = find_zero(Sfun, 2.3)
    wp2 = Pi2_e + Pi2_i
    Oe2 = Omega_e^2  # dominant electron contribution
    ω_UHR_approx = sqrt(wp2 + Oe2)
    @test ω_UHR ≈ ω_UHR_approx rtol = 1.0e-2  # ion mass correction ~1/mp_me

    # electrostatic_det (k.eps.k) should also vanish there as kperp->0 along x.
    kx = Wavenumber(1.0e-6, 0.0)
    @test abs(VM.electrostatic_det(plasma, complex(ω_UHR), kx)) < 1.0e-6
end
