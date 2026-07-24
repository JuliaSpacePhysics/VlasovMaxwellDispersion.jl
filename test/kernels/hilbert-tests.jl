@testitem "parallel Hilbert primitive (hilbert_pwpoly)" begin
    using VlasovMaxwellDispersion: cell_hilbert, cell_hilbert_landau, hilbert_pwpoly,
                                   hilbert_landau_pwpoly, Z

    # --- self-contained high-order Gauss-Legendre (no QuadGK dep) ---
    # Integrand p(v)/(v-ζ) is analytic on a cell when Im ζ ≠ 0 ⇒ GL converges
    # spectrally. Newton on Legendre P_n for nodes/weights.
    function gauss_legendre(n)
        x = zeros(n); w = zeros(n)
        for i in 1:n
            z = cos(pi * (i - 0.25) / (n + 0.5))
            local dp
            for _ in 1:100
                p0 = 1.0; p1 = 0.0
                for j in 1:n
                    p2 = p1; p1 = p0
                    p0 = ((2j - 1) * z * p1 - (j - 1) * p2) / j
                end
                dp = n * (z * p0 - p1) / (z^2 - 1)
                dz = p0 / dp; z -= dz
                abs(dz) < 1e-15 && break
            end
            x[i] = z
            w[i] = 2 / ((1 - z^2) * dp^2)
        end
        x, w
    end
    GLX, GLW = gauss_legendre(64)

    polyval(c, v) = (acc = zero(complex(promote_type(eltype(c), typeof(v)))); for k in length(c):-1:1; acc = c[k] + v * acc; end; acc)
    polyr(c, x) = (acc = 0.0; for k in length(c):-1:1; acc = real(c[k]) + x * acc; end; acc)

    # numerical ∫_{vl}^{vr} p(v)/(v-ζ) dv on [-1,1] map (ζ complex, off-axis)
    function cell_num(c, vl, vr, ζ)
        h = (vr - vl) / 2; m = (vr + vl) / 2
        s = zero(complex(float(ζ)))
        for (xi, wi) in zip(GLX, GLW)
            v = m + h * xi
            s += wi * polyval(c, v) / (v - ζ)
        end
        h * s
    end
    pwpoly_num(coeffs, nodes, ζ) = sum(cell_num(coeffs[i], nodes[i], nodes[i+1], ζ) for i in eachindex(coeffs))

    # principal value of ∫ p(v)/(v-x) dv (x real interior) via singularity
    # subtraction: ∫[p(v)-p(x)]/(v-x) dv (smooth) + p(x)·log|(vr-x)/(vl-x)|.
    function cell_pv(c, vl, vr, x)
        h = (vr - vl) / 2; m = (vr + vl) / 2
        s = 0.0; px = polyr(c, x)
        for (xi, wi) in zip(GLX, GLW)
            v = m + h * xi
            s += v == x ? 0.0 : wi * (polyr(c, v) - px) / (v - x)
        end
        h * s + px * log(abs((vr - x) / (vl - x)))
    end

    # --- sample piecewise cubics on a 3-cell grid ---
    nodes = [-2.0, -0.5, 0.7, 1.8]
    coeffs = [
        [0.3, 1.0, -0.4, 0.2],    # 0.3 + v - 0.4v² + 0.2v³
        [1.1, -0.7, 0.9, -0.3],
        [-0.2, 0.5, 0.6, 0.15],
    ]
    cellof(x) = findlast(<=(x), nodes[1:end-1])
    g(x) = polyr(coeffs[cellof(x)], x)
    pvtot(x) = sum(cell_pv(coeffs[i], nodes[i], nodes[i+1], x) for i in eachindex(coeffs))

    # ---- (a) exactness vs numerical Cauchy integral, several complex ζ ----
    ζs = [1.5 + 0.8im, -1.0 + 0.3im, 0.2 - 0.6im, 3.0 + 2.0im,
          -3.5 - 1.0im, 0.0 + 0.25im, 0.5 - 1.5im]
    @testset "exactness vs numerical" begin
        for ζ in ζs
            exact = hilbert_pwpoly(coeffs, nodes, ζ)
            num = pwpoly_num(coeffs, nodes, ζ)
            @test isapprox(exact, num; rtol=1e-9, atol=1e-10)
        end
    end

    @testset "single cell" begin
        c = [2.0, -1.0, 0.5, 0.1]
        for ζ in (0.9 + 0.4im, -1.3 - 0.7im)
            @test isapprox(cell_hilbert(c, -2.0, 1.8, ζ), cell_num(c, -2.0, 1.8, ζ); rtol=1e-9)
        end
    end

    # ---- (b) branch-cut continuity & Plemelj / Landau-causal limit ----
    # Plemelj for H(ζ)=∫g/(v-ζ): pole at v=ζ, residue +2πi g as ζ enters the
    # support from above ⇒  H(x+i0) - H(x-i0) = 2πi g(x)  (jump across cut), and
    # the causal (Landau, Im ζ→0⁺) split  H(x+i0) = PV ∫g/(v-x) dv + iπ g(x).
    @testset "Plemelj jump across axis" begin
        for x in (-1.0, 0.2, 1.0)        # interior, off cell boundaries
            ε = 1e-7
            above = hilbert_pwpoly(coeffs, nodes, x + ε * im)
            below = hilbert_pwpoly(coeffs, nodes, x - ε * im)
            @test isapprox(above - below, 2π * im * g(x); rtol=1e-5, atol=1e-6)
        end
    end

    @testset "Landau-causal Plemelj decomposition (Im ζ→0⁺)" begin
        for x in (-1.0, 0.2, 1.0)
            ε = 1e-8
            above = hilbert_pwpoly(coeffs, nodes, x + ε * im)
            expected = pvtot(x) + im * π * g(x)
            @test isapprox(above, expected; rtol=1e-5, atol=1e-6)
        end
    end

    # ---- continuity-error report across the axis ----
    # Off-support the value is single-valued (jump→0); on-support the jump equals
    # the Plemelj residue. Report (i) max residual of the jump identity, (ii) max
    # |Im of (above+below)/2| which must vanish since the average is the real PV.
    @testset "continuity error report" begin
        maxjumperr = 0.0
        maxpverr = 0.0
        for x in range(-1.9, 1.7; length=51)
            any(n -> isapprox(x, n; atol=1e-9), nodes) && continue
            ε = 1e-7
            above = hilbert_pwpoly(coeffs, nodes, x + ε * im)
            below = hilbert_pwpoly(coeffs, nodes, x - ε * im)
            maxjumperr = max(maxjumperr, abs((above - below) - 2π * im * g(x)))
            maxpverr = max(maxpverr, abs(imag((above + below) / 2)))
        end
        @info "Hilbert branch-cut continuity across Im ζ=0" maxjumperr maxpverr
        @test maxjumperr < 1e-4
        @test maxpverr < 1e-4
    end

    @testset "MPDES Landau continuation" begin
        for x in (-1.0, 0.2, 1.0)
            ε = 1e-7
            above = hilbert_pwpoly(coeffs, nodes, x + ε * im)
            below_causal = hilbert_landau_pwpoly(coeffs, nodes, x - ε * im)
            @test isapprox(below_causal, above; rtol=1e-5, atol=1e-6)
        end

        c = [1.2, -0.4, 0.3]
        x = 0.1
        ε = 1e-7
        @test isapprox(cell_hilbert_landau(c, -1.0, 2.0, x - ε * im),
                       cell_hilbert(c, -1.0, 2.0, x + ε * im);
                       rtol=1e-5, atol=1e-6)
    end

    @testset "MPDES polynomial Maxwellian matches fast Z" begin
        nodes_m = collect(range(-12.0, 12.0; length=1201))
        coeffs_m = map(1:length(nodes_m)-1) do i
            vl, vr = nodes_m[i], nodes_m[i + 1]
            fl = exp(-vl^2) / sqrt(pi)
            fr = exp(-vr^2) / sqrt(pi)
            slope = (fr - fl) / (vr - vl)
            [fl - slope * vl, slope]
        end
        for ζ in (0.3 + 0.4im, -1.2 + 0.2im, 0.7 - 0.3im)
            @test isapprox(hilbert_landau_pwpoly(coeffs_m, nodes_m, ζ), Z(ζ);
                           rtol=5e-2, atol=1e-3)
        end
    end
end
