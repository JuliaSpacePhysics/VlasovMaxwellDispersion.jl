module VlasovMaxwellDispersionRootsExt

using VlasovMaxwellDispersion: _seedsolver, _seed_offset, _split_options, ReturnCode
import VlasovMaxwellDispersion: findroot
import PlasmaDispersionBase: seedsolver
import Roots

seedsolver(prob, alg::Roots.AbstractUnivariateZeroMethod; kw...) = _seedsolver(prob, alg, _split_options(kw)...)

# find_zero terminates on step size, so it both throws on the iteration cap AND can return a finite step-converged non-root.
# Guard with an explicit residual check so a finite ω genuinely means |f(ω)|≈0
function findroot(alg::Roots.AbstractUnivariateZeroMethod, D, ω0; atol = 1.0e-10, maxevals = 100, kw...)
    nevals = Ref(0)
    f = ω -> (nevals[] += 1; D(ω))
    ω = try
        z = ComplexF64(Roots.find_zero(f, ω0 + _seed_offset(ω0) * im, alg; atol, maxevals, kw...))
        abs(f(z)) <= sqrt(atol) * max(abs(f(ω0)), 1) ? z : ComplexF64(NaN, NaN)
    catch
        ComplexF64(NaN, NaN)
    end
    return ω, nevals[], isfinite(ω) ? ReturnCode.Success : ReturnCode.Failure
end

end
