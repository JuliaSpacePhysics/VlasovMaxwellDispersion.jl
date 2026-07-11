# López et al. 2014 (PoP 21,092107) L-mode dispersion Λ_L, transcribed from
# Eqs (17,18,21-26) + Appendix A (A7-A20). Normalized: x=ω/Ωc, y=ck/Ωc,
# z=x/y=ω/(ck), t=Ωc/(ck)=1/y. Electron uses t_e=-t (A8/A9 sign), positron t_p=+t.
using QuadGK, SpecialFunctions

# --- S limits (A8,A9 electron; A15,A16 positron), t>0 ---
S1e(γ, t) = sqrt(1 - 1/γ^2) + t/γ
S2e(γ, t) = sqrt(1 - 1/γ^2) - t/γ
S1p(γ, t) = sqrt(1 - 1/γ^2) - t/γ
S2p(γ, t) = sqrt(1 - 1/γ^2) + t/γ

σofI(I) = I < 0 ? 2.0 : (I == 0 ? 1.0 : 0.0)
Θ(x) = x > 0 ? 1.0 : 0.0

# electron Heaviside boundaries (A13)
γ1e(R, t) = (R*t + sqrt(complex(t^2 + 1 - R^2))) / (1 - R^2) |> real
γ2e(R, t) = (R*t - sqrt(complex(t^2 + 1 - R^2))) / (1 - R^2) |> real
# positron (A20)
γ1p(R, t) = (-R*t + sqrt(complex(t^2 + 1 - R^2))) / (1 - R^2) |> real
γ2p(R, t) = (-R*t - sqrt(complex(t^2 + 1 - R^2))) / (1 - R^2) |> real

function θe(γ, R, I, t)
    σ = σofI(I)
    if R <= -sqrt(1 + t^2)
        return 0.0
    elseif R < -1
        return π*σ*Θ(γ - γ1e(R,t))*Θ(γ2e(R,t) - γ)
    elseif R < 1
        return π*σ*Θ(γ - γ1e(R,t))
    else
        return 0.0
    end
end

function θp(γ, R, I, t)
    σ = σofI(I)
    if R <= -1
        return 0.0
    elseif R <= 1
        return π*σ*Θ(γ - γ1p(R,t))
    elseif R < sqrt(1 + t^2)
        return π*σ*Θ(γ - γ1p(R,t))*Θ(γ2p(R,t) - γ)
    else
        return 0.0
    end
end

# closed-form continuation (A7 electron / A20 positron)
function Je(γ, z, t)
    R, I = real(z), imag(z)
    s1, s2 = S1e(γ,t), S2e(γ,t)
    re = 0.5*log(((R-s2)^2 + I^2)/((R+s1)^2 + I^2))
    im_ = atan((s2-R)/I) + atan((s1+R)/I) + θe(γ,R,I,t)
    complex(re, im_)
end
function Jp(γ, z, t)
    R, I = real(z), imag(z)
    s1, s2 = S1p(γ,t), S2p(γ,t)
    re = 0.5*log(((R-s2)^2 + I^2)/((R+s1)^2 + I^2))
    im_ = atan((s2-R)/I) + atan((s1+R)/I) + θp(γ,R,I,t)
    complex(re, im_)
end

# direct un-continued integral (Im z>0): J = ∫_{-S1}^{S2} dξ/(ξ-z)=Log((S2-z)/(-S1-z))
Je_direct(γ, z, t) = log((S2e(γ,t) - z)/(-S1e(γ,t) - z))
Jp_direct(γ, z, t) = log((S2p(γ,t) - z)/(-S1p(γ,t) - z))

# --- Λ_L (Eq 26), ωpe²/Ωc²=1 ---
function ΛL(x, y, μ; Jefun=Je, Jpfun=Jp)
    z = x/y; t = 1/y
    K2 = besselk(2, μ)
    term = 1 - y^2/x^2 - μ/y^2
    pref = (μ^2/(4*K2)) / (x*y^3)
    Ie, _ = quadgk(γ -> exp(-μ*γ)*Jefun(γ,z,t)*((y^2-x^2)*γ^2 - 2*x*γ - (1+y^2)), 1, Inf; rtol=1e-9)
    Ip, _ = quadgk(γ -> exp(-μ*γ)*Jpfun(γ,z,t)*((y^2-x^2)*γ^2 + 2*x*γ - (1+y^2)), 1, Inf; rtol=1e-9)
    term + pref*(Ie + Ip)
end

# complex Muller
function muller(f, x0, x1, x2; tol=1e-12, maxit=100)
    f0,f1,f2 = f(x0),f(x1),f(x2)
    for _ in 1:maxit
        q = (x2-x1)/(x1-x0)
        A = q*f2 - q*(1+q)*f1 + q^2*f0
        B = (2q+1)*f2 - (1+q)^2*f1 + q^2*f0
        C = (1+q)*f2
        den1 = B + sqrt(B^2 - 4*A*C); den2 = B - sqrt(B^2 - 4*A*C)
        den = abs(den1) > abs(den2) ? den1 : den2
        x3 = x2 - (x2-x1)*(2C/den)
        x0,x1,x2 = x1,x2,x3
        f0,f1,f2 = f1,f2,f(x3)
        abs(x2-x1) < tol && break
    end
    x2
end
