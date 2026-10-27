# =============================================================================
# nca_lite.jl  — minimal, transparent NCA for the PopPK-predicted profiles.
#
# Computes the three BE endpoints the PBPK side reports, with matching names:
#   Cmax     — observed peak concentration
#   AUC_tEnd — AUC from 0 to last time (linear-up / log-down trapezoid)
#   AUC_inf  — AUC_tEnd + Clast/lambda_z  (terminal slope from log-linear fit)
#
# Designed to run on a DENSE predicted profile (e.g. 0:0.05:tmax h) so the
# trapezoid is effectively exact and lambda_z is well determined. If you would
# rather use the licensed Pumas NCA engine, swap `nca_metrics` for a
# `Pumas.NCA` call returning (Cmax, AUC_tEnd, AUC_inf) — the downstream BE/VBE
# code only needs those three numbers per subject.
#
# `include`d (not a module) after `using Statistics`.
# =============================================================================

"""
    auc_trap(t, c) -> Float64

AUC by the linear-up / log-down trapezoidal rule over the supplied grid.
Falls back to the linear rule on any interval where the log step is undefined
(c <= 0, or a rising/equal segment).
"""
function auc_trap(t::AbstractVector, c::AbstractVector)
    a = 0.0
    @inbounds for i in 1:length(t)-1
        dt = t[i+1] - t[i]
        c1 = c[i]; c2 = c[i+1]
        if c1 > 0 && c2 > 0 && c2 < c1          # falling -> log trapezoid
            a += dt * (c1 - c2) / (log(c1) - log(c2))
        else                                     # rising / flat / zero -> linear
            a += dt * (c1 + c2) / 2
        end
    end
    a
end

"""
    lambda_z(t, c; nmin = 3, nmax = 8) -> (lz, clast, tlast, npts, r2)

Terminal elimination rate from log-linear regression on the last points of a
profile. Scans candidate windows of the final `nmin..nmax` points (after the
peak) and keeps the one with the best adjusted R^2, mirroring standard NCA
practice. Returns `lz = NaN` if no valid descending tail is found.
"""
function lambda_z(t::AbstractVector, c::AbstractVector; nmin = 3, nmax = 8)
    imax = argmax(c)
    # candidate terminal points: strictly after the peak, positive concentration
    idx = [i for i in (imax+1):length(t) if c[i] > 0]
    length(idx) < nmin && return (NaN, c[end], t[end], 0, NaN)
    best = (lz = NaN, clast = c[idx[end]], tlast = t[idx[end]], npts = 0, r2adj = -Inf)
    for k in nmin:min(nmax, length(idx))
        sel = idx[end-k+1:end]
        x = t[sel]; y = log.(c[sel])
        n = length(x)
        mx = sum(x)/n; my = sum(y)/n
        sxx = sum((x .- mx).^2); sxy = sum((x .- mx).*(y .- my))
        sxx == 0 && continue
        slope = sxy / sxx
        slope >= 0 && continue                    # need decline
        intercept = my - slope*mx
        ss_res = sum((y .- (intercept .+ slope.*x)).^2)
        ss_tot = sum((y .- my).^2)
        ss_tot == 0 && continue
        r2 = 1 - ss_res/ss_tot
        r2adj = 1 - (1 - r2)*(n - 1)/(n - 2)
        if r2adj > best.r2adj
            best = (lz = -slope, clast = c[sel[end]], tlast = x[end], npts = n, r2adj = r2adj)
        end
    end
    (best.lz, best.clast, best.tlast, best.npts, best.r2adj)
end

"""
    nca_metrics(t, c) -> NamedTuple

Cmax, Tmax, AUC_tEnd, AUC_inf, and the percent extrapolated, from a dense
concentration-time profile (`t` in hours, `c` in µmol/L).
"""
function nca_metrics(t::AbstractVector, c::AbstractVector)
    cmax, im = findmax(c)
    tmax = t[im]
    auc_t = auc_trap(t, c)
    lz, clast, _, npts, r2adj = lambda_z(t, c)
    auc_inf = (isfinite(lz) && lz > 0) ? auc_t + clast/lz : NaN
    pct_extrap = isfinite(auc_inf) ? 100 * (auc_inf - auc_t)/auc_inf : NaN
    (Cmax = cmax, Tmax = tmax, AUC_tEnd = auc_t, AUC_inf = auc_inf,
     lambda_z = lz, lz_npts = npts, lz_r2adj = r2adj, pct_extrap = pct_extrap)
end
