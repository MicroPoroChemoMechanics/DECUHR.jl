# ***************************************************************************
# * All the software  contained in this library  is protected by copyright. *
# * Permission  to use, copy, modify, and  distribute this software for any *
# * purpose without fee is hereby granted, provided that this entire notice *
# * is included  in all copies  of any software which is or includes a copy *
# * or modification  of this software  and in all copies  of the supporting *
# * documentation for such software.                                        *
# ***************************************************************************
# * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
# * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
# * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
# * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
# * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
# * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
# ***************************************************************************
# * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
# * ABOVE STATEMENT.                                                        *
# ***************************************************************************
#
# Reference: T.O. Espelid and A. Genz, "DECUHR: An Algorithm for Automatic
# Integration of Singular Functions over a Hyperrectangular Region",
# Numerical Algorithms 8 (1994), pp. 201-220.
#
# alpha_estimation.jl — Port of DECALP
#
# Estimates the singularity exponent ALPHA and detects logarithmic components
# by evaluating the integrand along a line toward the singular vertex and
# extrapolating the ratio of successive function values.
#
# References:
#   Espelid & Genz, Numerical Algorithms 8 (1994) 201-220
#   Bjorstad, Grosse & Dahlquist, BIT 21 (1981)

# ============================================================
# _estimate_alpha! — port of DECALP
#
# Arguments:
#   ndim      — number of dimensions
#   a, b      — lower/upper integration bounds (length ndim)
#   alpha_ref — Ref{Float64}: on exit contains estimated alpha
#   singul    — dimension of the singularity (first `singul` coordinates)
#   logf_ref  — Ref{Int}: on exit 1 if log term detected, 0 otherwise
#   funsub    — (x::Vector{Float64}, funvls::Vector{Float64}) -> nothing
#   x_work    — work vector length ≥ ndim
#   funs_work — work vector length ≥ numfun
# ============================================================
function _estimate_alpha!(
        ndim::Int, a::AbstractVector, b::AbstractVector,
        alpha_ref::Ref{Float64}, singul::Int,
        logf_ref::Ref{Int},
        funsub,
        x_work::Vector{Float64},
        funs_work::Vector{Float64}
    )
    bdg = 10
    l = 10
    lgtwo = log(2.0)
    ntot = 2 * bdg + l   # ratio table length; indices 0..ntot → Julia [1..ntot+1]
    t_arr = zeros(Float64, ntot + 1)   # t_arr[j+1] = Fortran T(j)

    # -------------------------------------------------------
    # Step 1 : build the ratio table T(0..2*BDG+L)
    #
    # Fortran structure (labels 10 / 60):
    #   H = 2; label 10: H = H/2; set ALL coords; call FUNSUB
    #   if sum|f|=0, go to 10 (keep halving)
    #   then for J=0..ntot: H = H/2; set SINGUL coords; call FUNSUB
    #     if sum|f|=0, go to 10 (restart from current H)
    #   T(J) = log(prev/curr) / log2
    # -------------------------------------------------------
    h = 1.0   # Fortran: H=2 then H=H/2 at label 10 → first value = 1
    sumy = 0.0
    table_built = false

    for _restart in 1:500   # outer restart loop (Fortran GOTO 10)

        # Find starting h giving non-zero |f|
        found_start = false
        for _step in 1:200
            for i in 1:ndim
                x_work[i] = a[i] + h * (b[i] - a[i])
            end
            funsub(x_work, funs_work)
            sumy = sum(abs, funs_work)
            if sumy > 0.0
                found_start = true
                break
            end
            h /= 2
        end
        found_start || break   # give up if no nonzero point found

        # Build ratio table (only singul coordinates continue halving)
        inner_failed = false
        for j in 0:ntot
            sumx = sumy
            h /= 2
            for i in 1:singul
                x_work[i] = a[i] + h * (b[i] - a[i])
            end
            funsub(x_work, funs_work)
            sumy = sum(abs, funs_work)
            if sumy <= 0.0
                inner_failed = true
                break   # current h restarts the outer loop
            end
            t_arr[j + 1] = log(sumx / sumy) / lgtwo
        end

        if !inner_failed
            table_built = true
            break
        end
        # h is already halved; outer loop restarts from here
    end

    table_built || return   # could not build a clean ratio table

    # -------------------------------------------------------
    # Step 2 : Aitken-like linear extrapolation
    # After this: t_arr[1..2*BDG+1] are all estimates of ALPHA
    # -------------------------------------------------------
    for i in 1:ntot
        n_val = 0
        for j_val in (i - 1):-1:max(0, i - l)
            n_val = 2 * n_val + 1
            # T(j_val) = T(j_val+1) + (T(j_val+1)-T(j_val)) / n_val
            t_arr[j_val + 1] = t_arr[j_val + 2] +
                (t_arr[j_val + 2] - t_arr[j_val + 1]) / n_val
        end
    end

    # -------------------------------------------------------
    # Step 3 : assign alpha and detect log term
    # T(2*BDG) = t_arr[2*bdg+1];  T(2*BDG-1) = t_arr[2*bdg]
    # -------------------------------------------------------
    if abs(t_arr[2 * bdg + 1] - t_arr[2 * bdg]) > 5.0e-6
        # Logarithmic term present
        logf_ref[] = 1

        # Step 3b : Bjorstad-Grosse-Dahlquist log-elimination
        # NB decreases by 2 from 2*BDG-2 until convergence or NB = 0
        signal = -1
        nb = 2 * bdg   # immediately decremented to 2*bdg-2 in loop
        while true
            nb -= 2
            itw = 2 * bdg - nb   # increments by 2 each step

            # Estimate convergence order KEST
            if nb > 0
                denom1 = t_arr[nb + 2] + t_arr[nb] - 2.0 * t_arr[nb + 1]
                denom2 = t_arr[nb + 3] + t_arr[nb + 1] - 2.0 * t_arr[nb + 2]
                term1 = abs(denom1) > 0.0 ?
                    (t_arr[nb + 2] - t_arr[nb + 1]) / denom1 : 0.0
                term2 = abs(denom2) > 0.0 ?
                    (t_arr[nb + 3] - t_arr[nb + 2]) / denom2 : 0.0
                kest = -1.0 - 1.0 / (term2 - term1)
            else
                kest = 0.0
            end

            # Update tableau entries T(0)..T(NB)
            for j_val in 0:nb
                denom = t_arr[j_val + 3] + t_arr[j_val + 1] - 2.0 * t_arr[j_val + 2]
                if abs(denom) > 0.0
                    t_arr[j_val + 1] = t_arr[j_val + 2] -
                        itw * (t_arr[j_val + 3] - t_arr[j_val + 2]) *
                        (t_arr[j_val + 2] - t_arr[j_val + 1]) /
                        (denom * (itw - 1))
                else
                    t_arr[j_val + 1] = t_arr[j_val + 2]
                end
            end

            # BGD stopping criterion
            if abs(kest - itw + 1) > 1.0 && signal < 0
                signal = max(0, nb - 2)
            end

            (nb == signal || nb <= 0) && break
        end

        index = max(0, signal)   # safety: default to T(0) if never set
        alpha_ref[] = t_arr[index + 1]

    else
        # No log term: use T(L) as the best estimate of ALPHA
        logf_ref[] = 0
        alpha_ref[] = t_arr[l + 1]   # T(L) in Fortran = t_arr[L+1]
    end

    return nothing
end
