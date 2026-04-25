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
# extrapolation.jl — Port of DEXTHR
#
# Linear extrapolation for the singular-region U-series to produce improved
# integral estimates and error bounds.
#
# References:
#   Espelid & Genz, Numerical Algorithms 8 (1994) 201-220
#   Espelid, BIT 34 (1994) 62-79

# ============================================================
# DEXTHR — extrapolation tableau update
#
# Index conventions (all 1-based Julia vs 0-based Fortran):
#   t[j, i+1]      = Fortran T(j, i)      for i = 0..emax
#   ne[i]           = Fortran NE(i)         for i = 1..emax
#   beta[i+1, j+1]  = Fortran BETA(i, j)   for i,j = 0..emax
#   uerr[j, k]      = Fortran UERR(j, k)   for k = 1..n_terms  (= E matrix)
#
# Arguments:
#   numfun   — number of integrand components
#   alpha    — degree of homogeneous singularity
#   logf     — 1 if logarithmic singularity, 0 otherwise
#   singul   — dimension of the singularity
#   exstep   — exponent sequence step size (always 1 from the driver)
#   n_terms  — number of U-terms accumulated so far (= NUMU)
#   t        — extrapolation tableau, size (numfun, emax+1); modified in-place
#   update   — 0 = new extrapolation step; >0 = correction of term `update`
#   unew     — new U-term or correction, length numfun
#   qnew     — new tail-correction (in), error accumulator (out), length numfun
#   qold     — previous tail-correction, length numfun
#   first_call — Ref{Bool}: true on first invocation, set to false by caller
#   emax     — maximum number of extrapolation steps
#   uerr     — error estimates for U-terms, size (numfun, n_terms) = E matrix
#   result   — current integral estimates, length numfun; updated if improved
#   abserr   — current error estimates, length numfun; updated if improved
#   exterr   — pure extrapolation error, length numfun (output)
#   ne       — denominator table, length emax; filled on first call
#   beta     — coefficient table, size (emax+1, emax+1); filled on first call
# ============================================================
function _extrapolate!(
        numfun::Int, alpha::Float64, logf::Int, singul::Int,
        exstep::Int, n_terms::Int,
        t::AbstractMatrix, update::Int,
        unew::AbstractVector, qnew::AbstractVector,
        qold::AbstractVector,
        first_call::Ref{Bool}, emax::Int,
        uerr::AbstractMatrix,
        result::AbstractVector, abserr::AbstractVector,
        exterr::AbstractVector,
        ne::Vector{Float64}, beta::Matrix{Float64}
    )

    const_val = 10.0  # heuristic constant for error estimation (CONST in DEXTHR)

    # -------------------------------------------------------
    # One-time initialisation (FIRST = .TRUE. in Fortran)
    # -------------------------------------------------------
    if first_call[]
        # T(J,0) = QOLD(J)  →  t[j,1] = qold[j]
        for j in 1:numfun
            t[j, 1] = qold[j]
        end

        # NE(1) = 2^(SINGUL+ALPHA) - 1
        ne[1] = 2.0^(singul + alpha) - 1.0

        if logf == 1
            # Odd indices: apply exstep doublings of 2*x+1 from NE(I-2)
            for i in 3:2:emax
                es = ne[i - 2]
                for _ in 1:exstep
                    es = 2.0 * es + 1.0
                end
                ne[i] = es
            end
            # Even indices: copy from previous odd
            for i in 2:2:emax
                ne[i] = ne[i - 1]
            end
        else
            # No log singularity: all consecutive
            for i in 2:emax
                es = ne[i - 1]
                for _ in 1:exstep
                    es = 2.0 * es + 1.0
                end
                ne[i] = es
            end
        end

        # Initialise beta coefficients
        # BETA(0,J)=1; BETA(I,J)=BETA(I,J-1)+(BETA(I,J-1)-BETA(I-1,J-1))/NE(J)
        # for I=1..J; BETA(I,J)=0 for I=J+1..EMAX
        for j_f in 0:emax           # Fortran J = 0..EMAX
            beta[1, j_f + 1] = 1.0   # BETA(0,J) = 1
            for i_f in 1:j_f
                beta[i_f + 1, j_f + 1] = beta[i_f + 1, j_f] +
                    (beta[i_f + 1, j_f] - beta[i_f, j_f]) / ne[j_f]
            end
            for i_f in (j_f + 1):emax
                beta[i_f + 1, j_f + 1] = 0.0
            end
        end
    end  # first_call

    steps = min(n_terms, emax)

    # -------------------------------------------------------
    # Update the extrapolation tableau
    # -------------------------------------------------------
    if update == 0
        # New extrapolation step (Fortran UPDATE=0)
        for j in 1:numfun
            save1 = t[j, 1] + (unew[j] + (qnew[j] - qold[j]))
            for i in 1:steps
                # SAVE2 = SAVE1 + (SAVE1 - T(J,I-1)) / NE(I)
                # T[j, i]   = Fortran T(J, I-1)
                save2 = save1 + (save1 - t[j, i]) / ne[i]
                t[j, i] = save1
                save1 = save2
            end
            t[j, steps + 1] = save1   # T(J, STEPS)
        end

    elseif update < n_terms - steps
        # Simple correction: add UNEW to all tableau entries
        for j in 1:numfun
            for i in 1:(steps + 1)   # i = Fortran I+1 for I=0..STEPS
                t[j, i] += unew[j]
            end
        end

    else
        # Weighted correction: T(J,I) += UNEW(J)*(1 - BETA(N-UPDATE+1, I))
        for j in 1:numfun
            for i_f in 0:steps
                # Fortran first-dim index: N-UPDATE+1 (0-based)
                # Clamp to valid range [0, EMAX] to avoid OOB
                fbi = min(n_terms - update + 1, emax)
                beta_val = beta[fbi + 1, i_f + 1]   # BETA(fbi, i_f)
                t[j, i_f + 1] += unew[j] * (1.0 - beta_val)
            end
        end
    end

    # -------------------------------------------------------
    # Error estimates
    # -------------------------------------------------------
    for j in 1:numfun
        exterr[j] = const_val * abs(t[j, steps + 1] - t[j, steps])
        qnew[j] = exterr[j]

        # Effect of U-term errors amplified by extrapolation
        for i in 1:steps
            # ABS(1 - BETA(I, STEPS)) * UERR(J, N+1-I)
            qnew[j] += abs(1.0 - beta[i + 1, steps + 1]) * uerr[j, n_terms + 1 - i]
        end
        # Remaining U-terms (beyond the extrapolation steps)
        for i in (steps + 1):n_terms
            qnew[j] += uerr[j, n_terms + 1 - i]
        end
    end

    # -------------------------------------------------------
    # Update result and abserr only when the new estimate is better
    # -------------------------------------------------------
    for j in 1:numfun
        if qnew[j] <= abserr[j]
            result[j] = t[j, steps + 1]
            abserr[j] = qnew[j]
        end
    end

    return nothing
end
