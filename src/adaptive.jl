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
# adaptive.jl — Port of DECHHR, DETRHR, DESBHR, DEADHR
#
# Parameter validation, binary max-heap management, sub-region subdivision,
# and the main adaptive integration loop.
#
# Index convention (Fortran 0-based → Julia 1-based):
#   Fortran index 0 (singular region) → Julia index 1
#   Fortran index k (k ≥ 1)          → Julia index k+1
#   list[pos] stores 1-based Julia subregion indices
#
# References:
#   Espelid & Genz, Numerical Algorithms 8 (1994) 201-220

"""
    ifail_message(ifail::Integer) -> String

Human-readable description of a raw DECUHR `IFAIL` code (as exposed in
`sol.stats.ifail`).  Not exported; call as `DECUHR.ifail_message(code)`.
"""
function ifail_message(ifail::Integer)
    ifail == 0 && return "normal exit: requested accuracy reached"
    ifail == 1 && return "maxiters (MAXPTS) exhausted before the requested accuracy was reached"
    ifail == 2 && return "KEY < 0 or KEY > 4"
    ifail == 3 && return "NDIM < 2 or NDIM > $_MAXDIM"
    ifail == 4 && return "KEY == 1 requires NDIM == 2"
    ifail == 5 && return "KEY == 2 requires NDIM == 3"
    ifail == 6 && return "NUMFUN < 1"
    ifail == 7 && return "A[j] >= B[j] for some j"
    ifail == 8 && return "maxiters (MAXPTS) < MINPTS"
    ifail == 9 && return "both ABSTOL < 0 and RELTOL < 0"
    ifail == 11 && return "SINGUL < 1 or SINGUL > NDIM"
    ifail == 12 && return "LOGF is not 0 or 1"
    ifail == 13 && return "maxiters (MAXPTS) too small for the first rule applications"
    ifail == 14 && return "ALPHA estimation failed or ALPHA <= -SINGUL (integral may not exist)"
    ifail == 15 && return "EMAX < 1"
    ifail == 17 && return "WRKSUB too small for the requested budget"
    return "unknown IFAIL code $ifail"
end

# ============================================================
# _check_params — port of DECHHR
#
# Returns (num, maxsub, keyf, wtleng, ifail).
# IFAIL codes:
#   2=bad KEY, 3=bad NDIM, 4=KEY=1 but NDIM≠2, 5=KEY=2 but NDIM≠3,
#   6=NUMFUN<1, 7=bad bounds, 8=MAXPTS<MINPTS, 9=both tols<0,
#   11=bad SINGUL, 12=bad LOGF, 13=MAXPTS too small, 15=EMAX<1, 17=WRKSUB too small
# ============================================================
function _check_params(
        maxdim::Int, ndim::Int, numfun::Int,
        a::AbstractVector, b::AbstractVector,
        singul::Int, logf::Int,
        minpts::Int, maxpts::Int,
        epsabs::Real, epsrel::Real,
        key::Int, emax::Int, wrksub::Int
    )
    if key < 0 || key > 4
        return 0, 0, 0, 0, 2
    end
    if ndim < 2 || ndim > maxdim
        return 0, 0, 0, 0, 3
    end
    if key == 1 && ndim != 2
        return 0, 0, 0, 0, 4
    end
    if key == 2 && ndim != 3
        return 0, 0, 0, 0, 5
    end

    keyf = if key == 0
        ndim == 2 ? 1 : ndim == 3 ? 2 : 3
    else
        key
    end

    num, wtleng = if keyf == 1
        65, 14
    elseif keyf == 2
        127, 13
    elseif keyf == 3
        pts = 1 + 2 * ndim + 6 * ndim^2 + 4 * ndim * (ndim - 1) * (ndim - 2) ÷ 3 + (1 << ndim)
        wt = ndim == 2 ? 8 : 9
        pts, wt
    else   # keyf == 4
        1 + 2 * ndim * (ndim + 2) + (1 << ndim), 6
    end

    if numfun < 1
        return num, 0, keyf, wtleng, 6
    end
    if any(a[j] >= b[j] for j in 1:ndim)
        return num, 0, keyf, wtleng, 7
    end
    if maxpts < minpts
        return num, 0, keyf, wtleng, 8
    end
    if epsabs < 0.0 && epsrel < 0.0
        return num, 0, keyf, wtleng, 9
    end
    if singul < 1 || singul > ndim
        return num, 0, keyf, wtleng, 11
    end
    if logf < 0 || logf > 1
        return num, 0, keyf, wtleng, 12
    end
    if maxpts < (singul + 2) * num
        return num, 0, keyf, wtleng, 13
    end
    if emax < 1
        return num, 0, keyf, wtleng, 15
    end

    # Compute MAXSUB (worst-case estimate)
    csing = (maxpts - num) ÷ ((singul + 1) * num)
    c1 = (maxpts - num * (1 + (singul + 1) * csing)) ÷ (2 * num)
    c1 = c1 == 0 ? -1 : c1
    maxsub = 2 + csing * singul + c1

    if maxsub > wrksub
        return num, maxsub, keyf, wtleng, 17
    end

    return num, maxsub, keyf, wtleng, 0
end

# ============================================================
# _heap_pop! — port of DETRHR DVFLAG=1
#
# Removes the top element (list[1]) from a max-heap of size `sbrgns`.
# `greate` is indexed by Julia subregion indices stored in `list`.
# Returns the new heap size (sbrgns - 1).
# ============================================================
function _heap_pop!(greate::AbstractVector, list::Vector{Int}, sbrgns::Int)::Int
    @inbounds begin
        great = greate[list[sbrgns]]   # value of last element (will replace top)
        new_sbrgns = sbrgns - 1
        subrgn = 1
        while true
            subtmp = 2 * subrgn
            subtmp > new_sbrgns && break
            # Choose the larger child
            if subtmp < new_sbrgns && greate[list[subtmp]] < greate[list[subtmp + 1]]
                subtmp += 1
            end
            # Sift down if necessary
            if great < greate[list[subtmp]]
                list[subrgn] = list[subtmp]
                subrgn = subtmp
            else
                break
            end
        end
        new_sbrgns > 0 && (list[subrgn] = list[new_sbrgns + 1])
    end
    return new_sbrgns
end

# ============================================================
# _heap_push! — port of DETRHR DVFLAG=2
#
# Inserts subregion `new_idx` (1-based Julia index) at heap position `new_pos`
# (= new heap size after insertion) and sifts up.
# ============================================================
function _heap_push!(
        greate::AbstractVector, list::Vector{Int},
        new_pos::Int, new_idx::Int
    )
    @inbounds begin
        great = greate[new_idx]
        subrgn = new_pos
        while true
            subtmp = subrgn ÷ 2
            subtmp < 1 && break
            if great > greate[list[subtmp]]
                list[subrgn] = list[subtmp]
                subrgn = subtmp
            else
                break
            end
        end
        list[subrgn] = new_idx
    end
    return nothing
end

# ============================================================
# _adaptive_integrate! — port of DEADHR
#
# Performs adaptive integration with singular-vertex extrapolation.
#
# Arguments:
#   ndim, numfun — problem dimensions
#   a, b         — bounds (length ndim)
#   maxsub       — maximum subregion count
#   funsub       — integrand: funsub(x, funvls) writes into funvls
#   singul       — singularity dimension
#   alpha        — singularity exponent (must be > -singul)
#   logf         — 1 if log term present, 0 otherwise
#   epsabs/rel   — tolerances
#   keyf         — resolved rule key (1..4)
#   num          — evaluations per subregion
#   wtleng       — number of rule generators
#   emax         — max extrapolation steps
#   minpts       — min evaluations before checking tolerance
#   maxpts       — max total evaluations
#
# Returns: (result, abserr, neval, ifail)
# ============================================================
function _adaptive_integrate!(
        ndim::Int, numfun::Int,
        a::AbstractVector, b::AbstractVector,
        maxsub::Int, funsub, ::Type{TV},
        singul::Int, alpha::Float64, logf::Int,
        epsabs::Float64, epsrel::Float64,
        keyf::Int, num::Int, wtleng::Int,
        emax::Int, minpts::Int, maxpts::Int
    ) where {TV}

    # ------------------------------------------------------------------
    # Initialise integration rule (DEINHR)
    # Rule coefficients are pure geometry: always Float64.
    # ------------------------------------------------------------------
    g = zeros(Float64, ndim, wtleng)
    w = zeros(Float64, 5, wtleng)
    errcof = zeros(Float64, 6)
    rulpts = zeros(Float64, wtleng)
    scales = zeros(Float64, 3, wtleng)
    norms = zeros(Float64, 3, wtleng)
    _init_rule!(ndim, keyf, wtleng, w, g, errcof, rulpts, scales, norms)

    # ------------------------------------------------------------------
    # Subregion data (Fortran index 0 → Julia 1, k → k+1)
    # Size maxsub+1 to hold singular (1) + maxsub non-singular regions.
    # Arrays holding integrand values are typed on TV; geometry stays Float64.
    # ------------------------------------------------------------------
    sz = maxsub + 1
    values = zeros(TV, numfun, sz)   # VALUES(NUMFUN, 0:MAXSUB)
    errors = zeros(TV, numfun, sz)   # ERRORS(NUMFUN, 0:MAXSUB)
    centrs = zeros(Float64, ndim, sz)   # CENTRS(NDIM,   0:MAXSUB)
    hwidts = zeros(Float64, ndim, sz)   # HWIDTS(NDIM,   0:MAXSUB)
    greate = zeros(TV, sz)           # GREATE(0:MAXSUB) — max error per region
    dir = zeros(Float64, sz)           # DIR(0:MAXSUB)    — direction index

    # Extrapolation data — value arrays on TV; coefficient arrays stay Float64
    q = zeros(TV, numfun, sz)        # Q(NUMFUN, 0:MAXSUB)
    u = zeros(TV, numfun, maxsub)    # U(NUMFUN, MAXSUB)
    e_mat = zeros(TV, numfun, maxsub)    # E(NUMFUN, MAXSUB)
    t_tab = zeros(TV, numfun, emax + 1)  # T(NUMFUN, 0:EMAX)
    ne = zeros(Float64, emax)              # NE(EMAX)
    beta = zeros(Float64, emax + 1, emax + 1) # BETA(0:EMAX, 0:EMAX)
    qold = zeros(TV, numfun)
    qnew = zeros(TV, numfun)
    unew = zeros(TV, numfun)
    uerr_s = zeros(TV, numfun)   # local UERR scratch (DEADHR's UERR)
    exterr = zeros(TV, numfun)

    # Heap and pointer arrays
    list = zeros(Int, maxsub)   # heap (1-based subregion indices)
    upoint = zeros(Int, sz)       # maps subregion → U-sequence index

    # Work arrays — x_work stays Float64 (integration points are coordinates)
    diff = zeros(TV, ndim)   # fourth-differences (may carry dual parts)
    order = collect(1:ndim)
    x_work = zeros(Float64, ndim)
    null_work = zeros(TV, numfun, 8)   # funsub output buffer
    fs_g = zeros(Float64, ndim)   # scratch for the fully-symmetric sum (generators)
    fs_x = zeros(Float64, ndim)   # scratch for the fully-symmetric sum (eval point)

    result = zeros(TV, numfun)
    abserr = zeros(TV, numfun)

    # ------------------------------------------------------------------
    # Initialise: singular region at Julia index 1 (Fortran index 0)
    # ------------------------------------------------------------------
    for j in 1:ndim
        centrs[j, 1] = (a[j] + b[j]) / 2.0
        hwidts[j, 1] = abs(b[j] - a[j]) / 2.0
        order[j] = j
    end
    dir[1] = Float64(-singul)   # DIR(0) = -SINGUL

    # Evaluate rule over the initial whole region
    let dref = Ref(dir[1]), gref = Ref{TV}(zero(TV))
        _eval_rule!(
            ndim, view(centrs, :, 1), view(hwidts, :, 1),
            wtleng, g, w, errcof, numfun, funsub, scales, norms,
            x_work, null_work, fs_g, fs_x,
            view(values, :, 1), view(errors, :, 1), dref, gref, diff, order
        )
        dir[1] = dref[]
        greate[1] = gref[]
    end

    for j in 1:numfun
        result[j] = values[j, 1]
        abserr[j] = errors[j, 1]
        q[j, 1] = values[j, 1]   # Q(J,0) = VALUES(J,0)
    end
    neval = num
    first_call = Ref(true)
    numu = 0
    sbrgns = 0

    # LIST(1) = 0 in Fortran → list[1] = 1 (singular region) as initial sentinel
    list[1] = 1

    # ------------------------------------------------------------------
    # Main adaptive loop (Fortran label 50)
    # ------------------------------------------------------------------
    ifail = 1   # default: budget exceeded

    while true
        # Top of heap (1-based Julia index)
        top = list[1]

        # Choose which region to subdivide
        # Singular region (index 1) if its error > heap top
        vacant = greate[top] < greate[1] ? 1 : top

        direct_i = round(Int, dir[vacant])     # integer direction value
        ndiv = max(2, -direct_i + 1)       # SINGUL+1 or 2

        # Budget check
        if neval + ndiv * num > maxpts
            ifail = 1
            break
        end

        # Next free Julia slot: Fortran POINTR=SBRGNS+1 → Julia sbrgns+2
        pointr = sbrgns + 2

        # Remove vacant from heap (only for non-singular regions)
        if vacant > 1
            sbrgns = _heap_pop!(greate, list, sbrgns)
        end

        # ------------------------------------------------------------------
        # Subdivide (DESBHR)
        # ------------------------------------------------------------------
        if direct_i > 0
            # Non-singular bisection in direction `direct_i`
            for j in 1:ndim
                hwidts[j, pointr] = hwidts[j, vacant]
                centrs[j, pointr] = centrs[j, vacant]
            end
            hwidts[direct_i, vacant] /= 2.0
            hwidts[direct_i, pointr] = hwidts[direct_i, vacant]
            centrs[direct_i, pointr] = centrs[direct_i, vacant] + hwidts[direct_i, pointr]
            centrs[direct_i, vacant] -= hwidts[direct_i, vacant]
            dir[pointr] = dir[vacant]   # DIR(POINTR) = DIRECT
        else
            # Singular subdivision: cut first |direct_i| directions
            neg_d = -direct_i
            for i in 1:neg_d
                nxt = order[i]
                hwidts[nxt, vacant] /= 2.0
                for j in 1:ndim
                    hwidts[j, pointr + i - 1] = hwidts[j, vacant]
                    centrs[j, pointr + i - 1] = centrs[j, vacant]
                end
                dir[pointr + i - 1] = Float64(nxt)
                centrs[nxt, pointr + i - 1] = centrs[nxt, vacant] + hwidts[nxt, vacant]
                centrs[nxt, vacant] -= hwidts[nxt, vacant]
            end
        end
        dir[vacant] = Float64(direct_i)   # DIR(VACANT) = DIRECT (always)

        # ------------------------------------------------------------------
        # Determine UPDATE index
        # ------------------------------------------------------------------
        if vacant == 1   # singular region → new U-term
            numu += 1
            update = numu
        else
            update = upoint[vacant]
        end

        # ------------------------------------------------------------------
        # Evaluate rule on VACANT position
        # ------------------------------------------------------------------
        if vacant == 1
            fill!(uerr_s, zero(TV))
            fill!(unew, zero(TV))
        else
            for j in 1:numfun
                uerr_s[j] = -errors[j, vacant]
                unew[j] = -values[j, vacant]
            end
        end

        let dref = Ref(dir[vacant]), gref = Ref{TV}(zero(TV))
            _eval_rule!(
                ndim, view(centrs, :, vacant), view(hwidts, :, vacant),
                wtleng, g, w, errcof, numfun, funsub, scales, norms,
                x_work, null_work, fs_g, fs_x,
                view(values, :, vacant), view(errors, :, vacant),
                dref, gref, diff, order
            )
            dir[vacant] = dref[]
            greate[vacant] = gref[]
        end

        if vacant == 1
            # Q(J, NUMU) = VALUES(J, 0)  →  q[j, numu+1] = values[j, 1]
            for j in 1:numfun
                q[j, numu + 1] = values[j, 1]
            end
        else
            upoint[vacant] = update
            for j in 1:numfun
                uerr_s[j] += errors[j, vacant]
                unew[j] += values[j, vacant]
            end
        end

        # ------------------------------------------------------------------
        # Evaluate rule on the remaining new regions (POINTR..POINTR+NDIV-2)
        # ------------------------------------------------------------------
        for i in 2:ndiv
            idx2 = pointr + i - 2
            let dref = Ref(dir[idx2]), gref = Ref{TV}(zero(TV))
                _eval_rule!(
                    ndim, view(centrs, :, idx2), view(hwidts, :, idx2),
                    wtleng, g, w, errcof, numfun, funsub, scales, norms,
                    x_work, null_work, fs_g, fs_x,
                    view(values, :, idx2), view(errors, :, idx2),
                    dref, gref, diff, order
                )
                dir[idx2] = dref[]
                greate[idx2] = gref[]
            end
            upoint[idx2] = update
            for j in 1:numfun
                uerr_s[j] += errors[j, idx2]
                unew[j] += values[j, idx2]
            end
        end

        # ------------------------------------------------------------------
        # Accumulate U and E series; set QNEW / QOLD for extrapolation
        # ------------------------------------------------------------------
        for j in 1:numfun
            qnew[j] = q[j, numu + 1]     # Q(J, NUMU)
            qold[j] = q[j, numu]        # Q(J, NUMU-1)
            u[j, update] += unew[j]
            e_mat[j, update] += uerr_s[j]
        end

        # ------------------------------------------------------------------
        # Extrapolation (DEXTHR)
        # UPDATE=0 signals a new step; else it's a correction of term `update`
        # ------------------------------------------------------------------
        ext_update = (vacant == 1) ? 0 : update
        _extrapolate!(
            numfun, alpha, logf, singul, 1, numu, t_tab, ext_update,
            unew, qnew, qold, first_call, emax, e_mat,
            result, abserr, exterr, ne, beta
        )
        first_call[] = false

        neval += ndiv * num

        # Update singular region's error estimate to the extrapolation error
        greate[1] = zero(TV)
        for j in 1:numfun
            greate[1] = max(greate[1], exterr[j])
            errors[j, 1] = exterr[j]
        end

        # ------------------------------------------------------------------
        # Insert new/updated regions into heap
        # ------------------------------------------------------------------
        if vacant > 1
            sbrgns += 1
            _heap_push!(greate, list, sbrgns, vacant)
        end
        for i in 2:ndiv
            idx2 = pointr + i - 2
            sbrgns += 1
            _heap_push!(greate, list, sbrgns, idx2)
        end

        # ------------------------------------------------------------------
        # Termination check
        # ------------------------------------------------------------------
        if neval < minpts
            continue   # not enough evaluations yet
        end

        converged = true
        for j in 1:numfun
            if abserr[j] > epsrel * abs(result[j]) && abserr[j] > epsabs
                converged = false
                break
            end
        end
        if converged
            ifail = 0
            break
        end
    end   # while true

    return result, abserr, neval, ifail
end
