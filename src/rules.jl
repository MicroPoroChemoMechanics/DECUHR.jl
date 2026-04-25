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
# rules.jl — Integration rules
# Port of: D132RE, D113RE, D09HRE, D07HRE, DEINHR, DEFSHR, DERLHR
#
# References:
#   Espelid & Genz, Numerical Algorithms 8 (1994) 201-220
#   Genz & Malik, SIAM J Numer. Anal. 20 (1983) 580-588

# ============================================================
# D132RE — 2D degree-13 rule (65 points, WTLENG = 14)
# ============================================================
function _init_rule_d13_2d!(w, g, errcof, rulpts)
    dim2g = [
        0.2517129343453109, 0.7013933644534266, 0.9590960631619962,
        0.9956010478552127, 0.5, 0.1594544658297559,
        0.3808991135940188, 0.6582769255267192, 0.8761473165029315,
        0.998243184053198, 0.9790222658168462, 0.6492284325645389,
        0.8727421201131239, 0.3582614645881228, 0.5666666666666666,
        0.2077777777777778,
    ]

    # dim2w[i,j]: generator i (1..14), rule j (1=basic, 2..5=null rules)
    dim2w = [
        3.37969236013446e-2   3.213775489050763e-1   3.372900883288987e-1  -8.264123822525677e-1   6.539094339575232e-1;
        9.508589607597761e-2  -1.767341636743844e-1  -1.644903060344491e-1   3.06583861409436e-1  -2.041614154424632e-1;
        1.176006468056962e-1   7.347600537466072e-2   7.707849911634622e-2   2.389292538329435e-3  -1.74698151579499e-1;
        2.65777458632695e-2  -3.638022004364754e-2  -3.80447835850631e-2  -1.343024157997222e-1   3.937939671417803e-2;
        1.70144177020064e-2   2.125297922098712e-2   2.223559940380806e-2   8.8333668405339e-2   6.974520545933992e-3;
        0.0                     1.460984204026913e-1   1.480693879765931e-1   0.0                     0.0;
        1.62659309863741e-2   1.747613286152099e-2   4.467143702185814e-6   9.786283074168292e-4   6.667702171778258e-3;
        1.344892658526199e-1   1.444954045641582e-1   1.50894476707413e-1  -1.319227889147519e-1   5.512960621544304e-2;
        1.328032165460149e-1   1.307687976001325e-4   3.647200107516215e-5   7.99001220015063e-3   5.443846381278607e-2;
        5.63747476999187e-2   5.380992313941161e-4   5.77719899901388e-4   3.391747079760626e-3   2.310903863953934e-2;
        3.9082790813105e-3   1.042259576889814e-4   1.041757313688177e-4   2.294915718283264e-3   1.506937747477189e-2;
        3.01279877743215e-2  -1.401152865045733e-3  -1.452822267047819e-3  -1.358584986119197e-2  -6.05702164890189e-2;
        1.030873234689166e-1   8.041788181514763e-3   8.338339968783705e-3   4.025866859057809e-2   4.225737654686337e-2;
        6.25e-2  -1.420416552759383e-1  -1.47279632923196e-1   3.760268580063992e-3   2.561989142123099e-2
    ]

    # Fortran: W(J,I) = DIM2W(I,J) → w[j,i] = dim2w[i,j]
    for i in 1:14, j in 1:5
        w[j, i] = dim2w[i, j]
    end

    fill!(g, 0.0)
    g[1, 2] = dim2g[1];  g[1, 3] = dim2g[2];  g[1, 4] = dim2g[3];  g[1, 5] = dim2g[4]
    g[1, 6] = dim2g[5]
    g[1, 7] = dim2g[6];  g[2, 7] = g[1, 7]
    g[1, 8] = dim2g[7];  g[2, 8] = g[1, 8]
    g[1, 9] = dim2g[8];  g[2, 9] = g[1, 9]
    g[1, 10] = dim2g[9];  g[2, 10] = g[1, 10]
    g[1, 11] = dim2g[10]; g[2, 11] = g[1, 11]
    g[1, 12] = dim2g[11]; g[2, 12] = dim2g[12]
    g[1, 13] = dim2g[13]; g[2, 13] = dim2g[14]
    g[1, 14] = dim2g[15]; g[2, 14] = dim2g[16]

    rulpts[1] = 1.0
    for i in 2:11
        rulpts[i] = 4.0
    end
    rulpts[12] = 8.0; rulpts[13] = 8.0; rulpts[14] = 8.0

    errcof[1] = 10.0; errcof[2] = 10.0; errcof[3] = 1.0
    errcof[4] = 5.0; errcof[5] = 0.5; errcof[6] = 0.25
    return nothing
end

# ============================================================
# D113RE — 3D degree-11 rule (127 points, WTLENG = 13)
# ============================================================
function _init_rule_d11_3d!(w, g, errcof, rulpts)
    dim3g = [
        0.19, 0.5, 0.75,
        0.8, 0.9949999999999999, 0.99873449983514,
        0.7793703685672423, 0.9999698993088767, 0.7902637224771788,
        0.4403396687650737, 0.4378478459006862, 0.9549373822794593,
        0.9661093133630748, 0.4577105877763134,
    ]

    dim3w = [
        7.923078151105734e-3   1.715006248224684e+0   1.936014978949526e+0   5.17082819560576e-1   2.05440450381852e+0;
        6.79717739278808e-2  -3.755893815889209e-1  -3.673449403754268e-1   1.445269144914044e-2   1.37775998849012e-2;
        1.086986538805825e-3   1.488632145140549e-1   2.929778657898176e-2  -3.601489663995932e-1  -5.76806291790441e-1;
        1.838633662212829e-1  -2.497046640620823e-1  -1.151883520260315e-1   3.628307003418485e-1   3.726835047700328e-2;
        3.362119777829031e-2   1.792501419135204e-1   5.086658220872218e-2   7.148802650872729e-3   6.814878939777219e-3;
        1.013751123334062e-2   3.44612675897389e-3   4.453911087786469e-2  -9.222852896022966e-2   5.723169733851849e-2;
        1.687648683985235e-3  -5.140483185555825e-3  -2.2878282571259e-2   1.719339732471725e-2  -4.493018743811285e-2;
        1.346468564512807e-1   6.536017839876425e-3   2.908926216345833e-2  -1.02141653746035e-1   2.729236573866348e-2;
        1.750145884600386e-3  -6.5134549392297e-4  -2.898884350669207e-3  -7.504397861080493e-3   3.54747395055699e-4;
        7.752336383837454e-2  -6.304672433547204e-3  -2.805963413307495e-2   1.648362537726711e-2   1.571366799739551e-2;
        2.461864902770251e-1   1.266959399788263e-2   5.638741361145884e-2   5.234610158469334e-2   4.990099219278567e-2;
        6.797944868483039e-2  -5.454241018647931e-3  -2.427469611942451e-2   1.445432331613066e-2   1.37791555266677e-2;
        1.419962823300713e-2   4.826995274768427e-3   2.148307034182882e-2   3.019236275367777e-3   2.878206423099872e-3
    ]

    for i in 1:13, j in 1:5
        w[j, i] = dim3w[i, j]
    end

    fill!(g, 0.0)
    g[1, 2] = dim3g[1]
    g[1, 3] = dim3g[2]
    g[1, 4] = dim3g[3]
    g[1, 5] = dim3g[4]
    g[1, 6] = dim3g[5]
    g[1, 7] = dim3g[6];  g[2, 7] = g[1, 7]
    g[1, 8] = dim3g[7];  g[2, 8] = g[1, 8]
    g[1, 9] = dim3g[8];  g[2, 9] = g[1, 9];  g[3, 9] = g[1, 9]
    g[1, 10] = dim3g[9];  g[2, 10] = g[1, 10]; g[3, 10] = g[1, 10]
    g[1, 11] = dim3g[10]; g[2, 11] = g[1, 11]; g[3, 11] = g[1, 11]
    g[1, 12] = dim3g[12]; g[2, 12] = dim3g[11]; g[3, 12] = g[2, 12]
    g[1, 13] = dim3g[13]; g[2, 13] = g[1, 13]; g[3, 13] = dim3g[14]

    rulpts[1] = 1.0; rulpts[2] = 6.0; rulpts[3] = 6.0; rulpts[4] = 6.0
    rulpts[5] = 6.0; rulpts[6] = 6.0; rulpts[7] = 12.0; rulpts[8] = 12.0
    rulpts[9] = 8.0; rulpts[10] = 8.0; rulpts[11] = 8.0; rulpts[12] = 24.0; rulpts[13] = 24.0

    errcof[1] = 4.0; errcof[2] = 4.0; errcof[3] = 0.5
    errcof[4] = 3.0; errcof[5] = 0.5; errcof[6] = 0.25
    return nothing
end

# ============================================================
# D09HRE — degree-9 rule, NDIM ≥ 2  (WTLENG = 9, or 8 for NDIM=2)
# ============================================================
function _init_rule_d9_nd!(ndim, wtleng, w, g, errcof, rulpts)
    twondm = Float64(1 << ndim)

    fill!(g, 0.0); fill!(w, 0.0)
    for j in 1:wtleng
        rulpts[j] = Float64(2 * ndim)
    end
    rulpts[wtleng] = twondm
    if ndim > 2
        rulpts[8] = Float64(4 * ndim * (ndim - 1) * (ndim - 2) ÷ 3)
    end
    rulpts[7] = Float64(4 * ndim * (ndim - 1))
    rulpts[6] = Float64(2 * ndim * (ndim - 1))
    rulpts[1] = 1.0

    # Squared generator parameters
    lam0 = 0.4707
    lam1 = 4.0 / (15.0 - 5.0 / lam0)
    ratio = (1.0 - lam1 / lam0) / 27.0
    lam2 = (5.0 - 7.0 * lam1 - 35.0 * ratio) / (7.0 - 35.0 * lam1 / 3.0 - 35.0 * ratio / lam0)
    ratio = ratio * (1.0 - lam2 / lam0) / 3.0
    lam3 = (7.0 - 9.0 * (lam2 + lam1) + 63.0 * lam2 * lam1 / 5.0 - 63.0 * ratio) /
        (9.0 - 63.0 * (lam2 + lam1) / 5.0 + 21.0 * lam2 * lam1 - 63.0 * ratio / lam0)
    lamp = 0.0625

    w8(r) = ndim > 2 ? r : 0.0   # helper: w[r,8] only exists for ndim > 2

    # Basic rule weights (row 1)
    w[1, wtleng] = 1.0 / (3.0 * lam0)^4 / twondm
    if ndim > 2
        w[1, 8] = (1.0 - 1.0 / (3.0 * lam0)) / (6.0 * lam1)^3
    end
    w[1, 7] = (1.0 - 7.0 * (lam0 + lam1) / 5.0 + 7.0 * lam0 * lam1 / 3.0) /
        (84.0 * lam1 * lam2 * (lam2 - lam0) * (lam2 - lam1))
    w[1, 6] = (1.0 - 7.0 * (lam0 + lam2) / 5.0 + 7.0 * lam0 * lam2 / 3.0) /
        (84.0 * lam1^2 * (lam1 - lam0) * (lam1 - lam2)) -
        w[1, 7] * lam2 / lam1 - 2.0 * (ndim - 2) * w8(w[1, 8])
    w[1, 4] = (
        1.0 - 9.0 * (
            (lam0 + lam1 + lam2) / 7.0
                - (lam0 * lam1 + lam0 * lam2 + lam1 * lam2) / 5.0
        )
            - 3.0 * lam0 * lam1 * lam2
    ) /
        (18.0 * lam3 * (lam3 - lam0) * (lam3 - lam1) * (lam3 - lam2))
    w[1, 3] = (
        1.0 - 9.0 * (
            (lam0 + lam1 + lam3) / 7.0
                - (lam0 * lam1 + lam0 * lam3 + lam1 * lam3) / 5.0
        )
            - 3.0 * lam0 * lam1 * lam3
    ) /
        (18.0 * lam2 * (lam2 - lam0) * (lam2 - lam1) * (lam2 - lam3)) - 2.0 * (ndim - 1) * w[1, 7]
    w[1, 2] = (
        1.0 - 9.0 * (
            (lam0 + lam2 + lam3) / 7.0
                - (lam0 * lam2 + lam0 * lam3 + lam2 * lam3) / 5.0
        )
            - 3.0 * lam0 * lam2 * lam3
    ) /
        (18.0 * lam1 * (lam1 - lam0) * (lam1 - lam2) * (lam1 - lam3)) -
        2.0 * (ndim - 1) * (w[1, 7] + w[1, 6] + (ndim - 2) * w8(w[1, 8]))

    # Null rule 1 (row 2)
    w[2, wtleng] = 1.0 / (108.0 * lam0^4) / twondm
    if ndim > 2
        w[2, 8] = (1.0 - 27.0 * twondm * w[2, wtleng] * lam0^3) / (6.0 * lam1)^3
    end
    w[2, 7] = (1.0 - 5.0 * lam1 / 3.0 - 15.0 * twondm * w[2, wtleng] * lam0^2 * (lam0 - lam1)) /
        (60.0 * lam1 * lam2 * (lam2 - lam1))
    w[2, 6] = (1.0 - 9.0 * (8.0 * lam1 * lam2 * w[2, 7] + twondm * w[2, wtleng] * lam0^2)) /
        (36.0 * lam1^2) - 2.0 * (ndim - 2) * w8(w[2, 8])
    w[2, 4] = (
        1.0 - 7.0 * (
            (lam1 + lam2) / 5.0 - lam1 * lam2 / 3.0
                + twondm * w[2, wtleng] * lam0 * (lam0 - lam1) * (lam0 - lam2)
        )
    ) /
        (14.0 * lam3 * (lam3 - lam1) * (lam3 - lam2))
    w[2, 3] = (
        1.0 - 7.0 * (
            (lam1 + lam3) / 5.0 - lam1 * lam3 / 3.0
                + twondm * w[2, wtleng] * lam0 * (lam0 - lam1) * (lam0 - lam3)
        )
    ) /
        (14.0 * lam2 * (lam2 - lam1) * (lam2 - lam3)) - 2.0 * (ndim - 1) * w[2, 7]
    w[2, 2] = (
        1.0 - 7.0 * (
            (lam2 + lam3) / 5.0 - lam2 * lam3 / 3.0
                + twondm * w[2, wtleng] * lam0 * (lam0 - lam2) * (lam0 - lam3)
        )
    ) /
        (14.0 * lam1 * (lam1 - lam2) * (lam1 - lam3)) -
        2.0 * (ndim - 1) * (w[2, 7] + w[2, 6] + (ndim - 2) * w8(w[2, 8]))

    # Null rule 2 (row 3)
    w[3, wtleng] = 5.0 / (324.0 * lam0^4) / twondm
    if ndim > 2
        w[3, 8] = (1.0 - 27.0 * twondm * w[3, wtleng] * lam0^3) / (6.0 * lam1)^3
    end
    w[3, 7] = (1.0 - 5.0 * lam1 / 3.0 - 15.0 * twondm * w[3, wtleng] * lam0^2 * (lam0 - lam1)) /
        (60.0 * lam1 * lam2 * (lam2 - lam1))
    w[3, 6] = (1.0 - 9.0 * (8.0 * lam1 * lam2 * w[3, 7] + twondm * w[3, wtleng] * lam0^2)) /
        (36.0 * lam1^2) - 2.0 * (ndim - 2) * w8(w[3, 8])
    w[3, 5] = (
        1.0 - 7.0 * (
            (lam1 + lam2) / 5.0 - lam1 * lam2 / 3.0
                + twondm * w[3, wtleng] * lam0 * (lam0 - lam1) * (lam0 - lam2)
        )
    ) /
        (14.0 * lamp * (lamp - lam1) * (lamp - lam2))
    w[3, 3] = (
        1.0 - 7.0 * (
            (lam1 + lamp) / 5.0 - lam1 * lamp / 3.0
                + twondm * w[3, wtleng] * lam0 * (lam0 - lam1) * (lam0 - lamp)
        )
    ) /
        (14.0 * lam2 * (lam2 - lam1) * (lam2 - lamp)) - 2.0 * (ndim - 1) * w[3, 7]
    w[3, 2] = (
        1.0 - 7.0 * (
            (lam2 + lamp) / 5.0 - lam2 * lamp / 3.0
                + twondm * w[3, wtleng] * lam0 * (lam0 - lam2) * (lam0 - lamp)
        )
    ) /
        (14.0 * lam1 * (lam1 - lam2) * (lam1 - lamp)) -
        2.0 * (ndim - 1) * (w[3, 7] + w[3, 6] + (ndim - 2) * w8(w[3, 8]))

    # Null rule 3 (row 4)
    w[4, wtleng] = 2.0 / (81.0 * lam0^4) / twondm
    if ndim > 2
        w[4, 8] = (2.0 - 27.0 * twondm * w[4, wtleng] * lam0^3) / (6.0 * lam1)^3
    end
    w[4, 7] = (2.0 - 15.0 * lam1 / 9.0 - 15.0 * twondm * w[4, wtleng] * lam0 * (lam0 - lam1)) /
        (60.0 * lam1 * lam2 * (lam2 - lam1))
    w[4, 6] = (1.0 - 9.0 * (8.0 * lam1 * lam2 * w[4, 7] + twondm * w[4, wtleng] * lam0^2)) /
        (36.0 * lam1^2) - 2.0 * (ndim - 2) * w8(w[4, 8])
    w[4, 4] = (
        2.0 - 7.0 * (
            (lam1 + lam2) / 5.0 - lam1 * lam2 / 3.0
                + twondm * w[4, wtleng] * lam0 * (lam0 - lam1) * (lam0 - lam2)
        )
    ) /
        (14.0 * lam3 * (lam3 - lam1) * (lam3 - lam2))
    w[4, 3] = (
        2.0 - 7.0 * (
            (lam1 + lam3) / 5.0 - lam1 * lam3 / 3.0
                + twondm * w[4, wtleng] * lam0 * (lam0 - lam1) * (lam0 - lam3)
        )
    ) /
        (14.0 * lam2 * (lam2 - lam1) * (lam2 - lam3)) - 2.0 * (ndim - 1) * w[4, 7]
    w[4, 2] = (
        2.0 - 7.0 * (
            (lam2 + lam3) / 5.0 - lam2 * lam3 / 3.0
                + twondm * w[4, wtleng] * lam0 * (lam0 - lam2) * (lam0 - lam3)
        )
    ) /
        (14.0 * lam1 * (lam1 - lam2) * (lam1 - lam3)) -
        2.0 * (ndim - 1) * (w[4, 7] + w[4, 6] + (ndim - 2) * w8(w[4, 8]))

    # Null rule 4 (row 5): only non-zero at position 2
    w[5, 2] = 1.0 / (6.0 * lam1)

    # Generators (sqrt to recover actual abscissas from squared values)
    lam0 = sqrt(lam0); lam1 = sqrt(lam1); lam2 = sqrt(lam2)
    lam3 = sqrt(lam3); lamp = sqrt(lamp)

    for i in 1:ndim
        g[i, wtleng] = lam0
    end
    if ndim > 2
        g[1, 8] = lam1; g[2, 8] = lam1; g[3, 8] = lam1
    end
    g[1, 7] = lam1; g[2, 7] = lam2
    g[1, 6] = lam1; g[2, 6] = lam1
    g[1, 5] = lamp
    g[1, 4] = lam3
    g[1, 3] = lam2
    g[1, 2] = lam1

    # Finalise: null rule weights = difference from basic; scale basic rule
    w[1, 1] = twondm
    for jj in 2:5
        for i in 2:wtleng
            w[jj, i] -= w[1, i]
            w[jj, 1] -= rulpts[i] * w[jj, i]
        end
    end
    for i in 2:wtleng
        w[1, i] *= twondm
        w[1, 1] -= rulpts[i] * w[1, i]
    end

    errcof[1] = 5.0; errcof[2] = 5.0; errcof[3] = 1.0
    errcof[4] = 5.0; errcof[5] = 0.5; errcof[6] = 0.25
    return nothing
end

# ============================================================
# D07HRE — degree-7 rule, NDIM ≥ 2  (WTLENG = 6)
# ============================================================
function _init_rule_d7_nd!(ndim, wtleng, w, g, errcof, rulpts)
    twondm = Float64(1 << ndim)

    fill!(g, 0.0); fill!(w, 0.0)
    for j in 1:wtleng
        rulpts[j] = Float64(2 * ndim)
    end
    rulpts[wtleng] = twondm
    rulpts[wtleng - 1] = Float64(2 * ndim * (ndim - 1))
    rulpts[1] = 1.0

    lam0 = 0.4707
    lamp = 0.5625
    lam1 = 4.0 / (15.0 - 5.0 / lam0)
    ratio = (1.0 - lam1 / lam0) / 27.0
    lam2 = (5.0 - 7.0 * lam1 - 35.0 * ratio) / (7.0 - 35.0 * lam1 / 3.0 - 35.0 * ratio / lam0)

    # Basic rule weights (row 1)
    w[1, 6] = 1.0 / (3.0 * lam0)^3 / twondm
    w[1, 5] = (1.0 - 5.0 * lam0 / 3.0) / (60.0 * (lam1 - lam0) * lam1^2)
    w[1, 3] = (1.0 - 5.0 * lam2 / 3.0 - 5.0 * twondm * w[1, 6] * lam0 * (lam0 - lam2)) /
        (10.0 * lam1 * (lam1 - lam2)) - 2.0 * (ndim - 1) * w[1, 5]
    w[1, 2] = (1.0 - 5.0 * lam1 / 3.0 - 5.0 * twondm * w[1, 6] * lam0 * (lam0 - lam1)) /
        (10.0 * lam2 * (lam2 - lam1))

    # Null rule 1 (row 2)
    w[2, 6] = 1.0 / (36.0 * lam0^3) / twondm
    w[2, 5] = (1.0 - 9.0 * twondm * w[2, 6] * lam0^2) / (36.0 * lam1^2)
    w[2, 3] = (1.0 - 5.0 * lam2 / 3.0 - 5.0 * twondm * w[2, 6] * lam0 * (lam0 - lam2)) /
        (10.0 * lam1 * (lam1 - lam2)) - 2.0 * (ndim - 1) * w[2, 5]
    w[2, 2] = (1.0 - 5.0 * lam1 / 3.0 - 5.0 * twondm * w[2, 6] * lam0 * (lam0 - lam1)) /
        (10.0 * lam2 * (lam2 - lam1))

    # Null rule 2 (row 3)
    w[3, 6] = 5.0 / (108.0 * lam0^3) / twondm
    w[3, 5] = (1.0 - 9.0 * twondm * w[3, 6] * lam0^2) / (36.0 * lam1^2)
    w[3, 3] = (1.0 - 5.0 * lamp / 3.0 - 5.0 * twondm * w[3, 6] * lam0 * (lam0 - lamp)) /
        (10.0 * lam1 * (lam1 - lamp)) - 2.0 * (ndim - 1) * w[3, 5]
    w[3, 4] = (1.0 - 5.0 * lam1 / 3.0 - 5.0 * twondm * w[3, 6] * lam0 * (lam0 - lam1)) /
        (10.0 * lamp * (lamp - lam1))

    # Null rule 3 (row 4)
    w[4, 6] = 1.0 / (54.0 * lam0^3) / twondm
    w[4, 5] = (1.0 - 18.0 * twondm * w[4, 6] * lam0^2) / (72.0 * lam1^2)
    w[4, 3] = (1.0 - 10.0 * lam2 / 3.0 - 10.0 * twondm * w[4, 6] * lam0 * (lam0 - lam2)) /
        (20.0 * lam1 * (lam1 - lam2)) - 2.0 * (ndim - 1) * w[4, 5]
    w[4, 2] = (1.0 - 10.0 * lam1 / 3.0 - 10.0 * twondm * w[4, 6] * lam0 * (lam0 - lam1)) /
        (20.0 * lam2 * (lam2 - lam1))

    lam0 = sqrt(lam0); lam1 = sqrt(lam1); lam2 = sqrt(lam2); lamp = sqrt(lamp)

    for i in 1:ndim
        g[i, wtleng] = lam0
    end
    g[1, wtleng - 1] = lam1; g[2, wtleng - 1] = lam1
    g[1, wtleng - 4] = lam2    # position 2 (wtleng=6)
    g[1, wtleng - 3] = lam1    # position 3
    g[1, wtleng - 2] = lamp    # position 4

    w[1, 1] = twondm
    for jj in 2:5
        for i in 2:wtleng
            w[jj, i] -= w[1, i]
            w[jj, 1] -= rulpts[i] * w[jj, i]
        end
    end
    for i in 2:wtleng
        w[1, i] *= twondm
        w[1, 1] -= rulpts[i] * w[1, i]
    end

    errcof[1] = 5.0; errcof[2] = 5.0; errcof[3] = 1.0
    errcof[4] = 5.0; errcof[5] = 0.5; errcof[6] = 0.25
    return nothing
end

# ============================================================
# DEINHR — select rule, compute SCALES and NORMS
# ============================================================
function _init_rule!(ndim, key, wtleng, w, g, errcof, rulpts, scales, norms)
    if key == 1
        _init_rule_d13_2d!(w, g, errcof, rulpts)
    elseif key == 2
        _init_rule_d11_3d!(w, g, errcof, rulpts)
    elseif key == 3
        _init_rule_d9_nd!(ndim, wtleng, w, g, errcof, rulpts)
    elseif key == 4
        _init_rule_d7_nd!(ndim, wtleng, w, g, errcof, rulpts)
    end

    we = zeros(wtleng)
    for k in 1:3
        for i in 1:wtleng
            scales[k, i] = abs(w[k + 1, i]) > 0.0 ? -w[k + 2, i] / w[k + 1, i] : 100.0
            for j in 1:wtleng
                we[j] = w[k + 2, j] + scales[k, i] * w[k + 1, j]
            end
            s = 0.0
            for j in 1:wtleng
                s += rulpts[j] * abs(we[j])
            end
            norms[k, i] = 2.0^ndim / s
        end
    end
    return nothing
end

# ============================================================
# DEFSHR — fully symmetric sum over all sign changes and permutations.
# Works on a local copy of g_col; the original is not modified.
# funsub(x, funvls) evaluates the integrand in-place.
# ============================================================
function _fully_symmetric_sum!(fulsms, center, hwidth, g_col, ndim, numfun, funsub, funvls)
    fill!(fulsms, 0.0)
    g = copy(g_col)     # local mutable copy
    x = zeros(Float64, ndim)

    while true
        # Set evaluation point from current permutation of g  (label 20)
        for i in 1:ndim
            x[i] = center[i] + g[i] * hwidth[i]
        end

        # Enumerate all sign changes for this permutation  (label 40 / DO 60)
        while true
            funsub(x, funvls)
            for j in 1:numfun
                fulsms[j] += funvls[j]
            end

            # Flip signs one coordinate at a time; on new negative → repeat evaluation
            jumped = false
            for i in 1:ndim
                g[i] = -g[i]
                x[i] = center[i] + g[i] * hwidth[i]
                if g[i] < 0.0
                    jumped = true
                    break
                end
            end
            jumped || break     # all signs restored → done with this permutation
        end

        # Find next distinct permutation in reverse-lexicographic order  (DO 80)
        found = false
        for i in 2:ndim
            if g[i - 1] > g[i]
                gi = g[i]
                ixchng = i - 1
                lxchng = 0
                for l in 1:((i - 1) ÷ 2)
                    gl = g[l]
                    g[l] = g[i - l]
                    g[i - l] = gl
                    if gl <= gi
                        ixchng -= 1
                    end
                    if g[l] > gi
                        lxchng = l
                    end
                end
                if g[ixchng] <= gi
                    ixchng = lxchng
                end
                g[i] = g[ixchng]
                g[ixchng] = gi
                found = true
                break
            end
        end
        found || break      # all permutations exhausted
    end
    return nothing
end

# ============================================================
# DERLHR — evaluate integration rule over one sub-region, estimate error.
#
# Conventions (matching Fortran):
#   direct_ref[] ≥ 0 → regular region, subdivision axis is written back
#   direct_ref[] < 0 → singular region, ORDER is sorted by DIFF magnitude
#   greate_ref[] ← max error over NUMFUN components
#   null_work     :: Matrix{Float64}(numfun, 8) — work storage
# ============================================================
function _eval_rule!(
        ndim, center, hwidth, wtleng, g, w, errcof,
        numfun, funsub, scales, norms, x, null_work,
        basval, rgnerr, direct_ref, greate_ref, diff, order
    )

    direct = direct_ref[]

    # Volume
    rgnvol = 1.0
    divaxn = 1
    for i in 1:ndim
        rgnvol *= hwidth[i]
        x[i] = center[i]
        if hwidth[i] > hwidth[divaxn]
            divaxn = i
        end
    end

    # f at center
    funsub(x, rgnerr)   # rgnerr used as temp here
    for j in 1:numfun
        basval[j] = w[1, 1] * rgnerr[j]
        for k in 1:4
            null_work[j, k] = w[k + 1, 1] * rgnerr[j]
        end
    end

    # Fourth differences along each axis (generators at positions 2 and 3)
    ratio = (g[1, 3] / g[1, 2])^2
    difmax = 0.0
    for i in 1:ndim
        x[i] = center[i] - hwidth[i] * g[1, 2]
        funsub(x, view(null_work, :, 5))
        x[i] = center[i] + hwidth[i] * g[1, 2]
        funsub(x, view(null_work, :, 6))
        x[i] = center[i] - hwidth[i] * g[1, 3]
        funsub(x, view(null_work, :, 7))
        x[i] = center[i] + hwidth[i] * g[1, 3]
        funsub(x, view(null_work, :, 8))
        x[i] = center[i]

        difsum = 0.0
        for j in 1:numfun
            frthdf = abs(
                2.0 * (1.0 - ratio) * rgnerr[j]
                    - (null_work[j, 7] + null_work[j, 8])
                    + ratio * (null_work[j, 5] + null_work[j, 6])
            )
            if abs(rgnerr[j]) + frthdf / 4.0 > abs(rgnerr[j])
                difsum += frthdf
            end
            for k in 1:4
                null_work[j, k] += w[k + 1, 2] * (null_work[j, 5] + null_work[j, 6]) +
                    w[k + 1, 3] * (null_work[j, 7] + null_work[j, 8])
            end
            basval[j] += w[1, 2] * (null_work[j, 5] + null_work[j, 6]) +
                w[1, 3] * (null_work[j, 7] + null_work[j, 8])
        end
        if difsum > difmax
            difmax = difsum
            divaxn = i
        end
        diff[i] = difsum
    end
    if direct >= 0
        direct_ref[] = Float64(divaxn)
    end

    # Remaining generators (positions 4..wtleng) via fully symmetric sums
    fulsms = view(null_work, :, 5)    # reuse as output buffer
    funvls = view(null_work, :, 6)    # reuse as per-point buffer
    for col in 4:wtleng
        _fully_symmetric_sum!(
            fulsms, center, hwidth,
            view(g, :, col), ndim, numfun, funsub, funvls
        )
        for j in 1:numfun
            basval[j] += w[1, col] * fulsms[j]
            for k in 1:4
                null_work[j, k] += w[k + 1, col] * fulsms[j]
            end
        end
    end

    # Error via null-rule linear combination (Genz-Malik)
    greate = 0.0
    for j in 1:numfun
        for i in 1:3
            search = 0.0
            for k in 1:wtleng
                search = max(
                    search,
                    abs(null_work[j, i + 1] + scales[i, k] * null_work[j, i]) * norms[i, k]
                )
            end
            null_work[j, i] = search
        end
        if errcof[1] * null_work[j, 1] <= null_work[j, 2] &&
                errcof[2] * null_work[j, 2] <= null_work[j, 3]
            rgnerr[j] = errcof[3] * null_work[j, 1]
        else
            rgnerr[j] = errcof[4] * max(null_work[j, 1], null_work[j, 2], null_work[j, 3])
        end
        rgnerr[j] *= rgnvol
        basval[j] *= rgnvol
        greate = max(greate, rgnerr[j])
    end
    greate_ref[] = greate

    # Sort ORDER[1..|direct|] by DIFF in descending order (for singular-region planning)
    if direct < 0
        neg_direct = -round(Int, direct)
        for i in 2:neg_direct
            cur = order[i]
            new_val = diff[cur]
            j = i - 1
            while j >= 1 && new_val > diff[order[j]]
                order[j + 1] = order[j]
                j -= 1
            end
            order[j + 1] = cur
        end
    end
    return nothing
end
