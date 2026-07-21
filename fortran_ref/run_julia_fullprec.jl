# Resultats Julia (Decuhr.jl) en PLEINE precision, memes cas/parametres que le
# driver Fortran (driver.f90), pour une comparaison bit-a-bit avec la reference.
import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))   # le package Decuhr.jl (dossier parent)
using Decuhr
using Printf

const B = 1_000_000
const W = 50_000
const E = 20

cases = [
    ("C1_xy_vertex", (x->(x[1]*x[2])^(-0.5)),            [0.0,0.0],     [1.0,1.0],     2, -0.5,      0, 4.0),
    ("C2_radial",    (x->1/sqrt(x[1]^2+x[2]^2)),         [0.0,0.0],     [1.0,1.0],     2, -1.0,      0, 2*log(1+sqrt(2.0))),
    ("C3_smooth",    (x->sin(x[1])*cos(x[2])),           [0.0,0.0],     [pi/2,pi/2],   1,  0.0,      0, 1.0),
    ("C4_3d_vertex", (x->(x[1]*x[2]*x[3])^(-1/3)),       [0.0,0.0,0.0], [1.0,1.0,1.0], 3, -1/3,     0, 27/8),
    ("C5_poly",      (x->x[1]^2+x[2]^2),                 [0.0,0.0],     [1.0,1.0],     1,  0.0,      0, 2/3),
    ("C6_log",       (x->-log(x[1]*x[2])),               [0.0,0.0],     [1.0,1.0],     2,  0.0,      1, 2.0),
]

println("case          result                exact                relerr      ifail neval")
for (name, f, lb, ub, singul, alpha, logf, exact) in cases
    ndim = length(lb)
    funsub = (x, fv) -> (fv[1] = f(x); nothing)
    r, e, neval, ifail = Decuhr._decuhr_driver(
        ndim, 1, lb, ub, 0, B, funsub, Float64, singul, alpha, logf, 1.0e-8, 1.0e-6, 0, W, E)
    relerr = abs(r[1]-exact)/abs(exact)
    @printf("%-12s  %.15E  %.13E  %.3E  %3d  %9d\n", name, r[1], exact, relerr, ifail, neval)
end
