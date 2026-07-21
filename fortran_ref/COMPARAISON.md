# Comparaison directe au Fortran de référence (Espelid & Genz, 1994)

Validation de ce portage Julia (`Decuhr.jl`) par comparaison **bit-à-bit** avec
la subroutine Fortran `DECUHR` d'origine, compilée localement (gfortran/MinGW),
sur **6 cas** (définis à l'identique dans `driver.f90` et
`run_julia_fullprec.jl`) avec **exactement les mêmes paramètres** :

```text
epsabs=1e-8  epsrel=1e-6  maxpts(budget)=1e6  wrksub=50000  emax=20
minpts=0  key=0  restar=0  numfun=1
```

## Contenu du dossier

| Élément | Rôle |
|---|---|
| `src_fortran/` | **Copie** des 14 sources Fortran de référence (`*.f`) |
| `Makefile` | Compile les 14 sources + le driver → `decuhr_cmp.exe` |
| `driver.f90` | Programme appelant `DECUHR` sur les 6 cas (C1–C6) |
| `run_julia_fullprec.jl` | Mêmes cas via le cœur Julia `Decuhr._decuhr_driver`, pleine précision |
| `logs/run_fortran.log` | Sortie du binaire Fortran |
| `logs/run_julia.log` | Sortie Julia pleine précision |

## Résultats — Fortran vs Julia (pleine précision)

| Cas | Fortran (référence) | Julia (`Decuhr.jl`) | Écart absolu | `IFAIL` (F / J) | `NEVAL` (F / J) |
|---|---|---|---|:---:|:---:|
| C1 `(x·y)^(-1/2)`     | `4.000060778932444` | `4.000060778932444` | **0** (identique) | 0 / 0 | 110045 / 110045 |
| C2 `1/√(x²+y²)`       | `1.762747166075242` | `1.762747166075242` | **0** | 0 / 0 | 715 / 715 |
| C3 `sin·cos`          | `1.000000000000012` | `1.000000000000012` | **0** | 0 / 0 | 195 / 195 |
| C4 `(x·y·z)^(-1/3)`   | `3.375030993803787` | `3.375030993936098` | ≈ 1.3·10⁻¹⁰ | **1 / 1** | 999871 / 999871 |
| C5 `x²+y²`            | `0.6666666666666667`| `0.6666666666666667`| **0** | 0 / 0 | 195 / 195 |
| C6 `-log(x·y)`        | `2.000000028790773` | `2.000000028790773` | **0** | 0 / 0 | 9035 / 9035 |

## Conclusions

1. **Le portage Julia reproduit le Fortran à l'identique.** 5 cas sur 6 sont
   **bit-à-bit identiques** (16 chiffres significatifs). Les `IFAIL` et `NEVAL`
   coïncident sur **les 6 cas**. C'est la validation la plus forte possible
   (l'ancienne réserve « pas de comparaison binaire au Fortran » est **levée**).

2. **C4 : écart de ≈ 10⁻¹⁰** (relatif ≈ 4·10⁻¹¹), au niveau du bruit de
   **réassociation des opérations flottantes** (ordre des sommations dans la
   règle 3-D / l'extrapolation, cumulé sur ~10⁶ évaluations à budget épuisé) —
   et **non** une divergence algorithmique : même `IFAIL=1`, même `NEVAL`,
   même erreur relative à l'exact (9.18·10⁻⁶).

3. **Réévaluation de la remarque « `retcode = MaxIters` malgré un résultat
   correct ».** Le Fortran de référence renvoie **exactement le même `IFAIL`**
   que Julia :
   - C1 : `IFAIL=0` (Succès) dès que le budget le permet (110045 < 10⁶). Le
     `MaxIters` observé via l'interface Integrals.jl venait uniquement du défaut
     `maxiters = 100000 < 110045` — augmenter `maxiters` suffit.
   - C4 : `IFAIL=1` (MaxIters) **côté Fortran aussi**, avec le même résultat
     correct non certifié.

   ⇒ L'estimateur d'erreur conservatif (`EXTERR = 10·|T(STEPS)−T(STEPS−1)|`)
   est une **caractéristique intrinsèque de l'algorithme DECUHR**, fidèlement
   reproduite. **On ne peut pas « faire mieux » sans s'écarter de la référence**
   (ce que le périmètre interdit). La réponse appropriée — déjà en place —
   est : documenter le comportement, exposer `sol.stats` (neval/ifail), et
   conseiller de se fier à `sol.u` / `sol.resid` plutôt qu'au seul `retcode`.

## Reproduire

Nécessite **gfortran** (MinGW) et **make** :

```bash
cd fortran_ref
make                                          # compile -> decuhr_cmp.exe (cf. Makefile)
./decuhr_cmp.exe            > logs/run_fortran.log   # resultats Fortran
julia run_julia_fullprec.jl > logs/run_julia.log     # resultats Julia (pleine precision)
```

Si `gfortran` n'est pas dans le `PATH` : `make FC=/chemin/vers/gfortran`.
La cible `make run` compile puis exécute ; `make clean` supprime l'exécutable.

> Les sources Fortran de `src_fortran/` portent le copyright Espelid & Genz
> reproduit en tête de chaque fichier (leur licence autorise explicitement la
> copie tant que cette notice est conservée — cf. `NOTICE` à la racine du
> dépôt). Ce dossier fait partie du dépôt `Decuhr.jl` comme banc de validation
> reproductible ; il n'est pas exécuté par la CI (`Pkg.test`).
