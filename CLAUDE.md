# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

TDSolve is a C++17 code for 1D / radially / spherically symmetric free-boundary tumor-growth PDE systems (cell species advected by a velocity field inside a moving radius `R(t)`, coupled to diffusing chemical species). It builds two executables:

- `TumorDynRun` (`main/TumorDyn.cc`) — a single simulation from one input deck.
- `ParAnalRun` (`main/ParAnal.cc`) — OpenMP-parallel parametric sweeps over many input decks, with optional least-squares response-surface fitting.

Version control is **Mercurial** (`hg`), not git. `build/`, the executables, and LaTeX by-products are in `.hgignore`.

## Build

`readme.txt` says "type make" — that is outdated; there is no top-level Makefile. Use CMake through the wrapper script:

```
sh build_cmake.sh          # configures into build/, make -j8, symlinks ./TumorDynRun and ./ParAnalRun
cd build && make -j8       # incremental rebuild after the first configure
```

- Requires BLAS, LAPACK, OpenMP (optional), and **Trilinos** (for Sacado automatic differentiation). The script locates Trilinos at `$CSUPPROOT/$TDSUPPNAME/trilinos`; `build_cmake_repo.sh` is the same but uses `$GENSUPPROOT/$GENSUPPNAME/trilinos` (shared `/repo` install) and does not create the `./ParAnalRun` symlink. These env vars come from the user's shell profile.
- `CMakeLists.txt` globs `src/*.cc`, so new source files need a CMake re-configure (rerun the script), not just `make`.
- Build type is hard-coded to Release. There is no lint target.
- Precision is `DType` (`include/Data_Type.h`): `double` unless `FLOAT_ALG` is defined. Use `DType`, not `double`, in solver code.

## Running

```
cd runs/test && ../../TumorDynRun -i input.hcl     # -h prints option help
```

Outputs are `outputDataFile` (`sol.dat`) and `radiusDataFile` (`rad.dat`) in the working directory. MATLAB scripts in `tools/` read/plot them (`*_tumordyn.m` for the FDM/FVM frameworks, `*_sacado.m` for the Sacado framework).

`ParAnalRun -i base.hcl -p input.par -n {PDGF|OPN|Potassium} [-l list.par] [-o fit.out -m LS2D] [-Sk 1 -Sd data]` — run with no/bad flags for the full usage text. `-p` gives a tensor-product sweep, `-l` an explicit list of parameter tuples. Results are appended to the data file (`-d`, default `regression.dat`). Note `-Sk 1` deletes and recreates the `-Sd` directory.

## Tests / benchmarks

There is no unit-test framework. Regression checking is by comparison against stored reference output:

- `benchmark/Sacado/{model,pdgf}/<case>/` each hold `input.hcl`, `run.sh`, and a `benchmark/` directory with reference `sol.dat`, `rad.dat`, `report`. `run.sh` needs `$TDS` set to the repo root; it runs the case into `data_r<hg revision>/`, which you then diff against `benchmark/`. Run one case with `cd benchmark/Sacado/pdgf/pdgf_case_be && sh run.sh`.
- `tests/Sacado*` are standalone Sacado/Trilinos smoke programs with their own CMake + `build_cmake.sh` (paths there still point at an old gcc-9.3.0/Trilinos-12.7 toolchain). `tests/potassium_ode/` is a MATLAB reference implementation for the potassium ODE model.

## Architecture

### Input decks (`TDInput`)

`.hcl` files are `name = value` lines (delimiters space, `=`, `:`, tab; comments start with `//` or `#`). `TDInput` registers every recognized key in three maps in `src/TDInput.cc` — `tdParameters` (DType), `tdOptions` (int), `tdToggles` (bool) — each pointing into the `TDCoef` / `TDMesh` / `TDSetup` / `TDLin` structs. **To add an input key, add the struct field, its default in the constructor, and a map entry.** `ParAnalRun` clones a base deck per run through the copy constructor, which copies `TDCoef`/`TDSetup`/`TDLin` wholesale but `TDMesh` field by field — new `TDMesh` members must be added there by hand. Input-file key names often differ from struct field names (e.g. `framework` → `tdSetup.frame`, `time_method` → `tmeth`, `m` → `alpha_a`; see the coefficient table at the top of `include/TDInput.h`). Unknown keys only produce a warning on stderr and the run continues with defaults, so a misspelled key does not stop a run.

A line consisting of a problem name (`model`, `pdgf`, `opn`, `oncolyticOld`, `potassium`, `oncolytic`, `onc_model`) opens a nested block of per-variable parameters terminated by `end`, stored in `tdCoef.probVarParVals[problem][variable][parameter]`.

### Solver hierarchy (`TDSolve`)

`TDSolve` (abstract, `include/TDSolve.h`) owns the shared state — solution arrays indexed by time level/stage (`data_*[nt]`), mesh, time-step control, radius update (including Radau variants), and linear-solver setup. The `framework` option selects the subclass:

| `framework` | Class | Notes |
|---|---|---|
| 0 | `TDSolveFDM` | original finite difference |
| 1 | `TDSolveFVM` | finite volume |
| 2 | `TDSolveFVMMod` (extends FVM) | modified FVM; MUSCL, Radau RK |
| 3 (and any other value) | `TDSolveSacado` | current, actively developed path |

Frameworks 0–2 hard-code the PDGF/OPN variables (`data_g`, `data_h`, `data_n`, `data_a`, `data_p1`, …) and are legacy. **`TDSolveSacado`** is generic: variables are stored in name-keyed maps, Jacobians come from `Sacado::Fad::DFad`, and each model is plugged in through member-function-pointer maps keyed by `prob_id` (−1 model, 0 pdgf, 1 opn, 2 oncolyticOld, 3 potassium, 4 oncolytic, 5 onc_model). Adding a model means writing `setup_var_*`, `setup_ic_*`, `setup_srce_*` (plus optional velocity/aux/extra-flux hooks), registering them in the maps in the `TDSolveSacado` constructor/`setup`, and adding the name to `probNames` in `TDInput`.

### Time-stepping lives in the drivers, not the solver

The time-integration loop is written out in `main/TumorDyn.cc` in terms of `TDSolve` primitives (`update_sol(fr,to,…)`, `copy_data`, `average_sol`, `update_rad_radau`, `compute_rhs`, `check_conv_*`, `reduce_dtn`). There are two branches: a short "classical" loop for the Sacado framework (`time_method` 0/1 = explicit/implicit Euler with Newton and automatic dt halving) and a longer legacy branch for frameworks 0–2 (TVD RK2, DIRK2, 2-stage Radau with direct or predictor/multi-corrector solves, selected by `time_method`, `vel_couple`, `use_pmc`). Convention: the accepted solution at t_{n+1} must end up in `data_*[1]` and is then copied to `[0]`. `manual/dev/radau.txt` documents which arrays each step of the Radau PMC loop reads and writes — consult it before touching that path.

`TDFactory` (`src/TDSolve_API.cc`) **duplicates** this loop for `ParAnalRun` (same `framework` switch as `TumorDyn.cc`; `run_pdgf` / `run_opn` / `run_potassium` each return model-specific output quantities via `get_gperc`, `get_total`, `get_argmax`, …). A change to stepping logic in `TumorDyn.cc` generally needs the matching change there.

### Supporting numerics (own implementations, no external linear algebra beyond BLAS/LAPACK)

- `HighDim_Alg` — `HVEC`, `HMAT`, tridiagonal/triangular matrices, tensors; `LowDim_Alg` — fixed-size `VEC3/4`, `MAT3/4`.
- `Alg_Solve` — `LinVec` / `LinMat` (full, tridiagonal, sparse) and `LinSolve` (direct, stationary iterative, Krylov, ILU preconditioning), configured by the `lsys*` input keys and the enums in `Comput_Choices.h`.
- `Alg_Limiters` — slope limiters; `Quad_Coef`, `RK_Coef` — quadrature rules and Butcher tables.
- `Fun_Estimator` + `LeastSquares_Solve` — polynomial least-squares fitting used by `ParAnalRun -o` (`-m LSnD` total degree, `LSnT` per-variable degree).

## Documentation

`manual/doc/manual.tex` is the user manual; `manual/numerics/*.tex` (and the PDFs in `manual/`) derive the numerical methods for the PDGF and oncolytic models; `manual/demo/input_sample.hcl` is an annotated-by-example input deck.
