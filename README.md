# Three-parameter-kite-central-configurations

This is the repository for the MATLAB code used in the computer-assisted proofs in the paper *"Three-parameter kite central configurations in the planar five-body problem"*.

## **Structure**

- Each `.m` file is implemented as a standard MATLAB function and begins by checking that the Symbolic Math Toolbox is available.
- The main verification logic appears first, followed by any file-local helper functions required by that verifier. Each verifier can be run independently. R1 and R3 additionally require their corresponding JSON certificate files.
- The `require` helper uses `isAlways(...,'Unknown','false')`, so a condition is accepted only when MATLAB can prove it to be true.
- Linear systems and algebraic residual checks follow the same symbolic workflow used in the reference examples, including `syms`, `subs`, `simplifyFraction`, `P0\b0`, and `equalZero`. Determinants are retained as exact scalar expressions.
- Global-domain verification stores symbolic interval endpoints and derivatives in ordinary MATLAB structs. Exact interval arithmetic is implemented by local helper functions at the end of each file.
- All mathematical inequalities are checked using exact rational arithmetic. Decimal approximations are used for display only. Strict interval bounds are rounded outward using symbolic `floor` and `ceil`.

The package contains ten `.m` files, including eight verifiers, `run_all.m`, and `selftest.m`. `run_all` supports both verifiers that return structured summaries and verifiers with no return value. Individual verifiers can still be invoked directly in the same way as before.

| ID | Main function | Verification scope |
|---|---|---|
| K0 | `verify_kite_quick_results` | Admissibility, positive masses, and the exact linear system at $`(\frac12,\frac12,\frac15)`$. |
| N1 | `verify_kite_noninjectivity` | Two inequivalent nonsingular configurations with the same positive masses $`(\frac{1}{50},1,\frac{1}{100},1,\frac{3}{10})`$, exact contraction certificates, and two local inverse branches. |
| K1 | `verify_kite_mass_orders` | All six strict axial mass orderings; exact linear systems, Cramer’s rule, and 60 Cartesian central-configuration residuals. |
| K2 | `verify_kite_regularity` | Jacobian bounds and a uniform contraction argument over an entire closed box, establishing a local arc where the three axial masses are equal. |
| S1 | `certify_singular_kite_ivt` | An intermediate value theorem certificate for singular kites using 324 subboxes, the positive-mass interval, and a local mass surface. |
| R1 | `verify_nondegeneracy` | A complete partition into 1,869 subboxes and sign conditions on $`H`$, $`H_u`$, and $`H_v`$, establishing regularity of the equal-mass locus. |
| R2 | `verify_region_intersections` | A unique intersection with the equal-mass locus along each of four positive-mass line segments in regions $`\mathcal{A}`$, $`\mathcal{B}`$, $`\mathcal{F}`$, and $`\mathcal{G}`$. |
| R3 | `verify_exclusion` | A complete partition into 4,021 subboxes and corner comparisons, excluding the relative closures of regions $`\mathcal{C}`$, $`\mathcal{D}`$, $`\mathcal{E}`$, $`\mathcal{H}`$, and $`\mathcal{I}`$. |

K2, R1, R2, and R3 retain proofs covering entire intervals or regions. In particular, K2 includes interval bounds on partial derivatives and a uniform contraction argument, rather than relying solely on a Jacobian evaluation at one point. The original local and global scope remains unchanged: the package makes no additional claims about dynamical nondegeneracy, global connectedness, or the total number of global branches.

MATLAB with Symbolic Math Toolbox is required. After downloading the project, extract the archive into a new folder and set MATLAB’s **Current Folder** to that folder (the folder containing `run_all.m`). To run all bundled verifiers, execute:

## **Running Verifications**

```matlab
% Run selected verifications
run_all('Only', {'K0'});
run_all('Only', {'N1'});
run_all('Only', {'K0','S1'});
report = run_all('Only', {'K2','R3'});

% List available items, check package files, or select an output directory
run_all('List', true)
run_all('CheckOnly', true)
run_all('OutputDir', fullfile(pwd,'my_results'))

% Use explicit certificate files
verify_nondegeneracy('path/to/nondegeneracy_certificate.json')
verify_exclusion('path/to/exclusion_certificate.json')

% Run quick tests for exact rational arithmetic,
% interval arithmetic, and invalid-input handling
selftest

% Recompute all eight certificates and compare the five available
% exact reference summaries: K1, K2, R1, R2, and R3
selftest('Full', true)
```

The runner saves individual logs and `summary.json` in the `results` folder.


## **Logging and Self-Tests**

`selftest('Full', true)` invokes `run_all` internally, so there is no need to run `run_all` again afterward.

By default, per-verifier logs are written to `results/`, and the combined machine-readable summary is written to:

```text
results/summary.json
```

Symbolic verification involves computations over large exact rational expressions. As a result, full-domain verification may take significantly longer than the single-point examples provided in the reference material.

K0 and S1 write their complete output to:

```text
results/K0.log
results/S1.log
```

These two verifiers do not return MATLAB values. Their corresponding `details` fields in `summary.json` are therefore empty arrays.

If an assertion fails, the corresponding verification is recorded as `FAIL`. `run_all` continues through the requested items and raises an error after they have completed. Structured summaries produced by the other verifiers are preserved unchanged.

N1 writes its textual output to:

```text
results/N1.log
```

Its full interval data, contraction bounds, and similarity invariants are stored in the N1 `details` entry in `summary.json`.

## **Environment Requirements**

The package targets MATLAB R2020b or later with the Symbolic Math Toolbox installed and available.

This compatibility range has not yet been validated across actual MATLAB installations.

All mathematical computations use symbolic (`sym`) arithmetic. 

## **Certificate Data**

The two JSON certificate files are:

```text
nondegeneracy_certificate.json
exclusion_certificate.json
```

Their byte contents are identical to those of the original attachments.

Then use `run_all` for all eight bundled checks, or `selftest('Full', true)`
for the arithmetic self-tests, all eight checks, and the five available
reference comparisons. `run_all` captures each verifier's output and writes
its log after that verifier returns, so a long `RUN ...` interval alone does
not establish that MATLAB has hung.

The reference time of implementation is as follows:
```matlab
RUN K0: Exact local positive-mass witness and mass Jacobian
PASS K0 (1.547 s)
RUN N1: Two inequivalent nonsingular kites with identical positive masses
PASS N1 (27.376 s)
RUN K1: Six strict axial mass orders
PASS K1 (3.155 s)
RUN K2: Local kite regularity and triple-equality arc
PASS K2 (41.025 s)
RUN S1: Singular kite IVT and local positive-mass surface
PASS S1 (3284.966 s)
RUN R1: Rhomboidal equal-mass locus regularity
PASS R1 (4767.922 s)
RUN R2: Positive-mass intersections in regions A, B, F and G
PASS R2 (27.947 s)
RUN R3: Exclusion of relative closures C, D, E, H and I
PASS R3 (3909.443 s)
```
