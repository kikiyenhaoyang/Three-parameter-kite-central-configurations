# Three-parameter-kite-central-configurations
This is the repository for the MATLAB code used in the computer-assisted proofs in the paper "Three-parameter kite central configurations in the planar five-body problem".


| ID | Main function | Verification scope |
|---|---|---|
| K0 | `verify_kite_quick_results` | Admissibility, positive masses, and the exact linear system at $(1/2,1/2,1/5)$. |
| N1 | `verify_kite_noninjectivity` | Two inequivalent nonsingular configurations with the same positive masses $(1/50,1,1/100,1,3/10)$, exact contraction certificates, and two local inverse branches. |
| K1 | `verify_kite_mass_orders` | All six strict axial mass orderings; exact linear systems, Cramer’s rule, and 60 Cartesian central-configuration residuals. |
| K2 | `verify_kite_regularity` | Jacobian bounds and a uniform contraction argument over an entire closed box, establishing a local arc where the three axial masses are equal. |
| S1 | `certify_singular_kite_ivt` | An intermediate value theorem certificate for singular kites using 324 subboxes, the positive-mass interval, and a local mass surface. |
| R1 | `verify_nondegeneracy` | A complete partition into 1,869 subboxes and sign conditions on $H,H_u,H_v$, establishing regularity of the equal-mass locus. |
| R2 | `verify_region_intersections` | A unique intersection with the equal-mass locus along each of four positive-mass line segments in regions $\mathcal A$, B, F, and G. |
| R3 | `verify_exclusion` | A complete partition into 4,021 subboxes and corner comparisons, excluding the relative closures of regions C, D, E, H, and I. |

K2, R1, R2, and R3 retain proofs covering entire intervals or regions. In particular, K2 includes interval bounds on partial derivatives and a uniform contraction argument, rather than relying solely on a Jacobian evaluation at one point. The original local and global scope remains unchanged: the package makes no additional claims about dynamical nondegeneracy, global connectedness, or the total number of global branches.

MATLAB with Symbolic Math Toolbox is required. After extracting the archive, set MATLAB’s **Current Folder** to `kite_matlab`. To run all bundled verifiers, execute:

```matlab
run_all
```

To run a single verifier, call its main function without the `.m` extension:

```matlab
verify_kite_noninjectivity
```

Alternatively, select a certificate by its ID through the runner:

```matlab
run_all('Only', {'N1'})
```

The runner saves individual logs and `summary.json` in the `results` folder.
