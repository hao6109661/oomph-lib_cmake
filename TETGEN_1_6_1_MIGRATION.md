# TetGen 1.6.1 migration and validation

## Scope

Repository: `hao6109661/oomph-lib_cmake`.
Branch: `fixed_size_vector_with_new_tetgen`.
Starting commit: `b0c418635e8245e046cad1ef86e7b1759d7389a1`.
The latest substantive user commit, `31261c7`, was explicitly treated as
unfinished work, rather than an independently validated baseline.

## Upstream version

The sources come from the official TetGen mirror, tag `v1.6.1`, commit
`205515773c72c0f2d0d8a16db200a32d748345eb`. The release was published on
17 September 2026. The primary development repository is
<https://codeberg.org/TetGen/TetGen>; its mirror is
<https://github.com/TetGen/TetGen>.

The three source files are copied unchanged. Their formatting and existing
operators belong to upstream TetGen; they have not been reformatted as
oomph-lib code. The upstream `LICENSE` is included alongside them.

SHA-256 checksums:

```
e6725da7cc9b2d9c36502c30f9161f7e0ede13f9e4e1e1246173089fda60537b  tetgen.cxx
d5a8e94d93e3423ee68439ff4e3c6c0517839e263433ac8de2cdaee322468af8  tetgen.h
8f1decc8526f572d9f1167827629fd256f519d7f02918865909ad907ed4ca4b3  predicates.cxx
```

## Changes and reasons

1. Replace the three source files in `external_src/oomph_tetgen_1.6` and
   add upstream's `LICENSE`. Preserve the directory and target names so
   existing include paths and consumers continue to work.
2. Declare `OOMPH_USE_NEW_TETGEN` as a documented CMake option, enabled by
   default. Setting it to `OFF` selects the bundled legacy TetGen 1.4.3.
3. Define `TETLIBRARY` with `PUBLIC` visibility on the newer TetGen target.
   The upstream header leaves this macro undefined. It controls both the
   `tetrahedralize(char*, ...)` declaration and inline error handling, so
   library sources and consuming translation units need consistent values.
   The conservative predicate-compilation setting also matches AppleClang,
   rather than only GNU and upstream Clang.
4. Reset `OOMPH_HAS_NEW_TETGEN` before selecting the implementation. This
   prevents a cached newer-version flag from surviving a switch to legacy.
5. Extend the existing `deep_copy_of_tetgenio()` to copy `Subdomain_Facets`,
   `numberofrefineelems`, `refine_elem_list` and `refine_elem_vol_list`.
   These are additional region/refinement inputs in 1.6.1. The vector is
   copied by value and the arrays receive independently owned storage.
6. Explicitly enable the newer implementation in both Ubuntu and macOS
   self-test workflows using `--oomph-OOMPH_USE_NEW_TETGEN=ON`.

No ternary operators are added to the oomph-lib adapter. No unrelated
repository or branch is within the modification scope.

## Review of the unfinished commit

The following areas were checked against the upstream interface:

- Boundary membership must be reconstructed from positive face markers.
  Point markers alone cannot describe all boundaries meeting at an edge or
  corner. The existing global edge-boundary lookup is therefore retained.
- Zero-based TetGen connectivity must be converted before using the
  scaffold's internal one-based indexing. The existing conversion is retained.
- The scaffold consumes four-node tetrahedra. oomph-lib constructs its own
  higher-order nodes; passing TetGen `-o2` output is unsupported.
- `carveholes()` assigns automatic region attributes above the largest
  explicitly processed attribute. The existing remapping to oomph-lib's
  default region zero remains necessary for the supplied valid region seeds.
- Existing fixes to neighbour copying, three-entry segment constraints and
  optional marker arrays are retained.
- Disabling optimisation for `predicates.cxx` is consistent with the
  upstream build, which also specifies `-O0` for geometric predicates.

## Important behavioural limits

In library mode, TetGen reports failures by throwing integer error codes.
In 1.6.1, `-YY` throws error 200 if boundary Steiner points remain after
suppression. It cannot make every prescribed PLC meshable without changing
its boundary. Successful linear meshing also does not guarantee positive
Jacobians after oomph-lib snaps higher-order nodes onto curved geometry.
Existing inversion checks must remain active.

## Validation

Local validation results are recorded below. GitHub Actions validation is
pending publication of this commit; no claim of a completed CI matrix is made.

## Geometry compatibility repairs

The numerical tests exposed two issues that cannot be fixed by replacing the
TetGen header alone:

1. The curved demo supplied 1,140 non-planar quadrilaterals. For example,
   boundary marker 11302 had a fourth-vertex distance of about 0.001283 from
   the plane of its first three vertices. TetGen 1.6.1 reported intersecting
   projected segments and skipped parts of these facets. `build_tetgenio()`
   now splits a warped quadrilateral along diagonal 0--2 into two planar
   triangles with the same boundary marker. A dimensionless scalar-triple-
   product test distinguishes warped quads from planar ones. Other polygons
   retain their existing treatment. The original oomph-lib surface objects
   are preserved, including their role in geometric snapping. The legacy
   implementation retains its original PLC representation.
2. Boundary projection can invert a small number of higher-order tetrahedra
   even when the input linear tetrahedra have positive Jacobians. The generic
   snapping routine now attempts a local interior-node repair when this is
   detected. It selects only non-boundary nodes of initially inverted cells
   and includes all incident elements in its objective. Boundary positions
   are never changed. The penalty is the sum of squared deficits below a
   normalised Jacobian of 1e-4, sampled at the element's integration knots.
   Derivatives are computed by central differences and a backtracking line
   search accepts only decreases in the patch penalty. The repair is bounded
   by 100 sweeps and 20 backtracking trials. Regardless of its result, the
   existing full-mesh Jacobian checks and failure output run afterwards.
   A mesh that still has inverted elements is rejected.

## Mesh-size controls and the revised acceptance criteria

The user accepts reasonable changes in the discrete solution caused by different
meshes. Preventing inverted elements and avoiding unnecessary growth in degrees
of freedom therefore take priority over matching the old regression snapshot.
No additional adaptation is performed just to match that snapshot.

For initial construction the newer switch string is
`pq2.0/0/150/172aAa<volume>`, followed by the existing `YY` only when the
prescribed boundaries prohibit splitting. The radius-edge target remains 2.0.
The minimum dihedral-angle target is left at the upstream default of zero.
The segment and facet encroachment angles are 150 and 172 degrees: these are
exactly the relaxed values that upstream 1.6.1 itself selects for refinement
in `tetgenmesh::delaunayrefinement()` in `tetgen.cxx`. They reduce unnecessary
boundary encroachment refinement while still allowing boundary splitting
where the original surface permits it. There is no added `D` switch.

For refinement the newer switch string is `Vq2.0ras/0O//0`. This retains the
per-element volume constraints and radius-edge target, but disables
post-refinement smoothing and mesh-improvement iterations. Those operations
can enlarge cells after volume-controlled refinement and previously caused
repeated size-transfer cycles. The existing size-accommodation criterion is
unchanged. It converged at iteration 1 in the final curved run.
The legacy strings remain `pq2.0aAa<volume>` and `Vqra`.

The initial and refinement controls were selected from geometric behaviour
and upstream source, without tuning against the solution norm. The demo's
geometry, volume request, error targets and number of adaptation steps are
unchanged. No unknowns are added by the local node repair.

## Investigation and discarded approaches

1. Rebuilt the legacy implementation against the branch's latest substantive
   user changes. Its adaptive-mesh CTest passed in 225.45 seconds, so the
   unfinished commit was reviewed rather than reverted wholesale.
2. Built official 1.6.1 independently: the cube-with-hole check produced 565
   nondegenerate tetrahedra and volume 26, equal to the analytic value.
3. Verified the newer header, library macro and independent ownership of the
   additional region/refinement data in the adapter.
4. Reproduced 60 inversions with a direct source replacement. Splitting the
   non-planar quads removed the projected-segment warnings but did not alone
   fix all inversions.
5. Checked repeated geometric projections of shared nodes; their maximum
   disagreement was about 3.33e-16. Inconsistent boundary projection was not
   the cause in this example.
6. Tried stronger quality controls and a harmonic-displacement prototype.
   Harmonic displacement alone left four inversions. The smaller local
   interior-node repair was retained instead of adding a global solve.
7. An earlier passing candidate used `D` and a five-degree target. Although
   its final norm passed the old curved snapshot, it required 185,143
   unknowns versus 76,459 with legacy. It is superseded and is not the code
   supplied in the revised patch.
8. Removing only those extra quality controls still required 100,952 initial
   unknowns versus legacy's 51,139. This run was deliberately interrupted
   during the first solve once the excessive initial size was established;
   it is not counted as a completed numerical validation.
9. Forcing `Y` even on splittable boundaries reduced the initial linear mesh
   to 28,751 tetrahedra but left 19 inversions after local repair. It was
   rejected. The final implementation does not force this restriction.
10. Applied upstream's relaxed encroachment angles to initial construction.
    The final curved and planar runs both completed normally, with all
    post-projection Jacobian checks passing and fewer unknowns than legacy.

## Final local results

Environment: Ubuntu 24.04, GCC/GFortran 13.3, CMake 3.28.3, Ninja 1.11.1,
CGAL 6.0.1, Boost 1.83, SuperLU 7.0.1, OpenBLAS 0.3.26. Serial library,
with `OPENBLAS_NUM_THREADS=1` for numerical runs.

| Quantity | Legacy 1.4.3 | Final 1.6.1 |
| --- | ---: | ---: |
| Curved initial tetrahedra | 44,325 | 32,620 |
| Curved initial unknowns | 51,139 | 36,366 |
| Curved final tetrahedra | 65,264 | 51,226 |
| Curved final unknowns | 76,459 | 58,324 |
| Curved first solved norm | 2.20633 | 2.203 |
| Curved final solved norm | 2.20054 | 2.20008 |
| Planar initial tetrahedra | 905 | 171 |
| Planar initial unknowns | 912 | 165 |
| Planar final tetrahedra | 8,784 | 3,285 |
| Planar final unknowns | 9,999 | 3,640 |
| Planar final solved norm | 1.70436 | 1.74548 |

The final curved unknown count is 23.72% below legacy. Its norm differs from
the actual legacy run by 0.02090%. The original curved regression comparison
passes at the unchanged 0.1% threshold. All curved construction and adaptation
stages report no inverted elements after snapping and local repair.

The planar final norm differs from legacy by 2.41264%, using the legacy norm
as denominator. `fpdiff.py` reports 2.35580% because it uses the larger value
as denominator, and exits with status 2 at the original 0.1% threshold. This
is the only failed entry in that comparison; its volume entry passes.
The regional volumes remain equal to the analytic values: 0.216, 17.568 and
17.784 for interior, exterior and total volume, with relative errors around
4e-12 percent. All planar Jacobian checks pass. This is a materially coarser
planar mesh, and the norm difference is reported rather than concealed.
A close norm alone is not proof of continuum accuracy or a small pointwise
solution error; this report does not claim either.

The original demo sources, validation script, CMake configuration and reference
files are retained. In particular, the earlier proposed version-specific planar
snapshot and its selection logic have been removed. There is no increase in
regression tolerance, skipped numerical check, or extra adaptation to force a
snapshot match. The standard adaptive CTest will therefore still report the
planar numerical comparison failure with this configuration.

| Other check | Result |
| --- | --- |
| Main-library build and installation | Passed before final header-only mesh-control adjustment |
| Final curved and planar release executables | Rebuilt; both exited with status 0 |
| Legacy complete adaptive CTest | Passed, 225.45 seconds |
| New installed-library file-reader CTest | Passed, 24.10 seconds |
| Native version and adapter deep-copy test | Passed |
| Native TetGen consistency and independent volume check | Passed |
| Formatting of final two mesh-control headers | Passed clang-format 18 dry-run |

The final geometry controls were verified by directly running both numerical
executables and the original `fpdiff.py` comparisons. The entire newer CTest
matrix has not been rerun. Gmsh-specific comparisons were not run locally:
the script requires Gmsh 3.0.6, which was not installed. The file-reader and
copy tests exercise unchanged compiled adapter code; the final changes are in
the mesh-construction headers instantiated by the rebuilt demo executables.

## Scope of the inversion evidence

No inverted element was accepted in either completed final demo. The original
full-mesh failure check remains active after the bounded local repair. The
check samples the existing integration knots; it is not a proof that the
higher-order Jacobian stays positive everywhere inside an element. A geometry
that cannot be repaired with fixed boundary positions still raises an error.
The discarded forced-Y experiment is direct evidence of this limitation.

The disk-in-tank circular, elliptical and rounded-square parameter sweeps have
not been run with this candidate. Their successful meshing, including protected
interface-node counts, must be verified before claiming that those production
cases or every possible curved geometry are resolved. The supplied evidence
supports the two oomph-lib demos and the stated interface checks only.

## Reproduction

Start with fresh main-library and demo build directories to avoid mixing
incompatible TetGen headers and libraries. For an existing third-party-library
installation generated by this repository, a typical sequence is:

```bash
cmake -S . -B build-tetgen161 -G Ninja \
  $(cat external_distributions/install/cmake_flags_for_oomph_lib.txt) \
  -DOOMPH_USE_NEW_TETGEN=ON \
  -DCMAKE_INSTALL_PREFIX="$PWD/install-tetgen161"
cmake --build build-tetgen161 --parallel 6
cmake --install build-tetgen161
cmake -S demo_drivers -B demo_drivers/build-tetgen161 -G Ninja \
  -Doomphlib_ROOT="$PWD/install-tetgen161"
OPENBLAS_NUM_THREADS=1 ctest --test-dir demo_drivers/build-tetgen161 \
  -R 'meshing.(adaptive_tet_meshes|mesh_from_tetgen)' --output-on-failure -j 1
```

Use the paths for the actual third-party installation if they differ. For the
complete test suite, remove the `-R` selection. Both updated self-test workflows
explicitly enable the newer TetGen. Use the usual `oomph_add_executable()`
wrapper for stand-alone applications so exported oomph-lib compile definitions
are applied consistently.


## Publication and CI status

The user explicitly authorised uploading this candidate to the branch
`fixed_size_vector_with_new_tetgen` on 9 November 2026, replacing the earlier
read-only restriction for this operation. This report accompanies that upload.
No other repository or branch is within its scope. The earlier rejected write
attempt made no remote changes.

The remote CI results are not yet available at the time this report is committed.
The known planar regression difference remains visible, as described above;
the full CI matrix must not be described as passing without checking its results.
