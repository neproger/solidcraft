# Analytical miter for planar Sweep: research charter

Status: active research. The AOIM Line-Arc hypothesis is already validated in
`research/aoim-line-arc`; it is not to be re-derived.
Created: 2026-11-10
Source research: `https://github.com/neproger/solidcore/tree/research/aoim-line-arc` at `9920e32`
Related code: `crates/solidcore/src/primitive/primitive_creation/p2/planar_path.rs`
Related reference: `docs/procedural-sweep.md`
Stage 0 result: `plan/analytical-miter-stage0.md`
Development journal: `plan/analytical-miter-journal.md`

## Why this file exists

The first attempt at this research drifted: the goal was an analytical miter for planar
Sweep, but a kernel limitation led to Boolean-union of nominal closed segment bodies. Union
worked, became the main effort, and was treated as if it answered the original question. It
did not.

This charter fixes the original goal, the expected result, the boundaries, the success
criteria and the changes that need sign-off. It is the reference every stage is measured
against. If the work no longer matches it, stop and agree on the change before writing code.

## Existing research: do not re-derive it

The mathematical work exists in the standalone branch `research/aoim-line-arc` (commit `9920e32`).
It is independent of this checkout and of the kernel crates. Read it before changing the plan:

- `README.md`: AOIM algebra, domain facts, hybrid and adaptive strategies, 2D and 3D checks.
- `BREP_PLAN.md`: direct shared-topology construction of the final shell, no Boolean and no
  internal seam face.
- `BREP_REPORT_RU.md`: the technical report with the blocked kernel paths.
- `src/`, `probes/brep_probe`: the implementation and the public-API probes.

Established and not open for re-checking: the AOIM offset-intersection curve and its existence
limit `d_min = -R cos^2(alpha/2)`; the finite-span coverage table; the external-bisector tangent;
the `HybridAtJunction` tangent continuation (`G1` at the joint, valid only where `t(d) < 0`);
the `AdaptivePath` bounded-miter construction; the continuing-turn counterexample (switching at
the joint is not universal); closed 2D regions, 3D prism closure, general simple polygon
profiles, and direct B-rep closure for Line-Arc at 30/90/150/179 degrees.

### Connections when the analytical seam does not close

These were already proposed in the research and belong to this plan, not to a new task:

| Strategy | When it applies | Status |
| --- | --- | --- |
| AOIM | `d >= d_min`; the offsets meet inside the finite spans | exact, preferred |
| `HybridAtJunction` | `d < d_min` on the inflected side, `t(d) < 0` | validated in 2D/3D; area grows |
| `AdaptivePath` (and event-driven variant) | long sharp wedge unacceptable and a path change is allowed | validated in 2D; `G0` at `A`, bounded miter and deviation |
| `HybridAtLimitTangent`, `HybridAtLimitBisector` | none | diagnostic only; do not use |

The strategy is chosen from geometry before construction: the existence limit, the sign of
`t(d)` and the continuing-turn curvature sign are events, not exceptions. Do not implement these
as try-then-fallback or as a chain of attempts: the kernel forbids that (`AGENTS.md`, "A
fallback is a code smell"). Decide the construction up front and build it once.

## Fixed fields

| Field | Question it answers |
| --- | --- |
| Original goal | Which mathematical problem are we solving? |
| Expected result | What exactly must exist when we are done? |
| Scope boundaries | What are we deliberately not building? |
| Success criteria | By which measurements do we accept the result? |
| Changes needing sign-off | Which direction changes are not ours to make alone? |

### Original goal

Compute the geometry of the miter joint of a planar Sweep analytically, **before any Sweep
body exists**, from the original profile and the mathematical description of the path.

### Expected result

A specification and an algorithm covering:

1. **Path Analysis**: spans, tangents, normals, curvature, joints, travel direction.
2. **Profile Projection**: place the original profile in the path's local frames without
   changing its geometry.
3. **Offset Generation / Sweep Footprint**: the profile's trace along neighbouring spans and
   the offset curves the joint needs.
4. **Joint Classification**: inner or outer corner, turn direction, tangency.
5. **Analytical Miter Construction**: the joint geometry itself, including curved paths.
6. **Miter Clipping**: which parts of the profile each Sweep span keeps or clips.
7. **B-rep Construction**: hand the prepared geometry to the existing SolidCore surface,
   trim, stitch and validity machinery.

The algorithm must preserve the profile's original shape regardless of how many lines,
arcs, rounds or other curves it contains, and must not depend on a particular profile type
or on which surface kind the result happens to produce (cylinder, torus, sphere, other).

The final result is a correct miter on the test matrix below, not merely a valid B-rep.

### Scope boundaries

| Not doing | Reason |
| --- | --- |
| Building Sweep bodies to discover the miter | The miter is computed before the bodies |
| Using `brep_intersect_mut` to recover the seam | That is post-hoc geometry |
| Building our own B-rep kernel, Wire/Face/Shell/Solid builders | SolidCore owns this |
| Extending Union Sweep to serve the analytical miter | Production Union already exists |
| Building a new Sweep or B-rep pipeline from scratch | SolidCore's existing pipeline is the base; only the joint construction changes |
| Rewriting the reference `create_sweep_along_planar_path` | It stays an independent behavioural reference |
| Treating a valid topology as proof of a correct miter | Validity is not geometric correctness |

### Success criteria

- Independent geometric oracle, not `check_brep(Level2)` or `is_manifold_solid` alone:
  manifoldness, closure, volume against an independent reference, gaps at joints,
  coincident boundaries, real self-intersections.
- The profile's original shape is preserved for any profile composition.
- Arc + line and arc + arc are correct where the reference falls back to a chord or produces
  an invalid result.
- Several curved spans, closed paths, multiple contours and holes behave correctly.
- No case passes by being fitted to one test; the matrix is fixed before the result is seen.

### Test matrix

1. Two line spans, 90 degrees, acute and obtuse angles.
2. Arc + line.
3. Arc + arc.
4. Several curved spans on one path.
5. Profile with many rounds and straight spans.
6. Closed paths of arcs, and of arcs plus lines.
7. Closed curved path with several profiles.
8. Multiple contours and holes; different `cap_start` / `cap_end` combinations.

Pilot scope: joints of type line-line and line-arc (inflected and continuing), including several
such joints and closed paths. Arc-arc joints, general profiles, holes and multiple contours are
later milestones.

### Changes that need sign-off

- Changing the original goal, expected result, scope boundaries or success criteria.
- Reusing Union (or any post-hoc body operation) as the analytical miter.
- Dropping or weakening the independent oracle.
- Using a connection strategy as a silent replacement for the analytical miter, or implementing
  it as try-then-fallback instead of a construction chosen from the geometry.
- Changing one of the fixed decisions above, adding a construction to the pilot, or mixing
  constructions inside one path.
- Replacing the existing SolidCore Sweep pipeline instead of changing its joint construction.

## Rules

1. **A technical limitation does not change the mathematical problem.** If SolidCore lacks
   an operation, that is a limit of the tool, not a reason to redefine the goal.
2. **Keep research and engineering workarounds separate.** A useful workaround is recorded
   as its own result; it is never announced as the answer to the original question.
3. **Re-check the direction after every stage.** Ask: are we still solving the original
   problem, or something else? If it is something else, stop and get agreement.

## Current status of the threads

| Thread | Status |
| --- | --- |
| Union Sweep (`create_sweep_along_planar_path_union`) | Complete separate engineering result and production planar Sweep. Not the answer to this research. Its known joint volume loss is `4/3 * r^3` per corner for a circular profile. |
| Hybrid (`create_sweep_along_planar_path_hybrid`) | Production alias of Union. |
| Analytical miter | This research. Line-Arc AOIM validated in `research/aoim-line-arc`; production integration not started. |
| Imprint-and-loop | Alternative direction, not in the main plan. |

## Corrected plan (continues the existing research)

The former blocker is closed: the model is not missing, it is the `research/aoim-line-arc`
branch. Do not re-derive the validated facts above. Union remains production; this work is
research.

### Fixed decisions (2026-11-10)

| Decision | Ruling |
| --- | --- |
| AdaptivePath | Research-only. It changes the path and the result area; production needs a separate API contract with explicit user opt-in. |
| AOIM and Union | Parallel construction. Union stays production. AOIM does not replace it now. No automatic Union fallback and no try-then-fallback. The construction is chosen from the input geometry and the declared strategy before building. |
| Profile generalization | Its own milestone, not part of the pilot. B-splines, several contours, holes, offset topology, self-intersections and split/merge need separate criteria. Line-Arc and the simple polygon remain the validated base. |
| 3D target | Path-frame planar Sweep: the profile carried along the planar path in its local frames. `Sweep3D` is an independent cross-check only; a fixed global-Z extrusion is not the production basis. |
| Joint Classifier | Mandatory Stage 0 result. It selects exactly one construction before geometry and returns an explicit error if that construction cannot be built. It never tries the next strategy. |
| Implementation base | SolidCore's existing planar Sweep pipeline (`create_sweep_along_planar_path`): profile placement, per-segment sweep, trim, stitch, classification and B-rep construction are reused. Only the joint construction changes. No new Sweep or B-rep pipeline. |

### Research still needed

The Line–Arc algebra is settled. The remaining research is:

1. **Joint Classifier** (Stage 0): joint type, AOIM domain validity, `sign(t(d))`, continuation
   direction, path-change permission, profile applicability; one construction per path or an
   explicit error.
2. **Path-frame 3D joint geometry** (Stage 2): express the seam and clipping in the profile's
   local frames along the path, not in a fixed global-Z extrusion.
3. **Seam realization in the existing pipeline** (Stage 0 decision, exercised in Stage 3):
   whether `trim_both` can take the AOIM seam, or the joint uses direct shared topology while
   the rest of the pipeline stays.
4. **Continuing-turn at an invalid AOIM domain**: no validated construction yet; `Unsupported`
   until a separate decision.
5. **Independent validation** of the integrated body (Stage 5).
6. **General profiles**: deferred milestone, not the pilot.

### Stage 0: disposition, classifier and boundary (no code)

- Classify each module of `research/aoim-line-arc/src` and of the existing planar Sweep work:
  math to keep, prototype to port, validation-only, or discard.
- Write the formal Joint Classifier:
  - input: joint type (inflected or continuing), AOIM domain validity (`d` against `d_min` /
    `d_max`), the sign of `t(d)`, continuation direction, path-change permission, profile
    applicability;
  - output: exactly one construction per path, or an explicit `Unsupported` error;
  - the decision is an event of the input geometry, never an exception or a retry.
- Decide the seam interface: whether `trim_both` can take the AOIM seam, or the joint uses direct
  shared topology while the rest of the pipeline stays. Record the minimal interface either way.
- Write the path-frame 3D formulation plan: how the joint seam and clipping are expressed in the
  profile's local frames.
- Fix the research/production boundary: the pilot stays research; Union stays production.
- Plan the transfer of mathematical components only. Do not transfer the standalone
  `Point2`/`Point3`, the mesh builders or the experimental region builders.
- List the SolidCore APIs the pilot needs and the blocked ones as dependencies, not as a reason
  to change the goal.

Preliminary disposition:

| Module | Disposition |
| --- | --- |
| `src/lib.rs` | AOIM algebra and domain facts: keep and port. `Point2` and errors: rewrite against the target types. |
| `src/continuing.rs` | Opposite-curvature model: keep and port; the classifier needs it. |
| `src/adaptive.rs` | Adaptive path: keep in research; not a production construction. |
| `src/sweep.rs` | 2D region builder and raster oracle: validation harness, not a production builder. |
| `src/sweep3d.rs` | Vertical-extrusion mesh: independent cross-check only. |
| `examples/` | Diagnostics and figures: stay in the research branch. |
| `probes/brep_probe` | Direct B-rep probe: keep as the integration reference. |

### Stage 1: port the validated mathematical model

- Preserve the AOIM algebra, domain facts and the independent tests.
- Rewrite the prototype's types and builders against the target types; do not copy them.
- Keep `Analytical`, `HybridAtJunction` and `ContinuingLineArc`.
- One construction is selected from geometry; no trial sequence.

### Stage 2: path-frame joint geometry in 3D

- Target is the profile carried along the planar path in its local frames, not a fixed global-Z
  extrusion.
- Express the AOIM seam and clipping in that geometry; compute them before building any sweep
  body.
- The research's 3D prism and `Sweep3D` are independent cross-checks only.
- Output is the seam/clipping contract the pilot consumes, not a B-rep.

### Stage 3: pilot in the existing planar Sweep

This is not a profile-generalization stage and not a new Sweep method. The base is SolidCore's
existing `create_sweep_along_planar_path`: profile placement, per-segment sweep, `trim_both`,
stitch, classification and B-rep construction stay. Change only the joint cutting geometry and
the minimal interface that passes the prepared seam to it.

- Do not rewrite `create_sweep_along_planar_path` in place and do not modify
  `create_sweep_along_planar_path_union`. The pilot is a separate research entry point that
  reuses the existing pipeline; the reference stays the control and Union stays untouched. Union
  is a separate result, not the pilot base.
- In kernel integration the change lands in the one construction, not as a permanent parallel
  copy. The parallel entry point is a research tool.
- The pilot consumes the Stage 2 seam contract; it does not recompute the seam with post-hoc
  intersections.
- Supported: planar paths of lines and arcs with inflected line-arc joints (and line-line). Use
  `Analytical` where the offset domain is valid and `HybridAtJunction` where it is not.
- Continuing-turn joints are covered only where the AOIM domain is valid. Where it is not, the
  pilot returns `Unsupported`: the research refuted the inflected hybrid there and has not yet
  chosen an alternative construction. Choosing one is a separate decision.
- One construction per path, chosen before building. Do not mix constructions inside one path.
- Unsupported joints return an explicit error. The pilot does not fall back to Union.
- Compare the pilot with Union on the same inputs: seam geometry, volume and the union lens,
  validity, manifoldness, gaps and self-intersections.
- No rewrite of the existing Sweep and no new profile support (B-splines, holes, several
  contours). Those are the deferred milestone.

### Stage 4: SolidCore integration

- Use the direct shared-topology path from `BREP_PLAN.md`: no Boolean and no internal seam face.
- The pilot may call SolidCore APIs from the research code; merging into the kernel is a
  separate gated step under the kernel's rules (no fallback chain, no flags, tests and goldens).
- Do not use `StitchFaces`; its path is blocked earlier by kernel storage work.

### Stage 5: validation

- Independent oracle: boundaries, volume, manifoldness, gaps, self-intersections, and every
  construction the classifier can select.
- Do not use the pilot's seam code as the oracle. For `HybridAtJunction`, the expected area and
  volume change is part of that construction; measure against the selected construction as well
  as the smooth pipe.
- Keep the research's "not yet validated" list: holes, the general-rectangle residual,
  twisted/moving frames, true trimmed NURBS, B-rep sewing.
- Do not claim volume or outward normals while `brep_volume` and `orient_trimmed_surfaces` are
  blocked.
- Separate engineering dependencies, tracked as their own work: the `topology.brep_p6` tail of
  `TestClosureSplitRegion`, `query.brep::CalculateBoundingBox`, and
  `intersect.surface_drop_lib_p3::SmSurface::DropCurve`. They block validation, not the direction.

### Deferred milestone: general profiles

Separate research with its own criteria: B-splines, several contours, holes, offset topology,
self-intersections, split/merge boundaries, geometric accuracy and topological correctness.
Not part of the pilot.

## Per-stage review checklist

- Which original sub-problem does this stage solve?
- Does it compute miter geometry before any Sweep body exists?
- Does it preserve the profile's original shape for an arbitrary profile?
- Does it avoid relying on post-hoc intersections of built bodies?
- Is the success measurement independent of the implementation?
- Does the change reuse the existing SolidCore pipeline, changing only the joint construction?
- Does the path use exactly one construction, chosen before building?
- Does an unsupported or ambiguous case return an explicit error instead of trying another
  strategy?
- If any answer is no, stop and request sign-off before continuing.
