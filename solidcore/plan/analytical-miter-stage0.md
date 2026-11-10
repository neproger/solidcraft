# Analytical miter: Stage 0 specification

Status: Stage 0 result (specification, no code)
Charter: `plan/analytical-miter.md`
Journal: `plan/analytical-miter-journal.md`

This file is the Stage 0 output: the disposition matrix, the formal Joint Classifier, the seam
realization rule and the open seam-contract questions. It does not change the goal or the fixed
decisions in the charter.

## 1. Disposition matrix

The Line–Arc algebra is settled in `research/aoim-line-arc`. Only the mathematics and its
independent tests transfer; the standalone types, mesh builders and experimental region builders
stay in the research branch.

| Source | What it is | Disposition |
| --- | --- | --- |
| `research/aoim-line-arc/src/lib.rs` | AOIM algebra, `LineArc`, `Point2`, `SeamLocation`, `GeometryError` | Keep the algebra, domain facts and tests. Port them to the target types; rewrite `Point2` and errors against SolidCore's types. |
| `research/aoim-line-arc/src/continuing.rs` | `ContinuingLineArc`, the opposite-curvature model | Keep and port; the classifier needs its domain and tangent sign. |
| `research/aoim-line-arc/src/adaptive.rs` | `AdaptivePath` and the event-driven variant | Keep in research only. Not ported and never selected by the pilot classifier. |
| `research/aoim-line-arc/src/sweep.rs` | 2D region builder, `Strategy`, raster area/gap/self-intersection checks | Validation harness. Not a production builder and not ported as a constructor. |
| `research/aoim-line-arc/src/sweep3d.rs` | Vertical-extrusion prism, seam-face sharing, spatial audit | Independent cross-check only. Not the target 3D geometry. |
| `research/aoim-line-arc/examples/` | Figures, metrics, adaptive and 3D diagnostics | Stay in the research branch. |
| `research/aoim-line-arc/probes/brep_probe` | Public-API direct B-rep probe | Keep as the reference for the shared-topology seam contract. |
| `research/aoim-line-arc/BREP_PLAN.md` | Direct shared-topology plan | Reference for Stage 4. It describes a vertical extrusion; the path-frame case is new work. |
| `research/aoim-line-arc/README.md`, `BREP_REPORT_RU.md` | Established facts and the blocked kernel paths | Reference; do not re-derive. |
| `p2/planar_path.rs` | The reference planar Sweep pipeline | Base. Reuse; do not rewrite in place. |
| `p2/planar_path.rs::trim_both` | Planar cutting curve + `project_and_trim` | Use for Line–Line joints only. Not for the AOIM seam. |
| `topology::brep::{make_shell_vertex, make_wire_edge, make_face_topology}` | SolidCore's shared-topology APIs | Use for curved AOIM joints; sufficiency is not yet proven. |
| `p2/planar_path_analytic_miter.rs` (superseded) | Post-hoc seam from `brep_intersect_mut` + `project_and_trim` | Do not use in the pilot. Historical prototype only. |
| `p2/planar_path_union.rs`, `p2/planar_path_hybrid.rs` | The Union result and its alias | Separate production result. Not the pilot base and not a fallback. |

## 2. Joint Classifier

### Purpose

Choose exactly one analytical construction for the whole path before any geometry is built. The
decision is a function of the input geometry; it is never a retry or an exception.

### Inputs

- `spans`: the path spans, each `Line` or `Arc`, with their joint angles and radii.
- `joints`: for each joint, its type and local geometry.
- `profile_coverage`: the set of signed offsets the profile presents at the joint. For the pilot
  this is a single interval `[d0, d1]` from a simple contour; other profile classes are
  `Unsupported` (the profile-generalization milestone is deferred).
- `path_change_allowed`: false in the pilot. `AdaptivePath` is research-only.
- `policy`: the fixed pilot policy. It never enables a fallback.

### Joint types

- `LineLine`
- `InflectedLineArc`
- `ContinuingLineArc`
- `ArcArc` (not in the pilot; later milestone)
- `Other` (not supported)

### Domains (from the research)

For an inflected joint with radius `R` and turn angle `alpha`:

- offset circle radius `R + d` must be positive;
- AOIM exists for `d >= d_min = -R cos^2(alpha/2)`;
- the missing side `d < d_min` is covered by `HybridAtJunction` only where
  `t(d) = d tan(alpha/2) < 0`.

For a continuing joint:

- offset circle radius `R - d` must be positive;
- AOIM exists for `d <= d_max = R cos^2(alpha/2)`;
- the missing side `d > d_max` has no validated construction; it is `Unsupported` in the pilot.

### Algorithm

1. If any joint is `ArcArc` or `Other`, return `Unsupported { joint, reason }`.
2. If `profile_coverage` is not a single simple interval, return
   `Unsupported { reason: "profile class" }` (pilot scope).
3. If every joint is `LineLine`, return `Reference`.
4. Otherwise the path has at least one curved joint:
   - For every `ContinuingLineArc` joint, require the whole coverage inside the AOIM domain.
     If it is not, return `Unsupported { joint, reason: "continuing-turn outside AOIM" }`.
   - For every `InflectedLineArc` joint:
     - if the whole coverage is inside the AOIM domain, the joint needs only `Analytical`;
     - else if every missing offset is in the hybrid domain (`d < d_min` and `t(d) < 0`), the
       joint needs `HybridAtJunction`;
     - else return `Unsupported { joint, reason: "inflected coverage" }`.
   - The path construction is:
     - `Analytical` if every curved joint needs only `Analytical`;
     - `HybridAtJunction` if at least one inflected joint needs the hybrid branch, and every
       missing offset is in the hybrid domain;
     - never a mix of `Analytical` and `HybridAtJunction` at path level.
5. `AdaptivePath` is never selected. A request that requires changing the path returns
   `Unsupported { reason: "AdaptivePath is research-only" }`.
6. Return the decision. There is no second attempt.

### Output

```text
enum PathConstruction {
    Reference,                       // every joint is line-line
    Analytical,                      // curved joints, all coverage inside AOIM
    HybridAtJunction,                // curved joints, inflected missing side in the hybrid domain
    Unsupported { joint, reason },
}

enum SeamRealization {
    PlanarTrim,                      // line-line: trim_both / project_and_trim
    SharedTopology,                  // curved: SolidCore topology/B-rep APIs
}
```

The classifier returns `PathConstruction`. `SeamRealization` is derived from each joint's type,
not chosen by a fallback: `LineLine` uses `PlanarTrim`, every curved joint uses
`SharedTopology`. This is a realization mechanism, not a second construction.

### Errors

An unsupported or ambiguous input returns an explicit `Unsupported` with the joint and the
reason. The caller does not try another construction and does not fall back to Union.

## 3. Seam contract: open questions

The direction is accepted, but the API sufficiency is not proven. The seam contract must show:

1. `make_shell_vertex`, `make_wire_edge` and `make_face_topology` can build **one** shared
   boundary for the AOIM seam, not two coincident independent edges.
2. The boundary represents the inflected AOIM seam and its `HybridAtJunction` tangent segment.
3. The resulting faces close with the adjacent sweep bodies, every edge is used consistently, and
   the shell is manifold.
4. The construction works in the path frame, not only in the research's vertical extrusion.
5. If an API is missing, name it and stop. Do not substitute a Boolean or `project_and_trim`.

Evidence: extend `probes/brep_probe` to the path-frame case, with `check_brep(Level2)`,
manifoldness, edge-use and shell checks. Volume and outward normals remain blocked by
`brep_volume` and `orient_trimmed_surfaces`; do not claim them.

## 4. Path-frame 3D formulation plan

- Place the profile with the existing `position_profile_curves_impl`, in the path's local
  frames.
- Compute the AOIM seam in the joint's local 2D frame from the research algebra, then map it into
  the profile's local frame at the joint.
- Parameterize the seam over the profile boundary. The pilot supports a simple contour with one
  offset interval; holes and several contours are the deferred milestone.
- Produce the seam as the shared boundary of the two adjacent sweep ends.
- Cross-check the same 2D inputs against the research's vertical extrusion and its raster/area
  oracle. `Sweep3D` is not the target geometry.

## 5. Research/production boundary

- The pilot is research. It starts on a dedicated `research/aoim-pilot` branch.
- Union stays the production planar Sweep and is not a fallback.
- The pilot may call SolidCore APIs from research code; merging into the kernel is a separate
  gated step under the kernel rules.
- No permanent parallel copy in the kernel: the change lands in the one construction.

## 6. APIs and blockers

Used by the pilot:

- `position_profile_curves_impl`, `sweep_profiles_along_path`;
- `trim_both` / `project_and_trim` for Line–Line;
- `make_shell_vertex`, `make_wire_edge`, `make_face_topology` for curved joints;
- `check_brep`, `is_manifold_solid`.

Blocked validation, tracked separately, not a reason to change the direction:

- `brep_volume` (`intersect.surface_drop_lib_p3::SmSurface::DropCurve`);
- `orient_trimmed_surfaces` (`query.brep::CalculateBoundingBox`);
- the `topology.brep_p6` tail of `TestClosureSplitRegion`.

## 7. Open questions

- Prove the shared-topology API sufficiency (section 3).
- Define the exact seam data passed from Stage 2 to Stage 3 (vertices, curves, edge uses).
- Continuing-turn outside AOIM stays `Unsupported` until a separate construction is chosen.
- General profiles (B-splines, holes, several contours) are a later milestone.
