# Analytical miter: development journal

Status: Stage 0 (specification, no code)
Started: 2026-11-10
Charter: `plan/analytical-miter.md`
Research source: `https://github.com/neproger/solidcore/tree/research/aoim-line-arc` at `9920e32`

This is the running log of the AOIM pilot. The goal, scope and fixed decisions live in the
charter; this file records what was done, what was learned and what failed. A problem is brought
to the owner and solved with them; it is not a reason to change the goal.

## Rules

- Do not change the goal or the plan because an implementation step is hard.
- Record the problem, the evidence and the options before choosing.
- A change to the charter needs explicit sign-off.
- One construction per path, chosen from geometry before building; no try-then-fallback.
- The existing SolidCore pipeline is the base; only the joint construction changes.

## Stage 0 checklist

- [x] disposition matrix
- [x] formal Joint Classifier
- [x] seam interface decision: `trim_both` for planar line-line, direct shared topology for curved joints
- [ ] seam contract: prove the listed topology APIs build one shared boundary
- [x] path-frame 3D formulation plan
- [x] research/production boundary
- [x] SolidCore API and blocker list

Stage 0 specification: `plan/analytical-miter-stage0.md`.

## 2026-11-10 — Stage 0 starts

Read before starting:

- `plan/analytical-miter.md` (charter and fixed decisions)
- `research/aoim-line-arc` (README, BREP_PLAN, BREP_REPORT_RU)
- `p2/planar_path.rs` (the reference pipeline)
- `p2/planar_path_analytic_miter.rs` on `research/planar-path-superseded` (the old prototype)

### First finding: `trim_both` cannot take the AOIM seam as it is

- `trim_both` (`planar_path.rs:193`) passes a `BSplineCurve` cutting curve, a plane normal and
  reference points to `project_and_trim`.
- `project_and_trim` imprints by parallel projection along one direction and classifies the
  deleted side against the cutting curve's plane. The superseded prototype measured that a
  non-planar seam loop cannot be trimmed this way: the two swept pieces do not abut and the union
  stays open.
- The old prototype then needed imprint-and-loop, which the charter keeps outside the main plan.

### Options

1. Keep `trim_both` for planar joints only (line-line), and build a curved joint with direct
   shared topology. This is the research's `BREP_PLAN.md` direction.
2. Extend `project_and_trim` to accept a general seam loop (imprint-and-loop). This is a kernel
   change and out of scope for the pilot.
3. Use a planar approximation of the seam. Rejected: it reintroduces the chord this task exists
   to remove.

### Accepted decision (owner, 2026-11-10)

Option 1 is accepted for the pilot:

1. `trim_both` stays for Line–Line, where the seam is planar.
2. For Line–Arc and other supported curved joints, the pilot builds the shared boundary
   directly.
3. The shared boundary must be one consistent topological entity, not two geometrically
   coincident but independent edges.
4. The construction uses SolidCore's APIs: `make_shell_vertex`, `make_wire_edge`,
   `make_face_topology`.
5. `project_and_trim` is not called for the AOIM seam.
6. If the topology cannot be built, the pilot returns an explicit error. No fallback to Union.

This is a seam-realization mechanism, not a change of the analytical construction: the
classifier still chooses one construction per path, and the mechanism follows each joint's
geometry. Direct shared topology uses SolidCore's topology/B-rep APIs, not a new kernel. Union
stays a separate production result and is not a fallback.

**Caveat.** The direction is accepted, but it is not yet proven that these APIs suffice for the
shared boundary. The seam contract must demonstrate it; until then this is a design hypothesis,
not a validated result.

### Next

- Prove the shared-topology API sufficiency in the seam contract.
- Define the exact seam data passed from Stage 2 to Stage 3.
- Keep continuing-turn outside AOIM `Unsupported` until a separate decision.

## 2026-11-10 — Stage 0 specification

Wrote `plan/analytical-miter-stage0.md`:

- disposition matrix: keep the algebra and tests, rewrite the standalone types, keep
  `AdaptivePath` in research, use the region/mesh code only as validation, and do not use the
  superseded post-hoc prototype;
- formal Joint Classifier: one construction family per path, `Unsupported` with a reason
  otherwise, `AdaptivePath` never selected;
- seam realization: `trim_both` for line-line, shared topology for curved joints, no
  `project_and_trim` on the AOIM seam;
- path-frame 3D formulation plan, research/production boundary, API and blocker list.

The seam contract remains open: the direction is accepted, but the listed topology APIs are not
yet proven sufficient. That is the next Stage 0 task, not a reason to change the goal.

## Problems and blockers

None yet.
