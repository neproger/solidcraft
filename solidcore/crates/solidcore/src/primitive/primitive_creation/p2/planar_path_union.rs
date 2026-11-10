//! Production planar-path sweep result: built by unioning nominal closed segment sweeps.
//!
//! This is a separate engineering result, not the analytical miter research. It sits
//! alongside [`super::planar_path`] (the reference) and builds each path segment as its
//! own capped closed body, letting the kernel's boolean union form the joint. No trimming
//! at a joint, no seam.
//!
//! For each segment the profile is positioned at the segment's own start and swept
//! along the segment's own (nominal, not extended) path, capped at both ends. Every
//! segment is therefore a closed solid, and `boolean_trees` unions them. Adjacent
//! bodies overlap in a small lens about the joint (their perpendicular end caps cross),
//! which the union counts once. That lens is the measured volume deficit against the
//! smooth-pipe volume `area * path length`; it is a geometric property of the
//! construction, not a tolerance or a build error.
//!
//! **Volume loss (established property of this construction).** Against the smooth-pipe
//! volume `area * path length`, one joint of a circular profile of radius `r` loses
//!
//! ```text
//! О”V = (4/3) r^3
//! ```
//!
//! for a kinked joint, independent of path length and curvature, and the total loss is
//! that per kinked joint (`N` joints: `О”V = (4/3) r^3 N`). It is exactly the volume of
//! the overlap lens where the two perpendicular end caps cross, counted once by the
//! union; for two perpendicular capped cylinders that overlap volume has the closed form
//! `в€«в€«в€«_{xВІ+yВІв‰¤rВІ, yВІ+zВІв‰¤rВІ, x,zв‰Ґ0} dx dy dz = 4/3 rВі`, so the deficit is the exact
//! geometry of the union of the two nominal solids, not a boolean or construction error.
//! Tangent-continuous joints abut on a coincident cap and lose nothing. For a
//! non-circular section the loss is the same kind of lens but its size follows the
//! section shape, not `r` (measured for a square of side `s`: about `(1/4) s^3` per
//! joint). It is a geometric property of the construction, never a tolerance or a
//! segment-build error; the matrix test asserts it for circular profiles.
//!
//! Deviation from the reference: no path extension and no mitre trim. Not upstream
//! (a solidcore experiment).

use super::{as_bspline, PrimitiveCreationAccess};
use crate::curve::curve::Curve;
use crate::error::{checked_or_ignored, Result, SolidError};
use crate::kernel::Kernel;
use crate::math::Vec3;
use crate::nurbs::b_spline_curve::p1::BSplineCurve;
use crate::primitive::primitive_creation::p1::{position_profile_curves_impl, PrimitiveCreation};
use crate::primitive::primitive_creation::SweepOptions;
use crate::topology::BrepId;

const OP: &str = "primitive::primitive_creation::p2::planar_path_union";

/// `Evaluate(t, 1, from_left, out)` (non-zero tangents) whose ignored status keeps the
/// previous position and tangent.
fn evaluate_keep(
    k: &Kernel,
    curve: &Curve,
    t: f64,
    from_left: bool,
    out: &mut [Vec3; 2],
) -> Result<()> {
    if let Some(values) = checked_or_ignored(
        curve.evaluate(k, t, 1, from_left, true),
        k.fixes.planar_path_sweep_statuses,
    )? {
        for (target, value) in out.iter_mut().zip(values) {
            *target = value;
        }
    }
    Ok(())
}

/// Not upstream (WR-400, tangency-labs/solidcore#837): whether the profile curves close.
fn profile_closes(k: &Kernel, profiles: &[Curve], tol: f64) -> Result<bool> {
    let mut ends = Vec::new();
    for curve in profiles {
        let domain = curve.natural_interval();
        ends.push(curve.evaluate_point(k, domain.min())?);
        ends.push(curve.evaluate_point(k, domain.max())?);
    }
    Ok(ends.iter().enumerate().all(|(i, end)| {
        ends.iter()
            .enumerate()
            .any(|(j, other)| j != i && end.distance(*other) < tol)
    }))
}

/// Whether `CreateCurveSweep` is asked for caps: always, unless the opt-in open-profile
/// fix asks only for a profile that closes (see `planar_path`'s `cap_profile_ends`).
fn cap_profile_ends(k: &Kernel, profiles: &[Curve], tol: f64) -> Result<bool> {
    if !k.fixes.planar_path_uncapped_open_profiles {
        return Ok(true);
    }
    profile_closes(k, profiles, tol)
}

/// Sweep the positioned profiles along one path segment, capped per `cap_start` /
/// `cap_end`, then remove the caps that were not wanted.
fn sweep_profiles_along_path(
    a: &mut impl PrimitiveCreationAccess,
    profiles: &[Curve],
    path: &BSplineCurve,
    tol: f64,
    cap_start: bool,
    cap_end: bool,
    move_profile_center_to_path: bool,
    brep: BrepId,
) -> Result<()> {
    let checked = a.kernel(brep)?.fixes.planar_path_sweep_statuses;
    let cap_ends = cap_profile_ends(a.kernel(brep)?, profiles, tol)?;

    // CreateCurveSweep somehow mangles the path curve
    let path = path.copy();

    let mut start_faces = Vec::new();
    let mut side_faces = Vec::new();
    let mut end_faces = Vec::new();
    let options = SweepOptions {
        cap_ends,
        translational_sweep: false,
        move_profile_center_to_path,
        orient_profile_perpendicular_to_path: false,
        move_to_profile: false,
        profile_at_path_start: false,
    };
    let region = a.infinite_region(brep)?;
    let mut primitive = PrimitiveCreation::new(Some(region), 0, false);
    let profiles: Vec<Option<Curve>> = profiles.iter().map(|c| Some(c.clone())).collect();
    checked_or_ignored(
        primitive.create_curve_sweep(
            a,
            brep,
            &profiles,
            &path,
            None,
            None,
            tol,
            options,
            false,
            &mut start_faces,
            &mut side_faces,
            &mut end_faces,
        ),
        checked,
    )?;

    if !cap_start {
        a.set_editing(brep, true)?;
        checked_or_ignored(a.remove_faces(brep, &start_faces), checked)?;
    }
    if !cap_end {
        a.set_editing(brep, true)?;
        checked_or_ignored(a.remove_faces(brep, &end_faces), checked)?;
    }
    Ok(())
}

impl PrimitiveCreation {
    /// Sweep the planar `profiles` along the ordered, connected planar `paths` by
    /// unioning one nominal, closed sweep per segment. See the module documentation.
    ///
    /// The outer ends follow `cap_start` / `cap_end`; every internal joint is capped so
    /// each segment is a closed body before the union. A closed path caps every joint.
    ///
    /// # Errors
    ///
    /// Those of [`PrimitiveCreation::create_sweep_along_planar_path`] for the input
    /// checks, the segment sweeps and the final union, stitch and classification.
    pub fn create_sweep_along_planar_path_union(
        a: &mut impl PrimitiveCreationAccess,
        context_of: BrepId,
        profiles: &[Curve],
        paths: &[Curve],
        tol: f64,
        cap_start: bool,
        cap_end: bool,
        move_profile_center_to_path: bool,
    ) -> Result<BrepId> {
        let mut owned = Vec::new();
        let result = Self::create_sweep_along_planar_path_union_owned(
            a,
            context_of,
            profiles,
            paths,
            tol,
            cap_start,
            cap_end,
            move_profile_center_to_path,
            &mut owned,
        );
        // Production cleanup: release every temporary on a refusal, independent of the
        // status flag (`delete_brep` is idempotent, so an already released id is inert).
        if result.is_err() {
            for id in owned {
                a.delete_brep(id);
            }
        }
        result
    }

    #[allow(clippy::too_many_arguments)]
    fn create_sweep_along_planar_path_union_owned(
        a: &mut impl PrimitiveCreationAccess,
        context_of: BrepId,
        profiles: &[Curve],
        paths: &[Curve],
        tol: f64,
        cap_start: bool,
        cap_end: bool,
        move_profile_center_to_path: bool,
        owned: &mut Vec<BrepId>,
    ) -> Result<BrepId> {
        let checked = a.kernel(context_of)?.fixes.planar_path_sweep_statuses;
        let first_path = paths
            .first()
            .ok_or_else(|| SolidError::invalid_input(OP, "no path curves"))?;
        let last_path = paths
            .last()
            .ok_or_else(|| SolidError::invalid_input(OP, "no path curves"))?;
        let mut bspline_paths = Vec::with_capacity(paths.len());
        for path in paths {
            bspline_paths
                .push(as_bspline(path).ok_or_else(|| {
                    SolidError::invalid_input(OP, "path curve is not a B-spline")
                })?);
        }

        // Ordered, connected path contract: a disconnected path would build several
        // disjoint segment bodies whose union is a non-manifold shell soup. Refuse it
        // before allocating anything.
        for pair in paths.windows(2) {
            let (Some(first), Some(second)) = (pair.first(), pair.get(1)) else {
                continue;
            };
            let end =
                first.evaluate_point(a.kernel(context_of)?, first.natural_interval().max())?;
            let start =
                second.evaluate_point(a.kernel(context_of)?, second.natural_interval().min())?;
            if end.distance(start) > tol {
                return Err(SolidError::invalid_input(
                    OP,
                    "path curves are not connected",
                ));
            }
        }

        // Determine whether path curves form a closed and periodic loop.
        let mut start_pv = [Vec3::UNDEF; 2];
        let mut end_pv = [Vec3::UNDEF; 2];
        evaluate_keep(
            a.kernel(context_of)?,
            first_path,
            first_path.natural_interval().min(),
            true,
            &mut start_pv,
        )?;
        evaluate_keep(
            a.kernel(context_of)?,
            last_path,
            last_path.natural_interval().max(),
            true,
            &mut end_pv,
        )?;
        let [start_point, _start_tangent] = start_pv;
        let [end_point, _end_tangent] = end_pv;
        let closed = start_point.distance(end_point) < tol;

        let last = paths.len().saturating_sub(1);
        let mut new_breps: Vec<BrepId> = Vec::with_capacity(paths.len());
        for (ii, path) in paths.iter().enumerate() {
            // Position the profile at this segment's own start (the nominal path, not
            // an extension): each segment is swept in its own frame.
            let mut new_profiles = Vec::new();
            position_profile_curves_impl(
                a.kernel(context_of)?,
                profiles,
                tol,
                path,
                &mut new_profiles,
                checked
                    || a.kernel(context_of)?
                        .fixes
                        .planar_path_profile_placement_statuses,
                paths,
            )?;

            // Every internal joint is capped; only an open path's outer ends follow the
            // caller. A closed path has no outer ends, so every end is a capped joint.
            let cap_this_start = ii > 0 || (ii == 0 && closed) || cap_start;
            let cap_this_end = ii < last || (ii == last && closed) || cap_end;

            let seg = a.new_brep(context_of)?;
            owned.push(seg);
            sweep_profiles_along_path(
                a,
                &new_profiles,
                bspline_paths
                    .get(ii)
                    .ok_or_else(|| SolidError::sm_err(OP))?,
                tol,
                cap_this_start,
                cap_this_end,
                move_profile_center_to_path,
                seg,
            )?;
            new_breps.push(seg);
        }

        // Union all segment bodies: the last first, then a union per piece (the
        // reference's final merge).
        let count = new_breps.len();
        let mut ordered = Vec::with_capacity(count);
        let mut tree: Vec<i64> = vec![0];
        ordered.extend(new_breps.pop());
        for ii in 1..count {
            ordered.extend(new_breps.pop());
            tree.push(i64::try_from(ii).map_err(|_| SolidError::sm_err(OP))?);
            tree.push(-1);
        }
        a.boolean_trees(&mut ordered, &tree)?;
        if ordered.len() != 1 {
            return Err(SolidError::sm_err(OP));
        }
        let sweep = *ordered.first().ok_or_else(|| SolidError::sm_err(OP))?;

        if let Err(e) = a.stitch_faces(sweep, tol) {
            a.delete_brep(sweep);
            return Err(e);
        }
        checked_or_ignored(a.remove_topological_edges_and_vertices(sweep), checked)?;
        checked_or_ignored(a.remove_topological_edges_and_vertices(sweep), checked)?;

        let solid_intent = closed || (cap_start && cap_end);
        let manifold = a.is_manifold_solid(sweep)?;
        // A solid was expected but the union/stitch left a non-manifold body (a
        // near-tangent joint below the reference's mitre threshold, or a thin sliver the
        // boolean could not resolve). Never hand that back: release it and refuse.
        if solid_intent && !manifold {
            a.delete_brep(sweep);
            return Err(SolidError::sm_err(OP));
        }
        if solid_intent && manifold {
            if let Err(e) = a.with_primitive_brep(sweep, |k, b| {
                crate::topology::brep::p1::set_region_is_void_flags_for_nested_solids(k, b)
            }) {
                a.delete_brep(sweep);
                return Err(e);
            }
        }

        let brep = a.brep(sweep)?;
        let mut has_material = false;
        for region in crate::topology::brep::p1::regions(brep, None)? {
            if brep.region(region)?.is_void() == 0 {
                has_material = true;
                break;
            }
        }
        if manifold && !has_material {
            a.delete_brep(sweep);
            return Err(SolidError::sm_err(OP));
        }
        Ok(sweep)
    }
}
