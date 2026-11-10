//! **Status: superseded research prototype (2026-11-10).** This recovers the seam from
//! the surface-surface intersection of already-built sweep bodies. The current charter
//! forbids that post-hoc geometry (`plan/analytical-miter.md`); keep this module for
//! comparison only, not as the AOIM pilot.
//!
//! Experimental planar-path sweep that replaces the kernel's two mitre-construction
//! sites with one analytic miter curve. This is a parallel, independent copy of
//! [`super::planar_path`]'s pipeline: the profile placement, the per-segment sweep,
//! `trim_both`, the union, the stitch and the classification are reused as-is. Only
//! the cutting geometry at a joint is computed differently.
//!
//! The reference implementation never changes; this module exists so the two can be
//! compared case by case. `create_sweep_along_planar_path_analytic_miter` differs from
//! the reference in exactly two places:
//!
//! - the junction between consecutive segments (including the non-linear branch,
//!   which the reference approximates by a straight chord through two B-rep
//!   intersection points);
//! - the closing junction of a closed path.
//!
//! The two cases are deliberately different:
//!
//! - Two straight segments: [`build_analytic_miter_curve`] builds the angle bisector
//!   of the one-sided tangents through the corner, parallel to the path plane normal.
//!   That is exactly the reference's `temp_cutting` line, so straight-segment results
//!   agree with the reference bit for bit, independent of the profile.
//! - Any curved segment: [`build_seam_loop_cutter`] keeps the true surface-surface
//!   intersection that the reference throws away. It gets the kernel's intersection
//!   curves between the two swept bodies, takes the two halves nearest the tangent
//!   bisector (together the seam loop about the tube), and joins them into one closed
//!   B-spline. This is the analytic seam, not a plane through a chord.
//!
//! **Known blocker (measured, not assumed).** `project_and_trim` imprints a curve by
//! parallel projection along one direction and classifies the deleted side against
//! that curve's plane. On a curved joint the seam is a non-planar loop, and a plane
//! cuts the two swept surfaces in two different cross-sections, so the trimmed pieces
//! cannot abut: the union stays open. A single arc sweep is exact (verified), so the
//! limitation is the trim primitive, not the sweep or the seam computation. Realizing
//! the seam needs an imprint-and-loop trim (UV trim curves merged on each face, side
//! chosen by the loop), which is a change to the trim primitive. Until then the
//! curved-joint comparison is recorded by the ignored matrix test.
//!
//! Not upstream (a solidcore experiment): the reference is `super::planar_path`.

use super::{as_bspline, PrimitiveCreationAccess};
use crate::curve::curve::Curve;
use crate::error::{checked_or_ignored, ignored, Result, SolidError};
use crate::geom::extent3d::Extent3d;
use crate::geom::geom_utility::p2::ContinuityType;
use crate::kernel::Kernel;
use crate::math::Vec3;
use crate::nurbs::b_spline_curve::p1::BSplineCurve;
use crate::nurbs::b_spline_curve::p2::BoundingBoxRequest;
use crate::nurbs::crv_adv::p1::CurveEnd;
use crate::primitive::primitive_creation::p1::{
    extend_path_curves_with_kernel, position_profile_curves_impl, PrimitiveCreation,
};
use crate::primitive::primitive_creation::SweepOptions;
use crate::topology::brep::p10::TrimType;
use crate::topology::BrepId;

const OP: &str = "primitive::primitive_creation::p2::planar_path_analytic_miter";

/// `Evaluate(t, 1, from_left, out)` (non-zero tangents) whose ignored status keeps
/// the previous position and tangent.
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

/// `CreateExtendedCurve(dist, end, SM_CT_G1, out)` whose ignored status leaves `out`.
fn extend_into(
    checked: bool,
    curve: &BSplineCurve,
    distance: f64,
    end: CurveEnd,
    out: &mut Option<BSplineCurve>,
) -> Result<()> {
    if let Some(Some(next)) = checked_or_ignored(
        curve.create_extended_curve_by_distance(distance, end, ContinuityType::G1),
        checked,
    )? {
        *out = Some(next);
    }
    Ok(())
}

/// `CreateLineSegment(3, start, end, out)` whose ignored status leaves `out`.
fn line_into(checked: bool, start: Vec3, end: Vec3, out: &mut Option<BSplineCurve>) -> Result<()> {
    if let Some(line) =
        checked_or_ignored(BSplineCurve::create_line_segment(3, start, end), checked)?
    {
        *out = Some(line);
    }
    Ok(())
}

/// Not upstream (WR-400, tangency-labs/solidcore#837): whether the profile curves
/// close, every curve end meeting another curve end (a closed curve's own other end)
/// within `tol`.
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

/// Whether `SweepProfilesAlongPath` asks `CreateCurveSweep` for caps: always upstream.
/// With `Kernel::fixes.planar_path_uncapped_open_profiles` (not upstream, WR-400,
/// tangency-labs/solidcore#837) only for a profile that closes.
fn cap_profile_ends(k: &Kernel, profiles: &[Curve], tol: f64) -> Result<bool> {
    if !k.fixes.planar_path_uncapped_open_profiles {
        return Ok(true);
    }
    profile_closes(k, profiles, tol)
}

/// `SweepProfilesAlongPath`: sweep the positioned profiles along one path segment,
/// always capped (see [`cap_profile_ends`]), then remove the unwanted caps.
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

/// Reference points on either side of a cutting line, then trim the B-reps on both
/// sides (`ProjectAndTrim(..., SM_TT_DELETE_POINT, &dTol, NULL, NULL)`, ignored).
fn trim_both(
    a: &mut impl PrimitiveCreationAccess,
    breps: [BrepId; 2],
    cutting: &BSplineCurve,
    normal: Vec3,
    references: [Vec3; 2],
    tol: f64,
) -> Result<()> {
    let [before, after] = breps;
    let checked = a.kernel(before)?.fixes.planar_path_sweep_statuses;
    let [reference1, reference2] = references;
    checked_or_ignored(
        a.project_and_trim(
            before,
            cutting,
            normal,
            reference2,
            TrimType::DeletePoint,
            tol,
        ),
        checked,
    )?;
    checked_or_ignored(
        a.project_and_trim(
            after,
            cutting,
            normal,
            reference1,
            TrimType::DeletePoint,
            tol,
        ),
        checked,
    )
    .map(|_| ())
}

/// Compute the normal only before the first operation that consumes it.
fn path_plane_normal(k: &Kernel, paths: &[Curve], cached: &mut Option<Vec3>) -> Result<Vec3> {
    if let Some(normal) = *cached {
        return Ok(normal);
    }
    let normal = crate::curve::curve::p3::compute_crvs_normal(k, paths)?;
    *cached = Some(normal);
    Ok(normal)
}

/// The analytic mitre cutter at a tangent discontinuity.
///
/// The cut is a plane parallel to the path plane normal through the corner; its trace
/// on the path plane is what this returns as the cutting curve. The trace direction is
/// the bisector of the incoming and outgoing one-sided tangents, i.e. the classical
/// miter plane of the two tangent lines. For two straight segments this plane is the
/// exact mitre (it is the reference's `temp_cutting` line), independent of the profile
/// shape.
///
/// `length` is the half-length the cutter is extended to on each side, so the plane
/// reaches past the profile.
///
/// Curved segments are the open case: their exact seam is not this plane. Keeping the
/// construction in one function is what lets that be replaced without touching the
/// surrounding pipeline.
fn build_analytic_miter_curve(
    checked: bool,
    corner: Vec3,
    incoming_tangent: Vec3,
    outgoing_tangent: Vec3,
    _path_normal: Vec3,
    length: f64,
    _tol: f64,
) -> Result<BSplineCurve> {
    // The bisector of the two tangents, oriented across the corner. A straight cutter
    // is symmetric, so the sign of the average does not change the line.
    let average = ((-incoming_tangent) + outgoing_tangent).unitize_or_z();
    let cutting_point = corner + (average * length);
    let mut line = None;
    line_into(checked, corner, cutting_point, &mut line)?;
    let line = line.ok_or_else(|| SolidError::null_pointer(OP))?;
    let mut extended = None;
    extend_into(checked, &line, length, CurveEnd::Start, &mut extended)?;
    extended.ok_or_else(|| SolidError::null_pointer(OP))
}

/// Evenly spaced samples of `curve` over its natural interval.
fn sample_curve(k: &Kernel, curve: &Curve, count: u32) -> Result<Vec<Vec3>> {
    let domain = curve.natural_interval();
    let mut points = Vec::new();
    for j in 0..=count {
        let f = f64::from(j) / f64::from(count);
        let t = domain.min() + (domain.max() - domain.min()) * f;
        points.push(curve.evaluate_point(k, t)?);
    }
    Ok(points)
}

/// Distance from `p` to the line through `origin` along the unit `direction`.
fn distance_to_line(p: Vec3, origin: Vec3, direction: Vec3) -> f64 {
    let w = p - origin;
    p.distance(origin + direction * w.dot(direction))
}

/// Build the analytic miter seam as a closed B-spline from the kernel's
/// surface-surface intersection curves.
///
/// The two curves whose samples lie closest to the tangent-bisector line are the two
/// halves of the seam loop about the tube, so their samples are joined into one closed
/// loop and interpolated. `None` when no two curves can be paired.
#[allow(clippy::cast_precision_loss)]
fn build_seam_loop_cutter(
    k: &Kernel,
    intersections: &[Curve],
    corner: Vec3,
    incoming_tangent: Vec3,
    outgoing_tangent: Vec3,
) -> Result<Option<BSplineCurve>> {
    let direction = ((-incoming_tangent) + outgoing_tangent).unitize_or_z();
    let mut scored: Vec<(f64, usize)> = Vec::new();
    for (index, curve) in intersections.iter().enumerate() {
        let points = sample_curve(k, curve, 8)?;
        let mut best = f64::INFINITY;
        for &p in &points {
            let dist = distance_to_line(p, corner, direction);
            if dist < best {
                best = dist;
            }
        }
        if best.is_finite() {
            scored.push((best, index));
        }
    }
    scored.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal));
    let Some(&(_, first_index)) = scored.first() else {
        return Ok(None);
    };
    let Some(&(_, second_index)) = scored.get(1) else {
        return Ok(None);
    };
    let first = sample_curve(
        k,
        intersections
            .get(first_index)
            .ok_or_else(|| SolidError::sm_err(OP))?,
        16,
    )?;
    let second = sample_curve(
        k,
        intersections
            .get(second_index)
            .ok_or_else(|| SolidError::sm_err(OP))?,
        16,
    )?;
    let first_end = *first.last().ok_or_else(|| SolidError::sm_err(OP))?;
    let second_start = *second.first().ok_or_else(|| SolidError::sm_err(OP))?;
    let second_end = *second.last().ok_or_else(|| SolidError::sm_err(OP))?;
    let mut points = first.clone();
    if second_start.distance(first_end) <= second_end.distance(first_end) {
        points.extend(second.iter().skip(1).copied());
    } else {
        points.extend(second.iter().rev().skip(1).copied());
    }
    // a closed interpolation must not repeat its first point at the end
    if points.len() > 1 {
        let start = points.first().copied();
        let end = points.last().copied();
        if let (Some(start), Some(end)) = (start, end) {
            if start.distance(end) < 1e-9 {
                points.pop();
            }
        }
    }
    if points.len() < 4 {
        return Ok(None);
    }
    let curve = BSplineCurve::interpolate_points(
        &points,
        None,
        3,
        None,
        None,
        true,
        crate::nurbs::crv_fit::p1::ParamMethod::ChordLength,
    )?;
    Ok(Some(curve))
}

impl PrimitiveCreation {
    /// [`PrimitiveCreation::create_sweep_along_planar_path`] with an analytic mitre
    /// cutter instead of the reference's straight chord through two intersection
    /// points. See the module documentation for what is shared and what differs.
    ///
    /// # Errors
    ///
    /// Those of [`PrimitiveCreation::create_sweep_along_planar_path`].
    pub fn create_sweep_along_planar_path_analytic_miter(
        a: &mut impl PrimitiveCreationAccess,
        context_of: BrepId,
        profiles: &[Curve],
        paths: &[Curve],
        tol: f64,
        cap_start: bool,
        cap_end: bool,
        move_profile_center_to_path: bool,
    ) -> Result<BrepId> {
        let checked = a.kernel(context_of)?.fixes.planar_path_sweep_statuses;
        let mut owned = Vec::new();
        let result = Self::create_sweep_along_planar_path_analytic_miter_owned(
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
        if checked && result.is_err() {
            for id in owned {
                a.delete_brep(id);
            }
        }
        result
    }

    #[allow(clippy::too_many_arguments)]
    fn create_sweep_along_planar_path_analytic_miter_owned(
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
        let [start_point, start_tangent] = start_pv;
        let [end_point, end_tangent] = end_pv;
        let closed = start_point.distance(end_point) < tol;
        let periodic = closed && start_tangent.is_parallel_to(end_tangent, 2.0);

        // Determine max size of profile
        let k = a.kernel(context_of)?;
        let mut total = Extent3d::new();
        let mut bounds = Extent3d::new();
        for profile in profiles {
            if let Some(boxes) = checked_or_ignored(
                profile.calculate_bounding_box(
                    k,
                    &profile.natural_interval(),
                    BoundingBoxRequest {
                        normal: true,
                        ..BoundingBoxRequest::default()
                    },
                    true,
                ),
                checked,
            )? {
                if let Some(normal) = boxes.normal {
                    bounds = normal;
                }
            }
            total.add_point(bounds.min());
            total.add_point(bounds.max());
        }
        let max_dimension = total.max_dimension();
        let mut extend = max_dimension * 2.0_f64.sqrt();
        extend *= 2.0;

        // Extend all the path curves at start and end to make sure sweeps go past corner
        let extended = if checked {
            super::super::p1::extend_path_curves_impl(
                &bspline_paths,
                closed,
                periodic,
                extend,
                true,
            )?
        } else {
            extend_path_curves_with_kernel(k, &bspline_paths, closed, periodic, extend)?
        };

        // A plane normal is only needed at a mitre. A single straight path has
        // no unique normal (`ComputeCrvsNormal` returns SM_ERR) and needs none.
        // Legacy mode retains the original eager call and undefined fallback.
        let mut path_normal = if checked {
            None
        } else {
            Some(
                ignored(crate::curve::curve::p3::compute_crvs_normal(k, paths))?
                    .unwrap_or(Vec3::UNDEF),
            )
        };

        // Prepare for loop where we position profile curves, sweep, and trim
        let mut new_breps: Vec<BrepId> = Vec::new();
        let mut before_brep: Option<BrepId> = None;
        let mut first_brep: Option<BrepId> = None;
        let mut before = [Vec3::UNDEF; 2];
        let mut after = [Vec3::UNDEF; 2];
        let mut first_pv = [Vec3::UNDEF; 2];

        for (ii, (path, extended_path)) in paths.iter().zip(&extended).enumerate() {
            // find start point and derivative of new path segment
            let range = path.natural_interval();
            evaluate_keep(a.kernel(context_of)?, path, range.min(), false, &mut after)?;
            after[1] = after[1].unitize_or_z();

            // position profile curves at the start of the path segment
            let mut new_profiles = Vec::new();
            position_profile_curves_impl(
                a.kernel(context_of)?,
                profiles,
                tol,
                &Curve::from(extended_path.clone()),
                &mut new_profiles,
                checked
                    || a.kernel(context_of)?
                        .fixes
                        .planar_path_profile_placement_statuses,
                paths,
            )?;

            // sweep the profiles along this path segment creating new brep
            let after_brep = a.new_brep(context_of)?;
            owned.push(after_brep);
            let cap_this_start = ii == 0 && !closed && cap_start;
            let cap_this_end = ii == paths.len() - 1 && !closed && cap_end;
            sweep_profiles_along_path(
                a,
                &new_profiles,
                extended_path,
                tol,
                cap_this_start,
                cap_this_end,
                move_profile_center_to_path,
                after_brep,
            )?;
            new_breps.push(after_brep);

            // The first segment only records its derivatives; later ones cut both
            // breps at the mitre.
            if ii == 0 {
                first_brep = Some(after_brep);
                first_pv = after;
            } else if !before[1].is_parallel_to(after[1], 2.0) {
                let path_normal =
                    path_plane_normal(a.kernel(context_of)?, paths, &mut path_normal)?;
                let before_id = before_brep.ok_or_else(|| SolidError::null_pointer(OP))?;
                let previous = ii
                    .checked_sub(1)
                    .and_then(|i| paths.get(i))
                    .ok_or_else(|| SolidError::sm_err(OP))?;
                let both_linear = path.is_linear(a.kernel(context_of)?, tol, false)?.linear
                    && previous
                        .is_linear(a.kernel(context_of)?, tol, false)?
                        .linear;
                let cutting = if both_linear {
                    build_analytic_miter_curve(
                        checked,
                        after[0],
                        before[1],
                        after[1],
                        path_normal,
                        extend,
                        tol,
                    )?
                } else {
                    let search =
                        a.with_breps_mut(before_id, after_brep, |k, before, after, _, _, _, _| {
                            crate::topology::brep::p1::brep_intersect_mut(
                                k, before, after, tol, 2.0, false,
                            )
                        });
                    let intersections = if checked {
                        Some(search?)
                    } else {
                        ignored(search)?
                    }
                    .map(|(curves, _)| curves)
                    .unwrap_or_default();
                    match build_seam_loop_cutter(
                        a.kernel(context_of)?,
                        &intersections,
                        after[0],
                        before[1],
                        after[1],
                    )? {
                        Some(curve) => curve,
                        None => build_analytic_miter_curve(
                            checked,
                            after[0],
                            before[1],
                            after[1],
                            path_normal,
                            extend,
                            tol,
                        )?,
                    }
                };
                // Reference points on either side of the cutting line.
                let reference1 = before[0] + (-before[1] * 0.1);
                let reference2 = after[0] + (after[1] * 0.1);
                trim_both(
                    a,
                    [before_id, after_brep],
                    &cutting,
                    path_normal,
                    [reference1, reference2],
                    tol,
                )?;
            }

            // Get ready for next path segment
            before_brep = Some(after_brep);

            // find end derivative of path before
            evaluate_keep(a.kernel(context_of)?, path, range.max(), true, &mut before)?;
            before[1] = before[1].unitize_or_z();

            // a closed path trims the end of the last and the beginning of the first
            if closed && ii == paths.len() - 1 {
                let after_id = first_brep.ok_or_else(|| SolidError::null_pointer(OP))?;
                after = first_pv;

                // tangent continuous paths need no trim
                if !before[1].is_parallel_to(after[1], 2.0) {
                    let path_normal =
                        path_plane_normal(a.kernel(context_of)?, paths, &mut path_normal)?;
                    before[1] = before[1].unitize_or_z();
                    after[1] = after[1].unitize_or_z();
                    let cutting = build_analytic_miter_curve(
                        checked,
                        after[0],
                        before[1],
                        after[1],
                        path_normal,
                        50.0,
                        tol,
                    )?;

                    let before_id = before_brep.ok_or_else(|| SolidError::null_pointer(OP))?;
                    let (reference1, reference2) = if before_id == after_id {
                        (
                            before[0] + (before[1] * 10.0),
                            after[0] + (-after[1] * 10.0),
                        )
                    } else {
                        (
                            before[0] + (-before[1] * 30.0),
                            after[0] + (after[1] * 30.0),
                        )
                    };
                    if before_id == after_id {
                        checked_or_ignored(a.stitch_faces(before_id, 0.01), checked)?;
                        checked_or_ignored(
                            a.project_and_trim(
                                before_id,
                                &cutting,
                                path_normal,
                                reference2,
                                TrimType::Split,
                                tol,
                            ),
                            checked,
                        )?;
                        let mut remove = Vec::new();
                        let brep = a.brep(before_id)?;
                        for edge in crate::topology::brep::p1::edges(brep, None)? {
                            if crate::topology::edge::is_lamina(brep, edge)? {
                                remove.extend(crate::topology::edge::faces(brep, edge, None)?);
                            }
                        }
                        checked_or_ignored(a.remove_faces(before_id, &remove), checked)?;
                        checked_or_ignored(
                            a.remove_topological_edges_and_vertices(before_id),
                            checked,
                        )?;
                    } else {
                        trim_both(
                            a,
                            [before_id, after_id],
                            &cutting,
                            path_normal,
                            [reference1, reference2],
                            tol,
                        )?;
                    }
                }
            }
        }

        // merge all results into one brep: the last first, then a union per piece
        let count = new_breps.len();
        let mut ordered = Vec::with_capacity(count);
        let mut tree: Vec<i64> = vec![0];
        ordered.extend(new_breps.pop());
        for ii in 1..count {
            ordered.extend(new_breps.pop());
            tree.push(i64::try_from(ii).map_err(|_| SolidError::sm_err(OP))?);
            // Convert to different form of Boolean numbering (0 - 1)
            tree.push(-1);
        }
        a.boolean_trees(&mut ordered, &tree)?;
        if ordered.len() != 1 {
            // Something went wrong
            return Err(SolidError::sm_err(OP));
        }
        let sweep = *ordered.first().ok_or_else(|| SolidError::sm_err(OP))?;

        if let Err(e) = a.stitch_faces(sweep, tol) {
            if !checked {
                a.delete_brep(sweep);
            }
            return Err(e);
        }
        checked_or_ignored(a.remove_topological_edges_and_vertices(sweep), checked)?;
        // A second attempt can remove additional edges
        checked_or_ignored(a.remove_topological_edges_and_vertices(sweep), checked)?;

        // Solid intent with a manifold boundary gets the nested-solid parity
        // classification.
        let solid_intent = closed || (cap_start && cap_end);
        let manifold = a.is_manifold_solid(sweep)?;
        if solid_intent && manifold {
            if let Err(e) = a.with_primitive_brep(sweep, |k, b| {
                crate::topology::brep::p1::set_region_is_void_flags_for_nested_solids(k, b)
            }) {
                if !checked {
                    a.delete_brep(sweep);
                }
                return Err(e);
            }
        }

        // Reject a manifold boundary whose regions are all void.
        let brep = a.brep(sweep)?;
        let mut has_material = false;
        for region in crate::topology::brep::p1::regions(brep, None)? {
            if brep.region(region)?.is_void() == 0 {
                has_material = true;
                break;
            }
        }
        if manifold && !has_material {
            if !checked {
                a.delete_brep(sweep);
            }
            return Err(SolidError::sm_err(OP));
        }
        Ok(sweep)
    }
}
