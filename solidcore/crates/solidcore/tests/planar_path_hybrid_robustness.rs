//! Robustness audit for the hybrid planar sweep. Every construction must satisfy the
//! invariant: a returned body is a valid manifold solid; an input the kernel cannot
//! build is a clean `Err`, never a partial or incorrect body. The cases cover near
//! tangent and degenerate segments, very short segments, profiles with holes and several
//! contours, a tolerance sweep, boolean failure/recovery, and long paths.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
use std::time::Instant;

use solidcore::curve::curve::Curve;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::p1::{BSplineCurve, NurbCircleParam};
use solidcore::primitive::primitive_creation::model::PrimitiveModel;
use solidcore::primitive::primitive_creation::p1::PrimitiveCreation;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::{Kernel, SolidError};

fn circle(radius: f64) -> Vec<Curve> {
    vec![Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
            radius,
            0.,
            360.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    )]
}

fn square(side: f64) -> Vec<Curve> {
    let h = side / 2.0;
    let c = [
        Vec3::new(-h, -h, 0.),
        Vec3::new(h, -h, 0.),
        Vec3::new(h, h, 0.),
        Vec3::new(-h, h, 0.),
        Vec3::new(-h, -h, 0.),
    ];
    c.windows(2).map(|w| line(w[0], w[1])).collect()
}

/// A square of `side` with a centered circular hole of `hole_radius`: two contours.
fn square_with_hole(side: f64, hole_radius: f64) -> Vec<Curve> {
    let mut curves = square(side);
    curves.extend(circle(hole_radius));
    curves
}

/// A closed square of `side` centered at `center` (in the XY plane).
fn square_at(center: Vec3, side: f64) -> Vec<Curve> {
    let h = side / 2.0;
    let p = |x: f64, y: f64| Vec3::new(center.x + x, center.y + y, 0.);
    let c = [p(-h, -h), p(h, -h), p(h, h), p(-h, h), p(-h, -h)];
    c.windows(2).map(|w| line(w[0], w[1])).collect()
}

/// Two disjoint squares side by side: two independent contours.
fn two_squares(side: f64, gap: f64) -> Vec<Curve> {
    let off = (side + gap) / 2.0;
    let mut curves = square_at(Vec3::new(-off, 0., 0.), side);
    curves.extend(square_at(Vec3::new(off, 0., 0.), side));
    curves
}

fn line(a: Vec3, b: Vec3) -> Curve {
    Curve::from(BSplineCurve::create_line_segment(3, a, b).unwrap())
}

/// A line, or `None` when the curve builder refuses it (e.g. a zero-length segment).
fn try_line(a: Vec3, b: Vec3) -> Option<Curve> {
    BSplineCurve::create_line_segment(3, a, b)
        .ok()
        .map(Curve::from)
}

fn arc_xy(center: Vec3, radius: f64, a0: f64, a1: f64) -> Curve {
    Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(center, Vec3::X, Vec3::Y),
            radius,
            a0,
            a1,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    )
}

fn build(
    k: &mut Kernel,
    profiles: &[Curve],
    paths: &[Curve],
    tol: f64,
    cap_start: bool,
    cap_end: bool,
) -> Result<solidcore::Brep, SolidError> {
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = PrimitiveCreation::create_sweep_along_planar_path_hybrid(
            model,
            PrimitiveModel::MAIN,
            profiles,
            paths,
            tol,
            cap_start,
            cap_end,
            true,
        )?;
        model.take_brep(id)
    })?;
    b.refresh_cached_topology();
    Ok(b)
}

fn build_union(
    k: &mut Kernel,
    profiles: &[Curve],
    paths: &[Curve],
    tol: f64,
) -> Result<solidcore::Brep, SolidError> {
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = PrimitiveCreation::create_sweep_along_planar_path_union(
            model,
            PrimitiveModel::MAIN,
            profiles,
            paths,
            tol,
            true,
            true,
            true,
        )?;
        model.take_brep(id)
    })?;
    b.refresh_cached_topology();
    Ok(b)
}

fn volume(k: &mut Kernel, b: &mut solidcore::Brep) -> Result<f64, SolidError> {
    let store = AttributeStore::default();
    // 1e-9: the mass integrator fails below this on some valid mixed-path bodies.
    brep_volume(
        k,
        b,
        1e-9,
        &PropertyAttributes {
            store: &store,
            key_of: &|_| ObjectKey(0),
        },
    )
}

/// Audit one case: `Err` is a clean refusal; `Ok` must be a valid manifold solid. A
/// `target` is reported (not asserted on refusal).
fn audit(
    name: &str,
    profiles: &[Curve],
    paths: &[Curve],
    tol: f64,
    cap_start: bool,
    cap_end: bool,
    target: Option<f64>,
) {
    let mut k = Kernel::default();
    let start = Instant::now();
    let built = build(&mut k, profiles, paths, tol, cap_start, cap_end);
    let millis = start.elapsed().as_millis();
    match built {
        Err(e) => eprintln!("{name:34} REFUSED ({millis} ms): {e:?}"),
        Ok(mut b) => {
            let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
            let manifold = is_manifold_solid(&b).unwrap_or(false);
            let v = volume(&mut k, &mut b);
            let faces = b.face_list().len();
            let vol = match &v {
                Ok(v) => format!("{v:.9}"),
                Err(e) => format!("<integrator {e:?}>"),
            };
            let err = target
                .zip(v.as_ref().ok().copied())
                .map(|(t, v)| (v - t) / t);
            eprintln!(
                "{name:34} valid={} manifold={manifold} faces={faces} volume={vol} err={} ({millis} ms)",
                check.is_valid(),
                err.map_or_else(|| "-".to_string(), |e| format!("{e:+.2e}")),
            );
            for info in check.lines().iter().take(3) {
                eprintln!("    ! {info}");
            }
            assert!(check.is_valid(), "{name}: returned an invalid body");
            assert!(manifold, "{name}: returned a non-manifold body");
        }
    }
}

#[test]
#[ignore = "diagnostic: near-tangent joint behaviour of reference vs union"]
fn probe_near_tangent_joints() {
    let r = 0.05;
    for deg in [0.1_f64, 0.5, 1.5, 2.5, 5.0] {
        let a = deg.to_radians();
        let paths = [
            line(Vec3::ZERO, Vec3::new(1., 0., 0.)),
            line(Vec3::new(1., 0., 0.), Vec3::new(1. + a.cos(), a.sin(), 0.)),
        ];
        for (label, union) in [("reference", false), ("union", true)] {
            let mut k = Kernel::default();
            let built = PrimitiveModel::build_new(&mut k, |_, model| {
                let id = if union {
                    PrimitiveCreation::create_sweep_along_planar_path_union(
                        model,
                        PrimitiveModel::MAIN,
                        &circle(r),
                        &paths,
                        1e-7,
                        true,
                        true,
                        true,
                    )?
                } else {
                    PrimitiveCreation::create_sweep_along_planar_path(
                        model,
                        PrimitiveModel::MAIN,
                        &circle(r),
                        &paths,
                        1e-7,
                        true,
                        true,
                        true,
                    )?
                };
                model.take_brep(id)
            });
            match built {
                Ok((_, b)) => {
                    let manifold = is_manifold_solid(&b).unwrap_or(false);
                    eprintln!("{deg:4} deg {label:9} manifold={manifold}");
                }
                Err(e) => eprintln!("{deg:4} deg {label:9} ERR {e:?}"),
            }
        }
    }
}

#[test]
fn near_tangent_and_short_segments() {
    let r = 0.05;
    // finite kink (5 degrees): the reference mitre handles it
    let a = 5.0_f64.to_radians();
    let near = [
        line(Vec3::ZERO, Vec3::new(1., 0., 0.)),
        line(Vec3::new(1., 0., 0.), Vec3::new(1. + a.cos(), a.sin(), 0.)),
    ];
    audit(
        "kink 5 deg",
        &circle(r),
        &near,
        1e-7,
        true,
        true,
        Some(std::f64::consts::PI * r * r * (1.0 + 1.0)),
    );

    // sub-threshold kink (0.5 degree): treated as tangent, no mitre
    let a = 0.5_f64.to_radians();
    let tiny = [
        line(Vec3::ZERO, Vec3::new(1., 0., 0.)),
        line(Vec3::new(1., 0., 0.), Vec3::new(1. + a.cos(), a.sin(), 0.)),
    ];
    audit(
        "kink 0.5 deg",
        &circle(r),
        &tiny,
        1e-7,
        true,
        true,
        Some(std::f64::consts::PI * r * r * 2.0),
    );

    // segment much shorter than the profile (L = 0.01, r = 0.05)
    let short = [
        line(Vec3::ZERO, Vec3::new(0.01, 0., 0.)),
        line(Vec3::new(0.01, 0., 0.), Vec3::new(0.01, 0.01, 0.)),
    ];
    audit(
        "short segments (L<r)",
        &circle(r),
        &short,
        1e-7,
        true,
        true,
        None,
    );

    // very short curved segment relative to the profile
    let a1 = 1.0_f64.to_radians();
    let arc_end = Vec3::new(a1.cos(), a1.sin(), 0.);
    let short_arc = [
        arc_xy(Vec3::ZERO, 1.0, 0.0, 1.0),
        line(arc_end, Vec3::new(arc_end.x, 1.0, 0.)),
    ];
    audit(
        "short arc segment",
        &circle(r),
        &short_arc,
        1e-7,
        true,
        true,
        None,
    );
}

#[test]
fn degenerate_segments() {
    let r = 0.05;
    // a zero-length segment: the curve builder refuses it before the sweep runs
    let p = Vec3::new(1., 0., 0.);
    match try_line(p, p) {
        None => eprintln!("zero-length segment: refused by the curve builder (clean)"),
        Some(zero) => {
            let zero = [line(Vec3::ZERO, p), zero, line(p, Vec3::new(1., 1., 0.))];
            audit(
                "zero-length segment",
                &circle(r),
                &zero,
                1e-7,
                true,
                true,
                None,
            );
        }
    }

    // a path whose segments do not connect: the union must refuse, not return a piece
    let disconnected = [
        line(Vec3::ZERO, Vec3::new(1., 0., 0.)),
        line(Vec3::new(5., 5., 0.), Vec3::new(6., 5., 0.)),
    ];
    audit(
        "disconnected path",
        &circle(r),
        &disconnected,
        1e-7,
        true,
        true,
        None,
    );

    // a single zero-length path: refused by the curve builder
    match try_line(Vec3::ZERO, Vec3::ZERO) {
        None => eprintln!("single degenerate path: refused by the curve builder (clean)"),
        Some(one) => audit(
            "single degenerate path",
            &circle(r),
            &[one],
            1e-7,
            true,
            true,
            None,
        ),
    }
}

#[test]
fn tolerance_sweep() {
    let r = 0.05;
    let arcline = [
        arc_xy(Vec3::ZERO, 5.0, 0.0, 90.0),
        line(Vec3::new(0., 5., 0.), Vec3::new(0., 10., 0.)),
    ];
    let len = 5.0 * std::f64::consts::FRAC_PI_2 + 5.0;
    let target = std::f64::consts::PI * r * r * len;
    for tol in [1e-8, 1e-7, 1e-6, 1e-5, 1e-4] {
        audit(
            &format!("arc+line tol={tol:e}"),
            &circle(r),
            &arcline,
            tol,
            true,
            true,
            Some(target),
        );
    }
}

#[test]
fn multi_contour_profiles() {
    let path = [line(Vec3::ZERO, Vec3::new(0., 0., 2.))];
    // square with a circular hole: two contours, one hole
    let side = 0.2;
    let hole = 0.05;
    let area = side * side - std::f64::consts::PI * hole * hole;
    audit(
        "square with a hole",
        &square_with_hole(side, hole),
        &path,
        1e-7,
        true,
        true,
        Some(area * 2.0),
    );
    // two disjoint square contours: must be refused or a valid body, never partial
    audit(
        "two disjoint contours",
        &two_squares(0.1, 0.05),
        &path,
        1e-7,
        true,
        true,
        None,
    );
}

#[test]
fn invalid_input_refused_and_kernel_recovers() {
    // A disconnected path is invalid input: the router must refuse it rather than let
    // the union return two disjoint shells as a non-manifold body. A later, valid build
    // on the same kernel must still succeed (no corrupted state). Boolean-step errors
    // are propagated by the mandatory `boolean_trees(..)?` in the union pipeline, whose
    // owned-body cleanup is covered by the reference's `audit_483` internal test.
    let mut k = Kernel::default();
    k.fixes.planar_path_sweep_statuses = true;
    let disconnected = [
        line(Vec3::ZERO, Vec3::new(1., 0., 0.)),
        line(Vec3::new(5., 5., 0.), Vec3::new(6., 5., 0.)),
    ];
    let failed = build(&mut k, &circle(0.05), &disconnected, 1e-7, true, true);
    assert!(failed.is_err(), "disconnected path must be refused");
    eprintln!("disconnected path refused: {:?}", failed.err());

    let elbow = [
        line(Vec3::ZERO, Vec3::new(0., 0., 1.)),
        line(Vec3::new(0., 0., 1.), Vec3::new(1., 0., 1.)),
    ];
    let mut b = build(&mut k, &circle(0.05), &elbow, 1e-7, true, true)
        .expect("kernel still usable after a refused boolean");
    assert!(is_manifold_solid(&b).unwrap_or(false));
    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    assert!(check.is_valid(), "{:?}", check.lines());
    eprintln!("recovery build after failure: valid");
}

/// A closed circle of `n` equal tangent-continuous arcs (a long curved path with many
/// segments).
fn closed_circle_arcs(n: usize, radius: f64) -> Vec<Curve> {
    let step = 360.0 / n as f64;
    (0..n)
        .map(|i| {
            let a0 = i as f64 * step;
            arc_xy(Vec3::ZERO, radius, a0, a0 + step)
        })
        .collect()
}

#[test]
fn performance_long_paths() {
    // long all-straight path (staircase of 90 degree kinks): routed to the reference.
    // Report the counts that build; a clean refusal is a valid audit outcome.
    for n in [8_usize, 16, 24] {
        let mut staircase = Vec::new();
        let mut p = Vec3::ZERO;
        for i in 0..n {
            let q = if i % 2 == 0 {
                p + Vec3::new(0.5, 0., 0.)
            } else {
                p + Vec3::new(0., 0.5, 0.)
            };
            staircase.push(line(p, q));
            p = q;
        }
        let start = Instant::now();
        let mut k = Kernel::default();
        let straight = build(&mut k, &circle(0.02), &staircase, 1e-7, true, true);
        let ms = start.elapsed().as_millis();
        match straight {
            Ok(b) => {
                assert!(is_manifold_solid(&b).unwrap_or(false));
                eprintln!("long straight path: {n} segments hybrid OK in {ms} ms");
            }
            Err(_) => {
                eprintln!("long straight path: {n} segments hybrid REFUSED in {ms} ms");
                // is the union a robust fallback for the same path?
                let start = Instant::now();
                let mut k = Kernel::default();
                let union = build_union(&mut k, &circle(0.02), &staircase, 1e-7);
                let ums = start.elapsed().as_millis();
                match union {
                    Ok(b) => {
                        let ok = is_manifold_solid(&b).unwrap_or(false);
                        eprintln!(
                            "  same path via union: {n} segments OK manifold={ok} in {ums} ms"
                        );
                    }
                    Err(e) => {
                        eprintln!("  same path via union: {n} segments REFUSED in {ums} ms ({e:?})")
                    }
                }
            }
        }
    }

    // long curved path: a closed circle split into more and more tangent-continuous
    // arcs, routed to the union, to show how the boolean union scales with the count
    for segments in [4_usize, 8, 12] {
        let path = closed_circle_arcs(segments, 5.0);
        let start = Instant::now();
        let mut k = Kernel::default();
        let curved = build(&mut k, &circle(0.02), &path, 1e-7, true, true);
        let ms = start.elapsed().as_millis();
        match curved {
            Ok(b) => {
                assert!(is_manifold_solid(&b).unwrap_or(false));
                eprintln!("curved closed path: {segments} segments OK in {ms} ms");
            }
            Err(e) => {
                eprintln!("curved closed path: {segments} segments REFUSED in {ms} ms ({e:?})");
            }
        }
    }
}
