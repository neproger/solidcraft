//! Matrix for `create_sweep_along_planar_path_union` (union of nominal closed segment
//! sweeps). This is a separate algorithm from the analytic miter and is measured on its
//! own: validity, manifold-solid, volume against the smooth-pipe volume `area * path
//! length`, face count, and build time. The measured deficit is the overlap lens near
//! each joint (the two perpendicular end caps cross; the union counts the lens once).
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

/// A closed square profile of side `side`, centered on the origin.
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

/// A closed rounded rectangle profile centered on the origin: four lines and four
/// quarter-circle corners (a profile with several roundings).
fn rounded_rectangle(width: f64, height: f64, radius: f64) -> Vec<Curve> {
    let (hw, hh) = (width / 2.0, height / 2.0);
    let p = |x: f64, y: f64| Vec3::new(x, y, 0.);
    let arc = |cx: f64, cy: f64, a0: f64, a1: f64| {
        Curve::from(
            BSplineCurve::create_circle_segment(
                3,
                &Axis2Placement::new(p(cx, cy), Vec3::X, Vec3::Y),
                radius,
                a0,
                a1,
                NurbCircleParam::Quadratic,
            )
            .unwrap(),
        )
    };
    vec![
        line(p(-hw + radius, -hh), p(hw - radius, -hh)),
        arc(hw - radius, -hh + radius, 270., 360.),
        line(p(hw, -hh + radius), p(hw, hh - radius)),
        arc(hw - radius, hh - radius, 0., 90.),
        line(p(hw - radius, hh), p(-hw + radius, hh)),
        arc(-hw + radius, hh - radius, 90., 180.),
        line(p(-hw, hh - radius), p(-hw, -hh + radius)),
        arc(-hw + radius, -hh + radius, 180., 270.),
    ]
}

fn line(a: Vec3, b: Vec3) -> Curve {
    Curve::from(BSplineCurve::create_line_segment(3, a, b).unwrap())
}

/// A circle arc in the XY plane through `center` of `radius`, degrees `a0..a1`.
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
    cap_start: bool,
    cap_end: bool,
) -> Result<solidcore::Brep, SolidError> {
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = PrimitiveCreation::create_sweep_along_planar_path_union(
            model,
            PrimitiveModel::MAIN,
            profiles,
            paths,
            1e-7,
            cap_start,
            cap_end,
            true,
        )?;
        model.take_brep(id)
    })?;
    b.refresh_cached_topology();
    Ok(b)
}

fn volume(k: &mut Kernel, b: &mut solidcore::Brep) -> f64 {
    let store = AttributeStore::default();
    brep_volume(
        k,
        b,
        1e-12,
        &PropertyAttributes {
            store: &store,
            key_of: &|_| ObjectKey(0),
        },
    )
    .unwrap()
}

#[derive(Default)]
struct Outcome {
    ok: bool,
    manifold: bool,
    valid: bool,
    error: f64,
    deficit: f64,
    volume: f64,
    faces: usize,
    millis: u128,
}

fn run(name: &str, profiles: &[Curve], paths: &[Curve], target: f64) -> Outcome {
    let mut k = Kernel::default();
    let start = Instant::now();
    let built = build(&mut k, profiles, paths, true, true);
    let millis = start.elapsed().as_millis();
    let mut outcome = Outcome {
        millis,
        ..Outcome::default()
    };
    let Ok(mut b) = built else {
        eprintln!("{name:24} BUILD FAILED");
        return outcome;
    };
    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    outcome.ok = true;
    outcome.valid = check.is_valid();
    outcome.manifold = is_manifold_solid(&b).unwrap_or(false);
    outcome.faces = b.face_list().len();
    let v = volume(&mut k, &mut b);
    outcome.error = (v - target) / target;
    outcome.deficit = target - v;
    outcome.volume = v;
    eprintln!(
        "{name:24} valid={} manifold={} faces={:3} volume={v:.9} target={target:.9} err={:+.3e} deficit={:.9} t={millis}ms",
        outcome.valid, outcome.manifold, outcome.faces, outcome.error, outcome.deficit,
    );
    if !outcome.valid {
        for info in check.lines().iter().take(3) {
            eprintln!("    ! {info}");
        }
    }
    outcome
}

#[test]
fn union_sweep_matrix() {
    let r = 0.05;
    let mut results: Vec<(&str, Outcome)> = Vec::new();

    // 1. straight + straight (right elbow)
    let elbow_paths = [
        line(Vec3::ZERO, Vec3::new(0., 0., 1.)),
        line(Vec3::new(0., 0., 1.), Vec3::new(1., 0., 1.)),
    ];
    let elbow_target = std::f64::consts::PI * 0.1 * 0.1 * 2.0;
    results.push((
        "line+line",
        run("line+line", &circle(0.1), &elbow_paths, elbow_target),
    ));

    // 2. arc + line (the case where the reference is invalid)
    let arcline_paths = [
        arc_xy(Vec3::ZERO, 5.0, 0.0, 90.0),
        line(Vec3::new(0., 5., 0.), Vec3::new(0., 10., 0.)),
    ];
    let arcline_len = 5.0 * std::f64::consts::FRAC_PI_2 + 5.0;
    results.push((
        "arc+line",
        run(
            "arc+line",
            &circle(r),
            &arcline_paths,
            std::f64::consts::PI * r * r * arcline_len,
        ),
    ));

    // 3. arc + arc (both curved, a 90 degree joint)
    let arcarc_paths = [
        arc_xy(Vec3::ZERO, 5.0, 0.0, 90.0),
        arc_xy(Vec3::new(3., 5., 0.), 3.0, 180.0, 270.0),
    ];
    let arcarc_len = 5.0 * std::f64::consts::FRAC_PI_2 + 3.0 * std::f64::consts::FRAC_PI_2;
    results.push((
        "arc+arc",
        run(
            "arc+arc",
            &circle(r),
            &arcarc_paths,
            std::f64::consts::PI * r * r * arcarc_len,
        ),
    ));

    // 4. several curved segments: line, arc, line
    let multi_paths = [
        line(Vec3::ZERO, Vec3::new(2., 0., 0.)),
        arc_xy(Vec3::ZERO, 2.0, 0.0, 90.0),
        line(Vec3::new(0., 2., 0.), Vec3::new(0., 4., 0.)),
    ];
    let multi_len = 2.0 + 2.0 * std::f64::consts::FRAC_PI_2 + 2.0;
    results.push((
        "line+arc+line",
        run(
            "line+arc+line",
            &circle(r),
            &multi_paths,
            std::f64::consts::PI * r * r * multi_len,
        ),
    ));

    // 5. closed path (square ring)
    let sq = 2.0;
    let ring_paths = [
        line(Vec3::ZERO, Vec3::new(sq, 0., 0.)),
        line(Vec3::new(sq, 0., 0.), Vec3::new(sq, sq, 0.)),
        line(Vec3::new(sq, sq, 0.), Vec3::new(0., sq, 0.)),
        line(Vec3::new(0., sq, 0.), Vec3::ZERO),
    ];
    let ring_target = std::f64::consts::PI * r * r * (4.0 * sq);
    results.push((
        "closed square ring",
        run("closed square ring", &circle(r), &ring_paths, ring_target),
    ));

    // 6. square profile (non-circular) on arc + line
    let square_target = 0.1 * 0.1 * arcline_len;
    results.push((
        "arc+line square",
        run(
            "arc+line square",
            &square(0.1),
            &arcline_paths,
            square_target,
        ),
    ));

    // 7. rounded-rectangle profile (several roundings) on arc + arc
    let rr_area = 0.2 * 0.1 - (4.0 - std::f64::consts::PI) * 0.02 * 0.02;
    results.push((
        "arc+arc rounded-rect",
        run(
            "arc+arc rounded-rect",
            &rounded_rectangle(0.2, 0.1, 0.02),
            &arcarc_paths,
            rr_area * arcarc_len,
        ),
    ));

    // Summary and the matrix's pass condition: every case builds, is a valid
    // manifold solid, and its volume is within the overlap-lens budget of the
    // smooth-pipe volume.
    eprintln!("---");
    for (name, o) in &results {
        eprintln!(
            "{name:24} build={} valid={} manifold={} err={:+.3e} t={}ms",
            o.ok, o.valid, o.manifold, o.error, o.millis
        );
    }
    for (name, o) in &results {
        assert!(o.ok, "{name}: build failed");
        assert!(o.valid, "{name}: invalid body");
        assert!(o.manifold, "{name}: not a manifold solid");
    }

    // The volume deficit is the overlap lens where the perpendicular end caps cross:
    // exactly (4/3) r^3 per joint for a circular profile, to under 1%, independent of
    // path length and curvature. This rules out a tolerance or segment-build error.
    let lens = |radius: f64| 4.0 / 3.0 * radius * radius * radius;
    let check_lens = |outcome: &Outcome, joints: usize, radius: f64| {
        let expected = joints as f64 * lens(radius);
        let rel = (outcome.deficit - expected).abs() / expected;
        eprintln!(
            "lens check: {joints} joint(s) r={radius}: deficit={:.9} expected={expected:.9} rel={rel:.2e}",
            outcome.deficit
        );
        assert!(rel < 1e-2, "deficit does not match the lens law: {rel:.3e}");
    };
    check_lens(&results[0].1, 1, 0.1);
    check_lens(&results[1].1, 1, r);
    check_lens(&results[2].1, 1, r);
    check_lens(&results[3].1, 2, r);
    check_lens(&results[4].1, 4, r);

    // Boolean stability: the same build twice gives the same body.
    let mut k2 = Kernel::default();
    let mut again = build(&mut k2, &circle(r), &arcline_paths, true, true).unwrap();
    let again_volume = volume(&mut k2, &mut again);
    assert!(
        (again_volume - results[1].1.volume).abs() <= 1e-12,
        "repeat build volume changed: {again_volume} vs {}",
        results[1].1.volume
    );
    assert_eq!(
        again.face_list().len(),
        results[1].1.faces,
        "repeat build faces differ"
    );
}
