//! Extended matrix for `create_sweep_along_planar_path_union`: closed paths of arcs and
//! of arcs mixed with lines, closed and open paths with several profile loops, and the
//! four `cap_start` / `cap_end` combinations. Separate tests so a case group can be run
//! on its own; the volume target is the smooth-pipe volume `area * path length`, and the
//! only expected deficit is the kinked-joint lens (none for tangent-continuous loops).
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
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

/// Two concentric circles: a hollow annular section (two profile loops).
fn annulus(outer: f64, inner: f64) -> Vec<Curve> {
    let ring = |radius: f64| {
        Curve::from(
            BSplineCurve::create_circle_segment(
                3,
                &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
                radius,
                0.,
                360.,
                NurbCircleParam::Quadratic,
            )
            .unwrap(),
        )
    };
    vec![ring(outer), ring(inner)]
}

fn line(a: Vec3, b: Vec3) -> Curve {
    Curve::from(BSplineCurve::create_line_segment(3, a, b).unwrap())
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

/// Build, validate and report one case; `expect_solid` is the expected manifold status
/// (`false` for an open outer end), and a solid is checked against `target`.
fn run_case(
    name: &str,
    profiles: &[Curve],
    paths: &[Curve],
    cap_start: bool,
    cap_end: bool,
    target: Option<f64>,
    expect_solid: bool,
) {
    let mut k = Kernel::default();
    let mut b = match build(&mut k, profiles, paths, cap_start, cap_end) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("{name:28} BUILD FAILED: {e:?}");
            panic!("{name} build failed");
        }
    };
    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    let manifold = is_manifold_solid(&b).unwrap_or(false);
    let v = volume(&mut k, &mut b);
    let faces = b.face_list().len();
    match target {
        Some(target) => {
            let err = (v - target) / target;
            eprintln!(
                "{name:28} valid={} manifold={manifold} faces={faces:3} volume={v:.9} target={target:.9} err={err:+.3e}",
                check.is_valid(),
            );
        }
        None => {
            eprintln!(
                "{name:28} valid={} manifold={manifold} faces={faces:3} volume={v:.9}",
                check.is_valid(),
            );
        }
    }
    if !check.is_valid() {
        for info in check.lines().iter().take(3) {
            eprintln!("    ! {info}");
        }
    }
    assert_eq!(manifold, expect_solid, "{name}: manifold status");
    if expect_solid {
        if let Some(target) = target {
            let err = (v - target) / target;
            assert!(err.abs() < 2e-2, "{name}: volume error {err:.3e}");
        }
    }
}

/// A closed loop of four quarter arcs (a circle of radius `radius`).
fn closed_arc_path(radius: f64) -> Vec<Curve> {
    [(0.0, 90.0), (90.0, 180.0), (180.0, 270.0), (270.0, 360.0)]
        .map(|(a0, a1)| arc_xy(Vec3::ZERO, radius, a0, a1))
        .to_vec()
}

/// A closed stadium: two straight sides and two semicircular ends.
fn closed_stadium_path() -> Vec<Curve> {
    vec![
        line(Vec3::new(-2., -1., 0.), Vec3::new(2., -1., 0.)),
        arc_xy(Vec3::new(2., 0., 0.), 1., -90., 90.),
        line(Vec3::new(2., 1., 0.), Vec3::new(-2., 1., 0.)),
        arc_xy(Vec3::new(-2., 0., 0.), 1., 90., 270.),
    ]
}

#[test]
fn closed_paths_of_arcs_and_mixed() {
    let r = 0.05;
    let ring = 5.0;
    let circle_len = 2.0 * std::f64::consts::PI * ring;
    run_case(
        "closed arcs (circle)",
        &circle(r),
        &closed_arc_path(ring),
        true,
        true,
        Some(std::f64::consts::PI * r * r * circle_len),
        true,
    );

    let stadium_len = 8.0 + 2.0 * std::f64::consts::PI;
    run_case(
        "closed mixed (stadium)",
        &circle(r),
        &closed_stadium_path(),
        true,
        true,
        Some(std::f64::consts::PI * r * r * stadium_len),
        true,
    );
}

#[test]
fn closed_curved_path_with_several_profiles() {
    let ring = 5.0;
    let circle_len = 2.0 * std::f64::consts::PI * ring;
    let area = std::f64::consts::PI * (0.06 * 0.06 - 0.04 * 0.04);
    run_case(
        "closed arcs, annulus",
        &annulus(0.06, 0.04),
        &closed_arc_path(ring),
        true,
        true,
        Some(area * circle_len),
        true,
    );
}

#[test]
fn open_path_with_several_profiles() {
    let arcline_paths = [
        arc_xy(Vec3::ZERO, 5.0, 0.0, 90.0),
        line(Vec3::new(0., 5., 0.), Vec3::new(0., 10., 0.)),
    ];
    let arcline_len = 5.0 * std::f64::consts::FRAC_PI_2 + 5.0;
    let area = std::f64::consts::PI * (0.06 * 0.06 - 0.04 * 0.04);
    run_case(
        "open arc+line, annulus",
        &annulus(0.06, 0.04),
        &arcline_paths,
        true,
        true,
        Some(area * arcline_len),
        true,
    );
}

#[test]
fn open_path_cap_combinations() {
    let paths = [
        line(Vec3::ZERO, Vec3::new(0., 0., 1.)),
        line(Vec3::new(0., 0., 1.), Vec3::new(1., 0., 1.)),
    ];
    for (cap_start, cap_end) in [(true, true), (true, false), (false, true), (false, false)] {
        let solid = cap_start && cap_end;
        let name = format!("caps({cap_start},{cap_end})");
        run_case(&name, &circle(0.1), &paths, cap_start, cap_end, None, solid);
    }
}
