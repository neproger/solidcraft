//! Status: historical comparison of the three prior planar-path implementations. Not the
//! AOIM pilot; see `plan/analytical-miter.md`.
//!
//! Side-by-side comparison of the three planar-path implementations, kept as separate
//! algorithms:
//!
//! - `reference` тАФ `create_sweep_along_planar_path` (the untouched control);
//! - `analytic` тАФ `create_sweep_along_planar_path_analytic_miter` (seam experiment);
//! - `union` тАФ `create_sweep_along_planar_path_union` (nominal closed segments, boolean
//!   union).
//!
//! Each row reports validity, manifold-solid, volume and build time. The volume target is
//! the smooth-pipe volume `area * path length`; only the reference (straight segments)
//! and the union (minus its lens) are expected to be close to it.
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

#[derive(Clone, Copy)]
enum Implementation {
    Reference,
    Analytic,
    Union,
}

impl Implementation {
    fn label(self) -> &'static str {
        match self {
            Implementation::Reference => "reference",
            Implementation::Analytic => "analytic",
            Implementation::Union => "union",
        }
    }
}

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
    implementation: Implementation,
    profiles: &[Curve],
    paths: &[Curve],
) -> Result<solidcore::Brep, SolidError> {
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = match implementation {
            Implementation::Reference => PrimitiveCreation::create_sweep_along_planar_path(
                model,
                PrimitiveModel::MAIN,
                profiles,
                paths,
                1e-7,
                true,
                true,
                true,
            )?,
            Implementation::Analytic => {
                PrimitiveCreation::create_sweep_along_planar_path_analytic_miter(
                    model,
                    PrimitiveModel::MAIN,
                    profiles,
                    paths,
                    1e-7,
                    true,
                    true,
                    true,
                )?
            }
            Implementation::Union => PrimitiveCreation::create_sweep_along_planar_path_union(
                model,
                PrimitiveModel::MAIN,
                profiles,
                paths,
                1e-7,
                true,
                true,
                true,
            )?,
        };
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

fn compare(case: &str, profiles: &[Curve], paths: &[Curve], target: f64) {
    eprintln!("case {case} (target {target:.9})");
    for implementation in [
        Implementation::Reference,
        Implementation::Analytic,
        Implementation::Union,
    ] {
        let mut k = Kernel::default();
        let start = Instant::now();
        let built = build(&mut k, implementation, profiles, paths);
        let millis = start.elapsed().as_millis();
        let Ok(mut b) = built else {
            eprintln!(
                "  {:<10} BUILD FAILED  ({millis} ms)",
                implementation.label()
            );
            continue;
        };
        let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
        let manifold = is_manifold_solid(&b).unwrap_or(false);
        let v = volume(&mut k, &mut b);
        eprintln!(
            "  {:<10} valid={:<5} manifold={:<5} faces={:3} volume={v:.9} err={:+.3e} ({millis} ms)",
            implementation.label(),
            check.is_valid(),
            manifold,
            b.face_list().len(),
            (v - target) / target,
        );
    }
}

#[test]
fn compare_three_implementations() {
    // straight 90 degree elbow
    let elbow = [
        line(Vec3::ZERO, Vec3::new(0., 0., 1.)),
        line(Vec3::new(0., 0., 1.), Vec3::new(1., 0., 1.)),
    ];
    compare(
        "straight elbow",
        &circle(0.1),
        &elbow,
        std::f64::consts::PI * 0.1 * 0.1 * 2.0,
    );

    // arc + line (the reference's invalid case)
    let arcline = [
        arc_xy(Vec3::ZERO, 5.0, 0.0, 90.0),
        line(Vec3::new(0., 5., 0.), Vec3::new(0., 10., 0.)),
    ];
    let arcline_len = 5.0 * std::f64::consts::FRAC_PI_2 + 5.0;
    compare(
        "arc + line",
        &circle(0.05),
        &arcline,
        std::f64::consts::PI * 0.05 * 0.05 * arcline_len,
    );

    // closed square ring
    let sq = 2.0;
    let ring = [
        line(Vec3::ZERO, Vec3::new(sq, 0., 0.)),
        line(Vec3::new(sq, 0., 0.), Vec3::new(sq, sq, 0.)),
        line(Vec3::new(sq, sq, 0.), Vec3::new(0., sq, 0.)),
        line(Vec3::new(0., sq, 0.), Vec3::ZERO),
    ];
    compare(
        "closed square ring",
        &circle(0.05),
        &ring,
        std::f64::consts::PI * 0.05 * 0.05 * (4.0 * sq),
    );
}
