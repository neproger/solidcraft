//! The analytic-mitre planar-path sweep against the reference on the straight 90°
//! elbow (the shape of `issue_2794_planar_path_mitred_elbow`). The two paths share the
//! whole pipeline and differ only in the cutting curve; for straight segments the
//! analytic miter plane is the exact bisector, so the results must agree.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
use solidcore::curve::curve::Curve;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::p1::{BSplineCurve, NurbCircleParam};
use solidcore::primitive::primitive_creation::model::PrimitiveModel;
use solidcore::primitive::primitive_creation::p1::PrimitiveCreation;
use solidcore::primitive::primitive_creation::p2::PrimitiveCreationAccess;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::Kernel;

#[derive(Clone, Copy)]
enum Miter {
    Reference,
    Analytic,
}

fn elbow(k: &mut Kernel, miter: Miter, radius: f64, leg: f64, tol: f64) -> solidcore::Brep {
    let profile = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
            radius,
            0.,
            360.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    let up = Vec3::new(0., 0., leg);
    let path = |a, b| Curve::from(BSplineCurve::create_line_segment(3, a, b).unwrap());
    let paths = [path(Vec3::ZERO, up), path(up, Vec3::new(leg, 0., leg))];
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = match miter {
            Miter::Reference => PrimitiveCreation::create_sweep_along_planar_path(
                model,
                PrimitiveModel::MAIN,
                std::slice::from_ref(&profile),
                &paths,
                tol,
                true,
                true,
                true,
            )?,
            Miter::Analytic => PrimitiveCreation::create_sweep_along_planar_path_analytic_miter(
                model,
                PrimitiveModel::MAIN,
                std::slice::from_ref(&profile),
                &paths,
                tol,
                true,
                true,
                true,
            )?,
        };
        model.take_brep(id)
    })
    .unwrap();
    b.refresh_cached_topology();
    b
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

fn square_ring(k: &mut Kernel, miter: Miter, radius: f64, side: f64, tol: f64) -> solidcore::Brep {
    let profile = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
            radius,
            0.,
            360.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    let c = [
        Vec3::ZERO,
        Vec3::new(side, 0., 0.),
        Vec3::new(side, side, 0.),
        Vec3::new(0., side, 0.),
        Vec3::ZERO,
    ];
    let paths: Vec<Curve> = c
        .windows(2)
        .map(|w| Curve::from(BSplineCurve::create_line_segment(3, w[0], w[1]).unwrap()))
        .collect();
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = match miter {
            Miter::Reference => PrimitiveCreation::create_sweep_along_planar_path(
                model,
                PrimitiveModel::MAIN,
                std::slice::from_ref(&profile),
                &paths,
                tol,
                false,
                false,
                true,
            )?,
            Miter::Analytic => PrimitiveCreation::create_sweep_along_planar_path_analytic_miter(
                model,
                PrimitiveModel::MAIN,
                std::slice::from_ref(&profile),
                &paths,
                tol,
                false,
                false,
                true,
            )?,
        };
        model.take_brep(id)
    })
    .unwrap();
    b.refresh_cached_topology();
    b
}

#[test]
fn analytic_miter_matches_reference_on_a_closed_square_ring() {
    let (radius, side, tol) = (1e-1, 1., 1e-7);
    let mut k = Kernel::default();

    let mut reference = square_ring(&mut k, Miter::Reference, radius, side, tol);
    let mut analytic = square_ring(&mut k, Miter::Analytic, radius, side, tol);

    let reference_check = check_brep(&mut k, &mut reference, CheckLevel::Level2).unwrap();
    assert!(reference_check.is_valid(), "{:?}", reference_check.lines());
    let analytic_check = check_brep(&mut k, &mut analytic, CheckLevel::Level2).unwrap();
    assert!(analytic_check.is_valid(), "{:?}", analytic_check.lines());

    let reference_volume = volume(&mut k, &mut reference);
    let analytic_volume = volume(&mut k, &mut analytic);
    assert!(
        (reference_volume - analytic_volume).abs() <= reference_volume.abs() * 1e-9,
        "{reference_volume} != {analytic_volume}"
    );
}

// Matrix diagnostic for a curved joint. The target is the smooth pipe volume
// `pi r^2 (arc length + line length)`; the reference is invalid here and the analytic
// path computes the true surface-surface seam, but neither closes a solid. The reason
// is measured, not guessed: the seam is a small non-planar loop (two 90-degree meridian
// arcs meeting at cusped poles), while `project_and_trim` only realizes a ruled/plane
// cut and classifies the deleted side against that plane. A plane cuts the torus and the
// cylinder in two different cross-sections, so their cut faces cannot abut. Realizing the
// seam needs an imprint-and-loop trim (UV trim curves merged on each face, side chosen by
// the seam loop rather than a plane), which is a change to the trim primitive, not to the
// cutter. Kept as an executable record until that lands.
#[test]
#[ignore = "blocker: project_and_trim is plane/ruled-only; a non-planar seam needs an imprint-and-loop trim"]
fn analytic_miter_builds_on_a_line_then_arc_path() {
    let radius = 0.05;
    let tol = 1e-7;
    let profile = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
            radius,
            0.,
            360.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    // a 90-degree arc in the XY plane, then a line leaving its end perpendicular: the
    // curved case the reference resolves by a two-point intersection chord
    let arc = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::IDENTITY,
            5.,
            0.,
            90.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    let line =
        Curve::from(BSplineCurve::create_line_segment(3, Vec3::Y * 5., Vec3::Y * 10.).unwrap());
    let paths = [arc, line];
    let mut k = Kernel::default();
    let build = |k: &mut Kernel, miter: Miter| {
        let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
            let id = match miter {
                Miter::Reference => PrimitiveCreation::create_sweep_along_planar_path(
                    model,
                    PrimitiveModel::MAIN,
                    std::slice::from_ref(&profile),
                    &paths,
                    tol,
                    true,
                    true,
                    true,
                )?,
                Miter::Analytic => {
                    PrimitiveCreation::create_sweep_along_planar_path_analytic_miter(
                        model,
                        PrimitiveModel::MAIN,
                        std::slice::from_ref(&profile),
                        &paths,
                        tol,
                        true,
                        true,
                        true,
                    )?
                }
            };
            model.take_brep(id)
        })
        .unwrap();
        b.refresh_cached_topology();
        b
    };

    let mut reference = build(&mut k, Miter::Reference);
    let mut analytic = build(&mut k, Miter::Analytic);
    let rc = check_brep(&mut k, &mut reference, CheckLevel::Level2).unwrap();
    let ac = check_brep(&mut k, &mut analytic, CheckLevel::Level2).unwrap();
    let rv = volume(&mut k, &mut reference);
    let av = volume(&mut k, &mut analytic);
    let rm = solidcore::topology::brep::p2::is_manifold_solid(&reference).unwrap_or(false);
    let am = solidcore::topology::brep::p2::is_manifold_solid(&analytic).unwrap_or(false);
    let target = std::f64::consts::PI * radius * radius * (5.0 * std::f64::consts::FRAC_PI_2 + 5.0);
    eprintln!(
        "line+arc: target volume={target:.9}\n  \
         reference valid={} manifold={rm} volume={rv:.9} err={:.3e} faces={}\n  \
         analytic  valid={} manifold={am} volume={av:.9} err={:.3e} faces={}",
        rc.is_valid(),
        (rv - target).abs() / target,
        reference.face_list().len(),
        ac.is_valid(),
        (av - target).abs() / target,
        analytic.face_list().len(),
    );
    let _ = (rc.lines(), ac.lines());
}

// Candidate next algorithm, measured here only: skip the mitre trim and let the final
// boolean union of the two nominal capped tubes form the joint. This sidesteps the
// plane-cut mismatch, at the cost of the small lens where the two tubes overlap.
#[test]
fn union_of_two_nominal_capped_tubes_is_a_solid() {
    let radius = 0.05;
    let tol = 1e-7;
    let profile = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
            radius,
            0.,
            360.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    let arc = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::IDENTITY,
            5.,
            0.,
            90.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    let line =
        Curve::from(BSplineCurve::create_line_segment(3, Vec3::Y * 5., Vec3::Y * 10.).unwrap());
    let mut k = Kernel::default();
    let built = PrimitiveModel::build_new(&mut k, |_, model| {
        let a = PrimitiveCreation::create_sweep_along_planar_path(
            model,
            PrimitiveModel::MAIN,
            std::slice::from_ref(&profile),
            std::slice::from_ref(&arc),
            tol,
            true,
            true,
            true,
        )?;
        let b = PrimitiveCreation::create_sweep_along_planar_path(
            model,
            PrimitiveModel::MAIN,
            std::slice::from_ref(&profile),
            std::slice::from_ref(&line),
            tol,
            true,
            true,
            true,
        )?;
        let mut ids = vec![a, b];
        model.boolean_trees(&mut ids, &[0, 1, -1])?;
        let id = ids.into_iter().next().unwrap();
        model.take_brep(id)
    });
    let (_, mut b) = built.unwrap();
    b.refresh_cached_topology();
    let manifold = solidcore::topology::brep::p2::is_manifold_solid(&b).unwrap_or(false);
    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    let v = volume(&mut k, &mut b);
    let target = std::f64::consts::PI * radius * radius * (5.0 * std::f64::consts::FRAC_PI_2 + 5.0);
    let err = (v - target).abs() / target;
    eprintln!(
        "union-of-nominal: manifold={manifold} valid={} volume={v:.9} target={target:.9} err={err:.3e}",
        check.is_valid(),
    );
    assert!(check.is_valid(), "{:?}", check.lines());
    assert!(manifold, "union of nominal capped tubes is not a solid");
    // the deficit is the lens the two tubes share, counted once instead of the
    // smooth pipe's redistribution; well under the profile size
    assert!(err < 5e-3, "volume error {err:.3e} too large");
}

#[test]
fn a_single_arc_sweep_is_a_solid() {
    let radius = 0.05;
    let tol = 1e-7;
    let profile = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
            radius,
            0.,
            360.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    let arc = Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::IDENTITY,
            5.,
            0.,
            90.,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    );
    let mut k = Kernel::default();
    let (_, mut b) = PrimitiveModel::build_new(&mut k, |_, model| {
        let id = PrimitiveCreation::create_sweep_along_planar_path(
            model,
            PrimitiveModel::MAIN,
            &[profile],
            &[arc],
            tol,
            true,
            true,
            true,
        )?;
        model.take_brep(id)
    })
    .unwrap();
    b.refresh_cached_topology();
    let manifold = solidcore::topology::brep::p2::is_manifold_solid(&b).unwrap_or(false);
    let check = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    let v = volume(&mut k, &mut b);
    let exact = std::f64::consts::PI * radius * radius * (5.0 * std::f64::consts::FRAC_PI_2);
    eprintln!(
        "single arc: manifold={manifold} valid={} volume={v:.9} exact={exact:.9} err={:.3e}",
        check.is_valid(),
        (v - exact).abs() / exact
    );
    assert!(
        manifold,
        "single arc sweep is not a solid: {:?}",
        check.lines()
    );
}

#[test]
fn analytic_miter_matches_reference_on_a_straight_elbow() {
    let (radius, leg, tol) = (1e-1, 1., 1e-7);
    let mut k = Kernel::default();

    let mut reference = elbow(&mut k, Miter::Reference, radius, leg, tol);
    let mut analytic = elbow(&mut k, Miter::Analytic, radius, leg, tol);

    let reference_check = check_brep(&mut k, &mut reference, CheckLevel::Level2).unwrap();
    assert!(reference_check.is_valid(), "{:?}", reference_check.lines());
    let analytic_check = check_brep(&mut k, &mut analytic, CheckLevel::Level2).unwrap();
    assert!(analytic_check.is_valid(), "{:?}", analytic_check.lines());

    let reference_volume = volume(&mut k, &mut reference);
    let analytic_volume = volume(&mut k, &mut analytic);
    assert!(
        (reference_volume - analytic_volume).abs() <= reference_volume.abs() * 1e-9,
        "{reference_volume} != {analytic_volume}"
    );

    let exact = std::f64::consts::PI * radius * radius * 2. * leg;
    assert!(
        (analytic_volume - exact).abs() < exact * 1e-4,
        "{analytic_volume} != {exact}"
    );

    assert_eq!(
        reference.face_list().len(),
        analytic.face_list().len(),
        "face counts differ"
    );
}
