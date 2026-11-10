//! The 90-degree elbow, measured step by step to establish whether the union's volume
//! loss is a construction/boolean defect or the geometry of the construction:
//!
//! 1. each nominal segment is a capped cylinder of volume `pi r^2 L` (construction and
//!    cap orientation correct);
//! 2. the reference mitre gives `2 pi r^2 L` (the smooth pipe);
//! 3. the union gives `2 pi r^2 L - V_overlap`, and `V_overlap = V(A) + V(B) - V(union)`;
//! 4. `V_overlap` equals the closed form `4/3 r^3` of two perpendicular capped cylinders,
//!    so the deficit is the construction's geometry, not a boolean error.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
use solidcore::curve::curve::Curve;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::p1::{BSplineCurve, NurbCircleParam};
use solidcore::primitive::primitive_creation::model::PrimitiveModel;
use solidcore::primitive::primitive_creation::p1::PrimitiveCreation;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::{Kernel, SolidError};

const R: f64 = 0.1;
const L: f64 = 1.0;

fn circle() -> Vec<Curve> {
    vec![Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::new(Vec3::ZERO, Vec3::X, Vec3::Y),
            R,
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

fn volume(k: &mut Kernel, b: &mut solidcore::Brep) -> f64 {
    let store = AttributeStore::default();
    brep_volume(
        k,
        b,
        1e-9,
        &PropertyAttributes {
            store: &store,
            key_of: &|_| ObjectKey(0),
        },
    )
    .unwrap()
}

fn build_union(k: &mut Kernel, paths: &[Curve]) -> Result<solidcore::Brep, SolidError> {
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = PrimitiveCreation::create_sweep_along_planar_path_union(
            model,
            PrimitiveModel::MAIN,
            &circle(),
            paths,
            1e-7,
            true,
            true,
            true,
        )?;
        model.take_brep(id)
    })?;
    b.refresh_cached_topology();
    Ok(b)
}

fn build_reference(k: &mut Kernel, paths: &[Curve]) -> Result<solidcore::Brep, SolidError> {
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = PrimitiveCreation::create_sweep_along_planar_path(
            model,
            PrimitiveModel::MAIN,
            &circle(),
            paths,
            1e-7,
            true,
            true,
            true,
        )?;
        model.take_brep(id)
    })?;
    b.refresh_cached_topology();
    Ok(b)
}

#[test]
fn elbow_volume_loss_is_the_exact_cap_overlap() {
    let up = Vec3::new(0., 0., L);
    let over = Vec3::new(L, 0., L);
    let segment_a = [line(Vec3::ZERO, up)];
    let segment_b = [line(up, over)];

    let mut k = Kernel::default();
    let mut a = build_union(&mut k, &segment_a).unwrap();
    let va = volume(&mut k, &mut a);

    let mut k = Kernel::default();
    let mut b = build_union(&mut k, &segment_b).unwrap();
    let vb = volume(&mut k, &mut b);

    let mut k = Kernel::default();
    let mut elbow = build_union(&mut k, &[segment_a[0].clone(), segment_b[0].clone()]).unwrap();
    let check = check_brep(&mut k, &mut elbow, CheckLevel::Level2).unwrap();
    let v_union = volume(&mut k, &mut elbow);

    let mut k = Kernel::default();
    let mut reference =
        build_reference(&mut k, &[segment_a[0].clone(), segment_b[0].clone()]).unwrap();
    let v_reference = volume(&mut k, &mut reference);

    let cylinder = std::f64::consts::PI * R * R * L;
    let smooth = 2.0 * cylinder;
    let overlap = va + vb - v_union;
    let exact_overlap = 4.0 / 3.0 * R * R * R;

    eprintln!("segment A volume  = {va:.12}  (pi r^2 L = {cylinder:.12})");
    eprintln!("segment B volume  = {vb:.12}  (pi r^2 L = {cylinder:.12})");
    eprintln!("union volume      = {v_union:.12}");
    eprintln!("reference volume  = {v_reference:.12}  (2 pi r^2 L = {smooth:.12})");
    eprintln!("overlap measured  = {overlap:.12}  (4/3 r^3 = {exact_overlap:.12})");

    // 1. segment construction / cap orientation: exactly a capped cylinder
    assert!(
        (va - cylinder).abs() < cylinder * 1e-9,
        "segment A not a cylinder"
    );
    assert!(
        (vb - cylinder).abs() < cylinder * 1e-9,
        "segment B not a cylinder"
    );

    // 2. the reference mitre is the smooth pipe (to the sweep/mitre numerics, ~1e-7)
    assert!(
        (v_reference - smooth).abs() < smooth * 1e-5,
        "reference is not the smooth pipe"
    );

    // 3. the union equals the smooth pipe minus the analytic overlap, directly
    assert!(
        (v_union - (smooth - exact_overlap)).abs() < smooth * 1e-6,
        "union is not smooth pipe minus the 4/3 r^3 overlap"
    );

    // 4. the measured overlap (V(A) + V(B) - V(union)) is the closed form 4/3 r^3
    assert!(
        (overlap - exact_overlap).abs() < exact_overlap * 1e-3,
        "overlap {overlap} != 4/3 r^3 {exact_overlap}"
    );

    assert!(check.is_valid(), "{:?}", check.lines());
}
