//! Production entry check: `create_sweep_along_planar_path_hybrid` is now the union of
//! nominal closed segment sweeps (the pilot's single production algorithm). It must be
//! identical to `create_sweep_along_planar_path_union` on every path shape — linear,
//! curved, closed and mixed — and return a valid manifold solid where one is expected.
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

#[derive(Clone, Copy, PartialEq, Eq)]
enum Entry {
    Union,
    Hybrid,
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
    entry: Entry,
    profiles: &[Curve],
    paths: &[Curve],
) -> Result<solidcore::Brep, SolidError> {
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = match entry {
            Entry::Union => PrimitiveCreation::create_sweep_along_planar_path_union(
                model,
                PrimitiveModel::MAIN,
                profiles,
                paths,
                1e-7,
                true,
                true,
                true,
            )?,
            Entry::Hybrid => PrimitiveCreation::create_sweep_along_planar_path_hybrid(
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
        1e-9,
        &PropertyAttributes {
            store: &store,
            key_of: &|_| ObjectKey(0),
        },
    )
    .unwrap()
}

fn same_as_union(name: &str, profiles: &[Curve], paths: &[Curve], check_level: bool) {
    let mut k = Kernel::default();
    let mut union = build(&mut k, Entry::Union, profiles, paths).unwrap();
    let union_volume = volume(&mut k, &mut union);
    let union_faces = union.face_list().len();

    let mut k = Kernel::default();
    let mut hybrid = build(&mut k, Entry::Hybrid, profiles, paths).unwrap();
    let manifold = is_manifold_solid(&hybrid).unwrap_or(false);
    let hybrid_volume = volume(&mut k, &mut hybrid);
    let hybrid_faces = hybrid.face_list().len();
    eprintln!(
        "{name:24} hybrid==union volume={hybrid_volume:.9} faces={hybrid_faces} manifold={manifold}"
    );

    assert!(
        (hybrid_volume - union_volume).abs() <= 1e-9 * union_volume.abs().max(1e-300),
        "{name}: hybrid volume differs from union"
    );
    assert_eq!(hybrid_faces, union_faces, "{name}: face count differs");
    assert!(manifold, "{name}: hybrid is not a manifold solid");
    if check_level {
        let check = check_brep(&mut k, &mut hybrid, CheckLevel::Level2).unwrap();
        assert!(check.is_valid(), "{name}: invalid: {:?}", check.lines());
    }
}

#[test]
fn production_entry_is_union_on_all_path_shapes() {
    // linear elbow
    let elbow = [
        line(Vec3::ZERO, Vec3::new(0., 0., 1.)),
        line(Vec3::new(0., 0., 1.), Vec3::new(1., 0., 1.)),
    ];
    same_as_union("linear elbow", &circle(0.1), &elbow, false);

    // arc + line (curved joint)
    let arcline = [
        arc_xy(Vec3::ZERO, 5.0, 0.0, 90.0),
        line(Vec3::new(0., 5., 0.), Vec3::new(0., 10., 0.)),
    ];
    same_as_union("arc+line", &circle(0.05), &arcline, true);

    // closed square ring
    let sq = 2.0;
    let ring = [
        line(Vec3::ZERO, Vec3::new(sq, 0., 0.)),
        line(Vec3::new(sq, 0., 0.), Vec3::new(sq, sq, 0.)),
        line(Vec3::new(sq, sq, 0.), Vec3::new(0., sq, 0.)),
        line(Vec3::new(0., sq, 0.), Vec3::ZERO),
    ];
    same_as_union("closed square ring", &circle(0.05), &ring, false);

    // mixed line + line + arc
    let mixed = [
        line(Vec3::ZERO, Vec3::new(2., 0., 0.)),
        line(Vec3::new(2., 0., 0.), Vec3::new(2., 2., 0.)),
        arc_xy(Vec3::new(2., 4., 0.), 2.0, -90.0, 0.0),
    ];
    same_as_union("line+line+arc", &circle(0.05), &mixed, false);
}
