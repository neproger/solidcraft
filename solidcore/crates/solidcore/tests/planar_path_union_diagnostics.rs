//! Independent geometric diagnostics for the union sweep, not relying on
//! `check_brep(Level2)` alone:
//!
//! - self-intersection: the kernel's own face-vs-face intersection of the result with
//!   itself must find no crossing curves;
//! - closedness: the result must be a closed manifold;
//! - boundary coincidence / gaps: an analytic membership oracle for a circular-profile
//!   pipe (the union of the two capped tubes, `dist-to-path <= r` within each segment's
//!   parameter range) is compared against the kernel's point-in-body classification on a
//!   dense sample about the joint. Any gap or overlap in the union boundary shows up as a
//!   disagreement.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]
use solidcore::curve::curve::Curve;
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_curve::p1::{BSplineCurve, NurbCircleParam};
use solidcore::primitive::primitive_creation::model::PrimitiveModel;
use solidcore::primitive::primitive_creation::p1::PrimitiveCreation;
use solidcore::query::api_queries::{edge_closest_point, get_edges};
use solidcore::query::brep::BrepQuery;
use solidcore::topology::brep::p1::brep_intersect_self;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::{EntityRef, Kernel};

const R: f64 = 0.05;
const ARC_RADIUS: f64 = 5.0;

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

fn profile() -> Vec<Curve> {
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

/// The arc (radius 5, first quadrant) then the line from its end along +Y; the joint is
/// at (0, 5, 0).
fn arcline() -> Vec<Curve> {
    vec![
        arc_xy(Vec3::ZERO, ARC_RADIUS, 0.0, 90.0),
        line(
            Vec3::new(0., ARC_RADIUS, 0.),
            Vec3::new(0., ARC_RADIUS + 5.0, 0.),
        ),
    ]
}

fn build(k: &mut Kernel) -> solidcore::Brep {
    let (_, mut b) = PrimitiveModel::build_new(k, |_, model| {
        let id = PrimitiveCreation::create_sweep_along_planar_path_union(
            model,
            PrimitiveModel::MAIN,
            &profile(),
            &arcline(),
            1e-7,
            true,
            true,
            true,
        )?;
        model.take_brep(id)
    })
    .unwrap();
    b.refresh_cached_topology();
    b
}

/// Self-intersections of the body's own faces, excluding the legitimate intersections
/// along shared edges: `brep_intersect_self` returns every face-face intersection,
/// including the shared boundary of adjacent faces. A curve counts as a genuine
/// self-intersection only when a sample of it is not close to any edge of the body.
fn genuine_self_intersections(k: &Kernel, b: &solidcore::Brep, curves: &[Curve]) -> (usize, usize) {
    let edges = get_edges(b).unwrap();
    let edge_tol = 1e-4;
    let mut genuine = 0;
    for curve in curves {
        let domain = curve.natural_interval();
        let mut on_edge = true;
        for j in 0..=4_u32 {
            let f = f64::from(j) / 4.0;
            let t = domain.min() + (domain.max() - domain.min()) * f;
            let Ok(p) = curve.evaluate_point(k, t) else {
                on_edge = false;
                break;
            };
            let mut best = f64::INFINITY;
            for &edge in &edges {
                if let Ok(cp) = edge_closest_point(k, b, edge, p) {
                    best = best.min(cp.closest.distance);
                }
            }
            if best > edge_tol {
                on_edge = false;
                break;
            }
        }
        if !on_edge {
            genuine += 1;
        }
    }
    (curves.len(), genuine)
}

#[test]
fn result_is_a_closed_manifold_without_self_intersections() {
    let mut k = Kernel::default();
    let mut b = build(&mut k);
    assert!(
        is_manifold_solid(&b).unwrap_or(false),
        "not a manifold solid"
    );
    assert!(
        solidcore::query::api_queries::brep_is_closed_manifold(&b).unwrap_or(false),
        "not a closed manifold"
    );
    let (curves, _) = brep_intersect_self(&k, &mut b, 1e-6, 2.0, false).unwrap();
    let (all, genuine) = genuine_self_intersections(&k, &b, &curves);
    eprintln!("self-intersection probes: {all} face-face curves, {genuine} genuine");
    assert_eq!(genuine, 0, "genuine self-intersections: {genuine} of {all}");
}

/// A tiny deterministic RNG so the sample is reproducible.
struct Lcg(u64);
impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Analytic membership of the union of the two capped circular-profile tubes: a point is
/// in the arc tube or the line tube, each `dist <= R` with its own parameter/cap range.
fn analytic_union(p: Vec3) -> (bool, f64) {
    // Arc: radius ARC_RADIUS about the origin in the XY plane, first quadrant.
    let rho = p.x.hypot(p.y);
    let ang = p.y.atan2(p.x);
    let arc_sd = if (0.0..=std::f64::consts::FRAC_PI_2).contains(&ang) {
        ((rho - ARC_RADIUS).powi(2) + p.z * p.z).sqrt() - R
    } else {
        f64::INFINITY
    };
    let arc_in = arc_sd <= 0.0;
    // Line: from (0, ARC_RADIUS, 0) along +Y for 5.
    let t = p.y - ARC_RADIUS;
    let line_sd = if (0.0..=5.0).contains(&t) {
        (p.x * p.x + p.z * p.z).sqrt() - R
    } else {
        f64::INFINITY
    };
    let line_in = line_sd <= 0.0;
    (arc_in || line_in, arc_sd.min(line_sd))
}

#[test]
fn boundary_matches_the_analytic_union_of_the_tubes() {
    let mut k = Kernel::default();
    let mut b = build(&mut k);
    let query = BrepQuery::of_brep(&b);

    let joint = Vec3::new(0., ARC_RADIUS, 0.);
    let mut rng = Lcg(0x5eed_1234_abcd_9876);
    let mut compared = 0usize;
    let mut inside = 0usize;
    let mut outside = 0usize;
    let mut mismatches = 0usize;
    let margin = 0.2 * R;

    let mut attempts = 0usize;
    while compared < 400 && attempts < 20_000 {
        attempts += 1;
        let v = Vec3::new(
            rng.next_unit() * 2.0 - 1.0,
            rng.next_unit() * 2.0 - 1.0,
            rng.next_unit() * 2.0 - 1.0,
        ) * 0.14;
        let p = joint + v;
        let (expected, signed) = analytic_union(p);
        // skip the boundary shell and points with no lateral foot (far from both tubes)
        if !signed.is_finite() || signed.abs() < margin {
            continue;
        }
        let pc = query.point_3d_classify(&k, &mut b, p, 1e-7, false).unwrap();
        let got = match pc.object() {
            Some(EntityRef::Region(region)) => {
                b.region(region).map(|r| r.is_void() == 0).unwrap_or(false)
            }
            _ => false,
        };
        compared += 1;
        if got {
            inside += 1;
        } else {
            outside += 1;
        }
        if got != expected {
            mismatches += 1;
            if mismatches <= 5 {
                eprintln!(
                    "  mismatch at ({:.5},{:.5},{:.5}): analytic={expected} kernel={got} sd={signed:.4e}",
                    p.x, p.y, p.z
                );
            }
        }
    }

    eprintln!(
        "membership: compared={compared} inside={inside} outside={outside} mismatches={mismatches}"
    );
    assert!(compared >= 200, "too few points compared: {compared}");
    assert!(
        inside > 20 && outside > 20,
        "sample not straddling the boundary"
    );
    assert_eq!(mismatches, 0, "boundary disagrees with the analytic union");
}
