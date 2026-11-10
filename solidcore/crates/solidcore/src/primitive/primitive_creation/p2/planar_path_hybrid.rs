//! Production planar-path sweep result: the boolean union of nominal closed segment sweeps.
//!
//! This is a separate engineering result, not the analytical miter research. It is a thin
//! alias of [`PrimitiveCreation::create_sweep_along_planar_path_union`]; the earlier
//! whole-path routing (all straight to the reference, curves to the union) was removed
//! once the union proved robust across straight, curved, mixed, closed and multi-profile
//! paths. The reference [`PrimitiveCreation::create_sweep_along_planar_path`] remains a
//! research implementation and test source; it is not a production entry.
//!
//! The invariants the production union guarantees (input connection, no incorrect solid,
//! temporaries released on a refusal) live in the union module, so there is one
//! implementation of each. The name is kept so the entry point stays stable.
//!
//! Not upstream (a solidcore construction).

use super::PrimitiveCreationAccess;
use crate::curve::curve::Curve;
use crate::error::Result;
use crate::primitive::primitive_creation::p1::PrimitiveCreation;
use crate::topology::BrepId;

impl PrimitiveCreation {
    /// Production planar sweep. Identical to
    /// [`PrimitiveCreation::create_sweep_along_planar_path_union`]; see its documentation
    /// and the union module for the volume-loss and input/result invariants.
    ///
    /// # Errors
    ///
    /// Those of [`PrimitiveCreation::create_sweep_along_planar_path_union`].
    pub fn create_sweep_along_planar_path_hybrid(
        a: &mut impl PrimitiveCreationAccess,
        context_of: BrepId,
        profiles: &[Curve],
        paths: &[Curve],
        tol: f64,
        cap_start: bool,
        cap_end: bool,
        move_profile_center_to_path: bool,
    ) -> Result<BrepId> {
        Self::create_sweep_along_planar_path_union(
            a,
            context_of,
            profiles,
            paths,
            tol,
            cap_start,
            cap_end,
            move_profile_center_to_path,
        )
    }
}
