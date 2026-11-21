//! Small, single-threaded browser entry point for the documentation.
use meshers_core::implicit::{self, Options, Region};
use meshers_core::{Mesh, cross, norm, sub};
use serde_json::json;
use std::{cell::RefCell, f64::consts::TAU};

thread_local! {
    static RESPONSE: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

fn minimum_face_angle(mesh: &Mesh) -> f64 {
    // The smallest angle is opposite the shortest side and is always acute.
    let min_sine = mesh
        .surface
        .iter()
        .map(|face| {
            let [a, b, c] = face.map(|i| mesh.points[i]);
            let ab = sub(b, a);
            let ac = sub(c, a);
            let mut sides = [norm(ab), norm(ac), norm(sub(c, b))];
            sides.sort_unstable_by(f64::total_cmp);
            norm(cross(ab, ac)) / (sides[1] * sides[2])
        })
        .fold(1.0_f64, f64::min);
    min_sine.asin().to_degrees()
}

// The worker copies the UTF-8 response before making another call.
#[unsafe(no_mangle)]
pub extern "C" fn generate(field: u32, cells: u32, width: f64) -> *const u8 {
    let result = if field <= 1 && (6..=16).contains(&cells) && (0.3..=0.8).contains(&width) {
        let field = |p: [f64; 3]| {
            let [x, y, z] = p.map(|v| TAU * v);
            if field == 0 {
                x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()
            } else {
                x.cos() + y.cos() + z.cos()
            }
        };
        let options = Options {
            cells: [cells as usize; 3],
            region: Region::Band {
                lower: -width,
                upper: width,
            },
            periodic: [true; 3],
            geometry_tolerance: 0.05,
            optimize_passes: 2,
            threads: 0,
            ..Options::default()
        };
        match implicit::generate(&field, options) {
            Ok(output) => json!({
                "points": output.mesh.points, "faces": output.mesh.surface,
                "tags": output.boundary_tags, "tets": output.mesh.tets.len(),
                "minQuality": output.diagnostics.minimum_mmg_quality,
                "minAngle": minimum_face_angle(&output.mesh)
            }),
            Err(error) => json!({"error": error.to_string()}),
        }
    } else {
        json!({"error": "Field or resolution outside the playground range."})
    };
    RESPONSE.with_borrow_mut(|response| {
        response.clear();
        serde_json::to_writer(&mut *response, &result).expect("writing JSON to a byte buffer");
        response.as_ptr()
    })
}

#[unsafe(no_mangle)]
pub extern "C" fn response_length() -> usize {
    RESPONSE.with_borrow(Vec::len)
}
