//! Small, single-threaded browser entry point for the documentation.
use meshers_core::implicit::{self, Options, Region};
use serde_json::json;
use std::{cell::RefCell, f64::consts::TAU};

thread_local! {
    static RESPONSE: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
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
            optimize_passes: 0,
            ..Options::default()
        };
        match implicit::generate(&field, options) {
            Ok(output) => json!({
                "points": output.mesh.points, "faces": output.mesh.surface,
                "tags": output.boundary_tags, "tets": output.mesh.tets.len()
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
