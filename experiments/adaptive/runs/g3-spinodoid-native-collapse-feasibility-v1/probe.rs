use super::*;
use crate::geometry::Geometry;

#[derive(serde::Deserialize)]
struct Input {
    points: Vec<Point>,
    tetrahedra: Vec<[usize; 4]>,
    surface: Vec<Face>,
    coefficients: Vec<[f64; 5]>,
    threshold: f64,
    origins: Vec<Point>,
}

impl Geometry for Input {
    fn value(&self, point: Point) -> f64 {
        self.threshold
            - self
                .coefficients
                .iter()
                .map(|mode| {
                    let angle = std::f64::consts::TAU
                        * (0..3)
                            .map(|axis| mode[axis] * (point[axis] + 0.5))
                            .sum::<f64>()
                        + mode[4];
                    mode[3] * angle.cos()
                })
                .sum::<f64>()
    }
    fn gradient(&self, point: Point) -> Point {
        let mut gradient = [0.; 3];
        for mode in &self.coefficients {
            let angle = std::f64::consts::TAU
                * (0..3)
                    .map(|axis| mode[axis] * (point[axis] + 0.5))
                    .sum::<f64>()
                + mode[4];
            for axis in 0..3 {
                gradient[axis] += std::f64::consts::TAU * mode[3] * mode[axis] * angle.sin();
            }
        }
        gradient
    }
    fn rank(&self, point: Point) -> [i64; 3] {
        point.map(|value| ((value + 0.5).rem_euclid(1.) * 1e11).round() as i64)
    }
    fn bounds(&self) -> [Point; 2] {
        [[-0.5; 3], [0.5; 3]]
    }
    fn band(&self) -> bool {
        false
    }
    fn minimum_quality(&self) -> f64 {
        0.1
    }
}

#[test]
#[ignore]
fn inspect_saved_pending_collapse() {
    let input: Input = serde_json::from_slice(
        &std::fs::read(std::env::var("MESHERS_COLLAPSE_INPUT").unwrap()).unwrap(),
    )
    .unwrap();
    let mesh = Mesh {
        points: input.points.clone(),
        tets: input.tetrahedra.clone(),
        surface: input.surface.clone(),
    };
    let mut rows = Vec::new();
    for origin in &input.origins {
        for guard in [0.3, 0.1] {
            let options = Options {
                tolerance: 0.002,
                min_length: 0.002,
                max_length: 0.1,
                quality_guard: guard,
                ..Default::default()
            };
            let mut state = Tets::new(&mesh, &input, 0., [true; 3]);
            let canonical = state.wrap(*origin);
            let vertex = state
                .points
                .iter()
                .position(|&point| norm(state.delta(canonical, point)) < 1e-11)
                .unwrap();
            let neighbours = state.adjacent(vertex);
            let size = state.sizes(&options, &state.edges()).unwrap();
            let mut accepted = Vec::new();
            let mut eligible = 0;
            let mut mask_compatible = 0;
            for target in neighbours {
                if input.value(state.points[target]).abs() > 1e-10
                    || state.masks[target] & WALL == 0
                {
                    continue;
                }
                eligible += 1;
                if state.masks[vertex] & !state.masks[target] == 0 {
                    mask_compatible += 1;
                }
                let target_point = state.points[target];
                let length = state.length(vertex, target);
                if state.collapse(vertex, target, &size, &options) {
                    let minimum = state
                        .tets
                        .iter()
                        .enumerate()
                        .filter(|(index, _)| state.alive[*index])
                        .map(|(_, &tet)| mmg(state.tet_points(tet, None)))
                        .fold(1., f64::min);
                    accepted.push(serde_json::json!({"target": target_point, "distance": length, "minimum_signed_quality_after": minimum, "live_tetrahedra": state.live_tetrahedra}));
                    state = Tets::new(&mesh, &input, 0., [true; 3]);
                }
            }
            rows.push(serde_json::json!({"origin": origin, "quality_guard": guard, "vertex_mask": state.masks[vertex], "vertex_frozen": state.frozen[vertex], "root_neighbours": eligible, "mask_compatible_neighbours": mask_compatible, "native_accepted_collapses": accepted}));
        }
    }
    let output = serde_json::json!({"rows": rows, "scope": "Native collapse predicate on reconstructed torus and exact spectral field; no live pending history or repaired-mesh admission inferred"});
    std::fs::write(
        std::env::var("MESHERS_COLLAPSE_OUTPUT").unwrap(),
        serde_json::to_vec_pretty(&output).unwrap(),
    )
    .unwrap();
}
