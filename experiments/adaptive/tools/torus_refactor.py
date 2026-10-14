"""One-off refactor of tet_adapt.rs for periodic torus adaptation."""

import pathlib

p = pathlib.Path("crates/meshers-core/src/tet_adapt.rs")
s = p.read_text()


def rep(old: str, new: str) -> None:
    global s
    assert old in s, old[:100]
    s = s.replace(old, new, 1)


rep(
    "//! Faces on periodic box sides are frozen, so an initially periodic mesh stays periodic.",
    """//! Periodic axes are adapted on the torus: copies of a vertex on opposite box faces are
//! merged, geometry uses minimal-image displacements, and the periodic box planes become
//! internal constrained faces. Unfolding the result recreates identical opposite faces.""",
)
rep(
    """    vb: Vec<Vec<Face>>,
    periodic: [bool; 3],
}""",
    """    vb: Vec<Vec<Face>>,
    /// Box length along periodic axes.
    period: [Option<f64>; 3],
}""",
)
start = s.index("    fn new(mesh: &Mesh, g: &'g G, threshold: f64, periodic: [bool; 3]) -> Self {")
end = s.index("    fn add_face(&mut self, f: Face, label: u8) {")
s = (
    s[:start]
    + """    fn new(mesh: &Mesh, g: &'g G, threshold: f64, periodic: [bool; 3]) -> Self {
        let [lo, hi] = g.bounds();
        let period = std::array::from_fn(|a| periodic[a].then(|| hi[a] - lo[a]));
        // Merge periodic copies; the representative lies on the lower plane.
        let mut index = HashMap::new();
        let mut points = Vec::new();
        let map: Vec<usize> = mesh
            .points
            .iter()
            .map(|&p| {
                let c: Point = std::array::from_fn(|a| {
                    if periodic[a] && p[a] == hi[a] { lo[a] } else { p[a] }
                });
                let key: [i64; 3] = std::array::from_fn(|a| {
                    ((c[a] - lo[a]) / (hi[a] - lo[a]) * 1e11).round() as i64
                });
                *index.entry(key).or_insert_with(|| {
                    points.push(c);
                    points.len() - 1
                })
            })
            .collect();
        let n = points.len();
        let tets: Vec<[usize; 4]> = mesh.tets.iter().map(|t| t.map(|i| map[i])).collect();
        let mut s = Self {
            g,
            threshold,
            points,
            masks: vec![0; n],
            normals: vec![[0.; 3]; n],
            frozen: vec![false; n],
            alive: vec![true; tets.len()],
            tets,
            vt: vec![Vec::new(); n],
            faces: HashMap::new(),
            vb: vec![Vec::new(); n],
            period,
        };
        for (t, tet) in s.tets.iter().enumerate() {
            for &v in tet {
                s.vt[v].push(t);
            }
        }
        for &f in &mesh.surface {
            let p = f.map(|i| mesh.points[i]);
            let mut label = 0;
            for a in 0..3 {
                if p.iter().all(|q| q[a] == lo[a]) {
                    label = 2 + 2 * a as u8;
                } else if p.iter().all(|q| q[a] == hi[a]) {
                    // Opposite periodic faces are one internal plane on the torus.
                    label = if periodic[a] { 2 } else { 3 } + 2 * a as u8;
                }
            }
            let f = f.map(|i| map[i]);
            if !s.faces.contains_key(&sorted(f)) {
                s.add_face(f, label);
            }
        }
        for v in 0..n {
            s.refresh(v);
        }
        s
    }

    /// Minimal-image displacement from `a` to `b`.
    fn delta(&self, a: Point, b: Point) -> Point {
        std::array::from_fn(|k| {
            let d = b[k] - a[k];
            self.period[k].map_or(d, |l| d - l * (d / l).round())
        })
    }

    /// Bring a point into the box along periodic axes.
    fn wrap(&self, mut p: Point) -> Point {
        let [lo, _] = self.g.bounds();
        for a in 0..3 {
            if let Some(l) = self.period[a] {
                p[a] = lo[a] + (p[a] - lo[a]).rem_euclid(l);
            }
        }
        p
    }

    /// Unwrapped coordinates of a simplex, relative to its first vertex.
    fn local<const N: usize>(&self, raw: [Point; N]) -> [Point; N] {
        let base = raw[0];
        raw.map(|p| {
            let d = self.delta(base, p);
            std::array::from_fn(|k| base[k] + d[k])
        })
    }

"""
    + s[end:]
)
rep(
    """        self.frozen[v] = (0..3).any(|a| self.periodic[a] && (mask >> (2 + 2 * a)) & 3 != 0)
            || mask.count_ones() >= 3;""",
    """        self.frozen[v] = mask.count_ones() >= 3;""",
)
rep(
    """            if r.abs() / g2.sqrt() <= scale * 1e-13 {
                return self.g.contains(p).then_some(p);
            }""",
    """            if r.abs() / g2.sqrt() <= scale * 1e-13 {
                let p = self.wrap(p);
                return self.g.contains(p).then_some(p);
            }""",
)
rep(
    """        if mask & WALL == 0 {
            return Some(p);
        }""",
    """        if mask & WALL == 0 {
            return Some(self.wrap(p));
        }""",
)
rep(
    """    fn tet_points(&self, t: [usize; 4], moved: Option<(usize, Point)>) -> [Point; 4] {
        t.map(|i| match moved {
            Some((j, q)) if i == j => q,
            _ => self.points[i],
        })
    }""",
    """    fn tet_points(&self, t: [usize; 4], moved: Option<(usize, Point)>) -> [Point; 4] {
        self.local(t.map(|i| match moved {
            Some((j, q)) if i == j => q,
            _ => self.points[i],
        }))
    }""",
)
rep(
    """        norm(sub(self.points[u], self.points[v]))""",
    """        norm(self.delta(self.points[u], self.points[v]))""",
)
rep(
    """        norm(sub(pu, pv)) * angle(nu, nv) / 6.""",
    """        norm(self.delta(pu, pv)) * angle(nu, nv) / 6.""",
)
rep(
    """        let mask = self.edge_mask(&bfaces);
        if (0..3).any(|a| self.periodic[a] && (mask >> (2 + 2 * a)) & 3 != 0) {
            return None;
        }
        let mid = std::array::from_fn(|a| 0.5 * (self.points[u][a] + self.points[v][a]));""",
    """        let mask = self.edge_mask(&bfaces);
        let d = self.delta(self.points[u], self.points[v]);
        let mid = self.wrap(std::array::from_fn(|a| self.points[u][a] + 0.5 * d[a]));""",
)
rep(
    """                let q = c.map(|i| if i == usize::MAX { p } else { self.points[i] });""",
    """                let q = self.local(c.map(|i| if i == usize::MAX { p } else { self.points[i] }));""",
)
rep(
    """        if mu != 0 && (bfaces.len() != 2 || self.edge_mask(&bfaces) != mu) {""",
    """        if mu != 0 && self.edge_mask(&bfaces) != mu {""",
)
rep(
    """            let [pa, pb, pc] = [a, b, c].map(|i| self.points[i]);
            let du = determinant([pa, pb, pc, self.points[u]]);
            let dv = determinant([pa, pb, pc, self.points[v]]);""",
    """            let du = determinant(self.tet_points([a, b, c, u], None));
            let dv = determinant(self.tet_points([a, b, c, v], None));""",
)
rep(
    """                    .map(|t| mmg(t.map(|w| self.points[w])))""",
    """                    .map(|&t| mmg(self.tet_points(t, None)))""",
)
rep(
    """            ts.map(|t| determinant(t.map(|w| self.points[w])).abs()).sum()""",
    """            ts.map(|t| determinant(self.tet_points(t, None)).abs()).sum()""",
)
rep(
    """            if determinant(tet.map(|k| self.points[k])) < 0. {""",
    """            if determinant(self.tet_points(tet, None)) < 0. {""",
)
rep(
    """            .map(|&t| mmg(t.map(|k| self.points[k])))""",
    """            .map(|&t| mmg(self.tet_points(t, None)))""",
)
rep(
    """                .map(|&t| determinant(t.map(|k| self.points[k])).abs())""",
    """                .map(|&t| determinant(self.tet_points(t, None)).abs())""",
)
rep(
    """        let centroid: Point = std::array::from_fn(|a| {
            adj.iter().map(|&w| self.points[w][a]).sum::<f64>() / adj.len() as f64
        });""",
    """        let mean: Point = adj.iter().fold([0.; 3], |m, &w| {
            let d = self.delta(p, self.points[w]);
            std::array::from_fn(|a| m[a] + d[a] / adj.len() as f64)
        });
        let centroid: Point = std::array::from_fn(|a| p[a] + mean[a]);""",
)
rep(
    """                let p = f.map(|i| self.points[i]);
                let c = std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.);""",
    """                let p = self.local(f.map(|i| self.points[i]));
                let c = std::array::from_fn(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.);""",
)
start = s.index("    fn into_mesh(self) -> Result<Mesh, String> {")
end = s.index("/// Adapt a valid boundary-conforming mesh")
s = (
    s[:start]
    + """    /// Unfold the torus into the box: each tetrahedron moves to the cell holding its
    /// centroid, and plane vertices get one copy per side they are used on.
    fn into_mesh(self) -> Result<Mesh, String> {
        let [lo, hi] = self.g.bounds();
        let mut map = HashMap::<(usize, [u8; 3]), usize>::new();
        let mut points = Vec::new();
        let mut tets = Vec::new();
        for (t, tet) in self.tets.iter().enumerate() {
            if !self.alive[t] {
                continue;
            }
            let mut p = self.tet_points(*tet, None);
            for a in 0..3 {
                if let Some(l) = self.period[a] {
                    let c = p.iter().map(|q| q[a]).sum::<f64>() / 4.;
                    let shift = ((c - lo[a]) / l).floor() * l;
                    for q in &mut p {
                        q[a] -= shift;
                    }
                }
            }
            let mut out = [0; 4];
            for (k, &v) in tet.iter().enumerate() {
                let mut image = [0u8; 3];
                let mut q = self.points[v];
                for a in 0..3 {
                    if let Some(l) = self.period[a] {
                        let m = ((p[k][a] - q[a]) / l).round();
                        if m == 1. && q[a] == lo[a] {
                            image[a] = 1;
                            q[a] = hi[a];
                        } else if m != 0. {
                            return Err("a tetrahedron crosses a periodic plane".into());
                        }
                    }
                }
                out[k] = *map.entry((v, image)).or_insert_with(|| {
                    points.push(q);
                    points.len() - 1
                });
            }
            tets.push(out);
        }
        let surface = crate::boundary(&tets)?;
        Ok(Mesh {
            points,
            tets,
            surface,
        })
    }
}

"""
    + s[end:]
)
rep(
    """    let mut s = Tets::new(mesh, g, threshold, periodic);""",
    """    let [lo, hi] = g.bounds();
    if (0..3).any(|a| periodic[a] && o.max_length > 0.25 * (hi[a] - lo[a])) {
        return Err("periodic adaptation needs max_length <= a quarter of the period".into());
    }
    let mut s = Tets::new(mesh, g, threshold, periodic);""",
)
p.write_text(s)
print("ok")
