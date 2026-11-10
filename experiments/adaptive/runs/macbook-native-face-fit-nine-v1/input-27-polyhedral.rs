//! Exact convex affine forms and bounded min/max combinations of those forms.
//! Unsupported graphs keep their original scalar meshing path.
use super::Node;
use meshers_core::implicit::intersection::MAX_CONSTRAINTS;

type Plane = [f64; 4];

fn constant(forms: &[Plane]) -> Option<f64> {
    (forms.len() == 1 && forms[0][..3] == [0.; 3]).then(|| forms[0][3])
}

fn scaled(forms: &[Plane], factor: f64) -> Option<Vec<Plane>> {
    // Negating a maximum makes a minimum, except for a single affine form.
    if factor < 0. && forms.len() != 1 {
        return None;
    }
    Some(forms.iter().map(|p| p.map(|x| x * factor)).collect())
}

fn affine_forms(nodes: &[Node]) -> Vec<Option<Vec<Plane>>> {
    let mut forms: Vec<Option<Vec<Plane>>> = Vec::with_capacity(nodes.len());
    for (op, args, value) in nodes {
        let result = (|| {
            let operand = |i: usize| forms[args[i]].as_deref();
            match op.as_str() {
                "constant" => Some(vec![[0., 0., 0., *value]]),
                "input" => {
                    let mut p = [0.; 4];
                    p[args[0]] = 1.;
                    Some(vec![p])
                }
                "add" => {
                    let (a, b) = (operand(0)?, operand(1)?);
                    if a.len() * b.len() > MAX_CONSTRAINTS {
                        return None;
                    }
                    Some(
                        a.iter()
                            .flat_map(|a| {
                                b.iter().map(move |b| std::array::from_fn(|i| a[i] + b[i]))
                            })
                            .collect(),
                    )
                }
                "sub" => {
                    let (a, b) = (operand(0)?, operand(1)?);
                    if b.len() != 1 {
                        return None;
                    }
                    Some(
                        a.iter()
                            .map(|a| std::array::from_fn(|i| a[i] - b[0][i]))
                            .collect(),
                    )
                }
                "mul" => {
                    let (a, b) = (operand(0)?, operand(1)?);
                    if let Some(factor) = constant(a) {
                        scaled(b, factor)
                    } else {
                        scaled(a, constant(b)?)
                    }
                }
                "div" => scaled(operand(0)?, 1. / constant(operand(1)?)?),
                "neg" => scaled(operand(0)?, -1.),
                "abs" => {
                    let a = operand(0)?;
                    if a.len() != 1 {
                        return None;
                    }
                    Some(vec![a[0], a[0].map(|x| -x)])
                }
                "maximum" => Some(operand(0)?.iter().chain(operand(1)?).copied().collect()),
                _ => None,
            }
        })();
        let result = result.and_then(|result| {
            let mut unique = Vec::new();
            for p in result {
                if p.iter().any(|v| !v.is_finite()) {
                    return None;
                }
                if !unique.contains(&p) {
                    unique.push(p);
                }
            }
            (unique.len() <= MAX_CONSTRAINTS).then_some(unique)
        });
        forms.push(result);
    }
    forms
}

pub(super) fn halfspaces(nodes: &[Node], output: usize) -> Option<Vec<Plane>> {
    let result = affine_forms(nodes)[output].take()?;
    // Zero-gradient branches need separate full/empty-domain semantics.
    // Leave those and single-plane fields on the existing scalar path.
    (result.len() >= 2 && result.iter().all(|p| p[..3] != [0.; 3])).then_some(result)
}

pub(super) type Regions = (Vec<Plane>, Vec<u32>);

/// Bounded min/max combinations of convex affine forms. Arithmetic outside a
/// nonconvex combination is deliberately unsupported rather than approximated.
pub(super) fn regions(nodes: &[Node], output: usize) -> Option<Regions> {
    let convex = affine_forms(nodes);
    let mut groups: Vec<Option<Vec<Vec<Plane>>>> = Vec::with_capacity(nodes.len());
    for (i, (op, args, _)) in nodes.iter().enumerate() {
        let result = if let Some(forms) = &convex[i] {
            Some(vec![forms.clone()])
        } else {
            (|| {
                let (a, b) = (
                    groups[*args.first()?].as_ref()?,
                    groups[*args.get(1)?].as_ref()?,
                );
                match op.as_str() {
                    "minimum" => Some(a.iter().chain(b).cloned().collect()),
                    "maximum" if a.len() * b.len() <= 32 => Some(
                        a.iter()
                            .flat_map(|a| {
                                b.iter().map(move |b| {
                                    let mut group = a.clone();
                                    for plane in b {
                                        if !group.contains(plane) {
                                            group.push(*plane);
                                        }
                                    }
                                    group
                                })
                            })
                            .collect(),
                    ),
                    _ => None,
                }
            })()
        };
        groups.push(result.filter(|groups| {
            groups.len() <= 32 && groups.iter().all(|group| group.len() <= MAX_CONSTRAINTS)
        }));
    }
    let groups = groups[output].take()?;
    let mut planes = Vec::new();
    let mut masks = Vec::new();
    for group in groups {
        let mut mask = 0;
        for plane in group {
            if plane[..3] == [0.; 3] {
                return None;
            }
            let index = if let Some(index) = planes.iter().position(|p| *p == plane) {
                index
            } else {
                if planes.len() == MAX_CONSTRAINTS {
                    return None;
                }
                planes.push(plane);
                planes.len() - 1
            };
            mask |= 1 << index;
        }
        if !masks.contains(&mask) {
            masks.push(mask);
        }
    }
    (planes.len() >= 2 && masks.len() >= 2).then_some((planes, masks))
}
