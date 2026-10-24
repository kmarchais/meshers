//! Exact max-of-affine decomposition, bounded by the constrained cutter's capacity.
//! Unsupported or nonconvex graphs keep their original scalar meshing path.
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

pub(super) fn halfspaces(nodes: &[Node], output: usize) -> Option<Vec<Plane>> {
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
    let result = forms[output].take()?;
    // Zero-gradient branches need separate full/empty-domain semantics.
    // Leave those and single-plane fields on the existing scalar path.
    (result.len() >= 2 && result.iter().all(|p| p[..3] != [0.; 3])).then_some(result)
}
