use meshers_core::{
    GyroidOptions, MeshingError, generate_gyroid, generate_gyroid_with_accelerator,
};

#[test]
fn invalid_configuration_is_a_typed_error() {
    for options in [
        GyroidOptions {
            resolution: 0,
            ..GyroidOptions::default()
        },
        GyroidOptions {
            threshold: f64::NAN,
            ..GyroidOptions::default()
        },
        GyroidOptions {
            threads: 257,
            ..GyroidOptions::default()
        },
    ] {
        assert!(matches!(
            generate_gyroid(options),
            Err(MeshingError::InvalidOptions(_))
        ));
    }
}

struct FailingAdapter;
impl meshers_core::accelerator::Factory for FailingAdapter {
    fn create(
        &self,
        _: &meshers_core::Mesh,
        _: &[meshers_core::accelerator::Group],
        _: &[Vec<usize>],
        _: &[usize],
    ) -> Result<Box<dyn meshers_core::accelerator::Optimizer>, String> {
        Err("adapter unavailable".into())
    }
}

#[test]
fn adapter_failure_reaches_the_caller() {
    let result = generate_gyroid_with_accelerator(
        GyroidOptions {
            resolution: 24,
            optimize_passes: 1,
            ..GyroidOptions::default()
        },
        Some(&FailingAdapter),
    );
    assert!(
        matches!(result, Err(MeshingError::GenerationFailed(message)) if message == "adapter unavailable")
    );
}

struct BadDimensions;
impl meshers_core::accelerator::Factory for BadDimensions {
    fn create(
        &self,
        _: &meshers_core::Mesh,
        _: &[meshers_core::accelerator::Group],
        _: &[Vec<usize>],
        _: &[usize],
    ) -> Result<Box<dyn meshers_core::accelerator::Optimizer>, String> {
        Ok(Box::new(BadDimensions))
    }
}
impl meshers_core::accelerator::Optimizer for BadDimensions {
    fn proposals(
        &mut self,
        _: &meshers_core::Mesh,
        _: &[usize],
        _: f64,
        _: f64,
        _: f64,
        _: usize,
    ) -> Result<Vec<meshers_core::Point>, String> {
        Ok(vec![])
    }
    fn synchronize_updates(&mut self, _: &meshers_core::Mesh, _: &[usize]) -> Result<(), String> {
        Ok(())
    }
}

#[test]
fn mismatched_adapter_result_is_rejected() {
    let result = generate_gyroid_with_accelerator(
        GyroidOptions {
            resolution: 24,
            optimize_passes: 1,
            ..GyroidOptions::default()
        },
        Some(&BadDimensions),
    );
    assert!(
        matches!(result, Err(MeshingError::GenerationFailed(message)) if message.contains("candidate dimensions"))
    );
}
