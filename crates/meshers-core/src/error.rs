/// Failure to configure or complete generation. No partial mesh is returned.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MeshingError {
    InvalidOptions(String),
    GenerationFailed(String),
}
impl std::fmt::Display for MeshingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidOptions(message) | Self::GenerationFailed(message) => f.write_str(message),
        }
    }
}
impl std::error::Error for MeshingError {}
