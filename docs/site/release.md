# Releasing version 0.1

Python `meshers` and Rust `meshers-core` use version 0.1.0. Python wheels support
CPython 3.10+; source builds require Rust 1.95 or newer. The first publication
target is PyPI. The Rust core can be published independently to crates.io.

## Required checks

Run CI on the exact source commit intended for release. A push to main runs
these workflows without requiring a pull request:

- Rust tests, formatting, Clippy and crate package verification.
- Python tests and executable examples on Linux, Windows and macOS.
- Python lint and types, solver validation and documentation construction.
- Coverage reports and repeated runtime, memory and quality measurements.
- Installed wheels for Linux x86-64, Windows x86-64, macOS ARM64 and macOS x86-64.
- An isolated source-distribution build and installation.

Inspect the workflow results and artifacts. A configured workflow is not evidence
that its checks passed. Local Linux wheels may require a newer glibc than the
portable Linux wheel built by the release workflow.

## Publish tested artifacts

1. Require all release checks to pass on the final commit.
2. Download that commit's tested wheels and source distribution. Verify package
   versions, platform tags, license notices and hashes.
3. Verify PyPI ownership and publishing authorization, then publish those files.
4. Tag the reviewed source and attach release notes describing supported behavior.

The workflows build and test packages; they do not upload to PyPI automatically.
The candidate has not yet been published. The source archive includes the Rust
core and does not require a prior crates.io release. A conda-forge package is
not part of v0.1.
