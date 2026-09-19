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

1. Merge release preparation and require all four validation workflows to pass
   on the exact main commit.
2. Run **Publish to PyPI** from main, supplying the successful **Release artifacts**
   run ID. The workflow verifies the source SHA and all check results, downloads
   the tested artifacts, and publishes through the `pypi` environment.
3. Verify installation from PyPI in a fresh environment, then tag the same source
   commit and create its GitHub release.

The publishing workflow is `.github/workflows/publish-pypi.yml`. Configure that
filename and the repository in PyPI Trusted Publishing. The GitHub environment
is `pypi`. Publishing is explicitly triggered and never runs on pull requests.
It reuses tested files rather than rebuilding them during upload.

The source archive includes the Rust core and does not require a prior crates.io
release. A conda-forge package is not part of v0.1.
