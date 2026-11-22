---
icon: lucide/package
---

# Release a wheel

<ol class="mesh-flow">
  <li><strong>Verified commit</strong>CI and numerical gates pass.</li>
  <li><strong>Build artifacts</strong>Test the wheels and source distribution.</li>
  <li><strong>PyPI release</strong>Publish those files and verify a fresh install.</li>
</ol>

Use the exact `main` commit you intend to publish. The source archive includes the Rust core. Source builds require Rust 1.95+.

1. **Check the commit.** Require successful CI, Coverage, Performance and Release artifacts runs. Inspect installed wheels for Linux x86-64, Windows x86-64 and macOS ARM64/x86-64, plus the isolated source-distribution build.
2. **Archive the measurements.** Run the [release benchmark suite](performance.md#release-suite-and-public-history) against a compatible baseline. Keep candidate and final-tag results in the [benchmark history](benchmark-history.md).
3. **Publish the tested files.** Dispatch **Publish to PyPI** from `main` with the **Release artifacts** run ID. The workflow verifies the SHA and downloads those files. Configure PyPI Trusted Publishing for `publish-pypi.yml` and the `pypi` GitHub environment.
4. **Verify the public install** in a fresh environment, then tag the same commit and create its GitHub release.

[Publishing workflow](https://github.com/kmarchais/meshers/blob/main/.github/workflows/publish-pypi.yml) · [Wheel requirements](getting-started.md)
