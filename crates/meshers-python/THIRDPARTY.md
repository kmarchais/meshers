# Dependency license notices

`THIRDPARTY.json` contains the Rust dependency notices shipped in Python wheels
and source distributions. It includes normal dependencies across platforms from
the Python extension's Cargo lockfile. Build-only dependencies are excluded.

Regenerate from this directory with cargo-bundle-licenses 4.2.0:

```sh
cargo bundle-licenses --format json --output THIRDPARTY.json
```

Review warnings and all `NOT FOUND` entries before committing. Four Wasmtime
crates omit the upstream license file from their published source archives:
cranelift-assembler-x64, cranelift-bitset, wasmtime-internal-core and
wasmtime-internal-jit-icache-coherence. Their Apache-2.0 WITH LLVM-exception
notice was copied from cranelift-codegen 0.135.2, from the same Wasmtime source
revision `e9f1ea232fd245aea338ab3eb7d73487ae75cab1`.

After changing dependencies, regenerate and review the notices, then inspect
the built wheel and source archive for `THIRDPARTY.json`.
