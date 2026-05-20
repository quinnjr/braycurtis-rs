# braycurtis-rs

Rust replacement for the **BrayCurtis** PluMA plugin — pairwise Bray-Curtis
dissimilarity between samples plus a 2-D classical-MDS (PCoA) ordination for
downstream plotting.

```text
    d_jk = Σ_i |x_ij − x_ik| / Σ_i (x_ij + x_ik)
```

Matches `vegan::vegdist(method = "bray")` and phyloseq's
`distance(physeq, method = "bray")`.

## Plugin contract

Implements
[`pluma-plugin-trait`](https://crates.io/crates/pluma-plugin-trait) and exports
the prefixed FFI shims (`BrayCurtis_plugin_create`, `_destroy`, `_input`,
`_run`, `_output`) that PluMA's `Rust::loadRustPlugin` looks up via `dlsym`.

## Parameter file

Same whitespace-delimited shape used by the upstream R plugin (`#` comments
ignored):

```
otufile    CSV/otu_table_normalized.csv
mapping    sample_data.csv            # optional
column     Description                # optional
```

## Output

PluMA passes `output()` an *output prefix*:

- `<prefix>.csv` — square samples×samples dissimilarity matrix.
- `<prefix>.json` — matrix + classical-MDS coordinates + per-sample group
  labels, ready for matplotlib/plotly rendering by `LLMSummarizer` or any
  other consumer.

## Building

```bash
cargo build --release
```

Output: `target/release/libbraycurtis.so`.

## Installing into a PluMA tree

```bash
cd /path/to/PluMA/plugins
mkdir -p BrayCurtis
ln -sfn ../../../braycurtis-rs/target/release/libbraycurtis.so \
        BrayCurtis/libBrayCurtisPlugin.so
ln -sfn ../../../braycurtis-rs/Cargo.toml BrayCurtis/Cargo.toml
```

PluMA's Rust loader globs `<plugins>/<name>/Cargo.toml`, then dlopens the
sibling `lib<name>Plugin.so` and resolves `BrayCurtis_plugin_*` symbols.

## Tests and benchmarks

```bash
cargo test
cargo bench
```

## References

- Bray, J. R. & Curtis, J. T. (1957). *Ecological Monographs* 27(4):325-349.

## License

MIT
