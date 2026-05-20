//! Bray-Curtis dissimilarity plugin for PluMA.
//!
//! ```text
//!     d_jk = Σ_i |x_ij − x_ik| / Σ_i (x_ij + x_ik)
//! ```
//!
//! Matches `vegan::vegdist(method = "bray")` and phyloseq's
//! `distance(physeq, method = "bray")`. The plugin reads the same
//! PluMA-style parameter file shape used by the upstream R plugin
//! (`otufile`, `mapping`, `column`), writes a square samples×samples
//! dissimilarity CSV, and a sidecar JSON containing both the matrix and a
//! classical-MDS 2-D ordination suitable for matplotlib / plotly rendering.

use pluma_plugin_trait::PluMAPlugin;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;
use std::ffi::CStr;
use std::os::raw::c_char;
use std::path::{Path, PathBuf};

/// Compute Bray-Curtis dissimilarity between two non-negative abundance
/// vectors.
#[inline]
pub fn bray_curtis(a: &[f64], b: &[f64]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let mut num = 0.0;
    let mut den = 0.0;
    for (&x, &y) in a.iter().zip(b.iter()) {
        num += (x - y).abs();
        den += x + y;
    }
    if den <= 0.0 {
        0.0
    } else {
        (num / den).clamp(0.0, 1.0)
    }
}

// ---------------------------------------------------------------------------
// Parameter / CSV I/O
// ---------------------------------------------------------------------------

/// Resolve a path referenced *inside* a parameter file against the PluMA
/// pipeline convention (`<prefix>/parameters/<name>.txt` → resolve against
/// `<prefix>`). Falls back to the parameter file's parent directory, then
/// the raw value.
pub fn resolve_input_path(param_file: &Path, value: &str) -> PathBuf {
    let p = Path::new(value);
    if p.is_absolute() {
        return p.to_path_buf();
    }
    let param_dir = param_file.parent();
    let prefix = param_dir.and_then(Path::parent);
    if let Some(pre) = prefix {
        let candidate = pre.join(p);
        if candidate.exists() {
            return candidate;
        }
    }
    if let Some(pd) = param_dir {
        let candidate = pd.join(p);
        if candidate.exists() {
            return candidate;
        }
    }
    prefix
        .map(|pre| pre.join(p))
        .unwrap_or_else(|| p.to_path_buf())
}

fn read_parameters<P: AsRef<Path>>(path: P) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let text = std::fs::read_to_string(path)?;
    let mut params = HashMap::new();
    for raw in text.lines() {
        let line = raw.split('#').next().map(|s| s.trim()).unwrap_or("");
        if line.is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        if let (Some(k), Some(v)) = (it.next(), it.next()) {
            params.insert(k.to_string(), v.to_string());
        }
    }
    Ok(params)
}

fn read_otu_csv<P: AsRef<Path>>(
    path: P,
) -> Result<(Vec<String>, Vec<String>, Vec<Vec<f64>>), Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;
    let headers = reader.headers()?.clone();
    let samples: Vec<String> = headers.iter().skip(1).map(|s| s.to_string()).collect();
    let num_samples = samples.len();
    let mut otus = Vec::new();
    let mut counts: Vec<Vec<f64>> = vec![Vec::new(); num_samples];
    for record in reader.records() {
        let record = record?;
        let mut it = record.iter();
        if let Some(otu) = it.next() {
            otus.push(otu.to_string());
            for (i, field) in it.enumerate() {
                if i < num_samples {
                    counts[i].push(field.trim().parse::<f64>().unwrap_or(0.0));
                }
            }
        }
    }
    Ok((samples, otus, counts))
}

fn read_metadata_column<P: AsRef<Path>>(
    path: P,
    column: &str,
) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_path(path)?;
    let headers = reader.headers()?.clone();
    let col_idx = match headers.iter().position(|h| h == column) {
        Some(i) => i,
        None => return Ok(HashMap::new()),
    };
    let mut out = HashMap::new();
    for record in reader.records() {
        let record = record?;
        if let Some(s) = record.get(0) {
            if let Some(v) = record.get(col_idx) {
                out.insert(s.to_string(), v.to_string());
            }
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// Plugin
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize)]
pub struct BrayCurtisPlotData {
    pub samples: Vec<String>,
    pub matrix: Vec<Vec<f64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mds: Option<Vec<[f64; 2]>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mds_eigenvalues: Option<[f64; 2]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub groups: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub group_column: Option<String>,
}

#[derive(Debug, Clone)]
pub struct BrayCurtisPlugin {
    samples: Vec<String>,
    otus: Vec<String>,
    counts: Vec<Vec<f64>>, // counts[sample][otu]
    dissimilarity: Vec<Vec<f64>>,
    mapping_path: Option<String>,
    group_column: Option<String>,
}

impl Default for BrayCurtisPlugin {
    fn default() -> Self {
        Self {
            samples: Vec::new(),
            otus: Vec::new(),
            counts: Vec::new(),
            dissimilarity: Vec::new(),
            mapping_path: None,
            group_column: None,
        }
    }
}

impl BrayCurtisPlugin {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_matrix(samples: Vec<String>, otus: Vec<String>, counts: Vec<Vec<f64>>) -> Self {
        Self {
            samples,
            otus,
            counts,
            ..Default::default()
        }
    }

    pub fn dissimilarity(&self) -> &Vec<Vec<f64>> {
        &self.dissimilarity
    }

    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    pub fn num_samples(&self) -> usize {
        self.samples.len()
    }

    pub fn num_otus(&self) -> usize {
        self.otus.len()
    }

    fn resolve_group_labels(&self) -> Option<Vec<String>> {
        let path = self.mapping_path.as_ref()?;
        let col = self.group_column.as_ref()?;
        let map = read_metadata_column(path, col).ok()?;
        if map.is_empty() {
            return None;
        }
        Some(
            self.samples
                .iter()
                .map(|s| map.get(s).cloned().unwrap_or_default())
                .collect(),
        )
    }

    fn write_csv<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn Error>> {
        let mut writer = csv::Writer::from_path(path)?;
        let mut header = vec!["".to_string()];
        header.extend(self.samples.iter().cloned());
        writer.write_record(&header)?;
        for (i, s) in self.samples.iter().enumerate() {
            let mut row = vec![s.clone()];
            for j in 0..self.samples.len() {
                row.push(format!("{:.6}", self.dissimilarity[i][j]));
            }
            writer.write_record(&row)?;
        }
        writer.flush()?;
        Ok(())
    }

    fn write_json<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn Error>> {
        let groups = self.resolve_group_labels();
        let (mds, eigs) = if self.samples.len() >= 3 {
            classical_mds_2d(&self.dissimilarity)
        } else {
            (None, None)
        };
        let payload = BrayCurtisPlotData {
            samples: self.samples.clone(),
            matrix: self.dissimilarity.clone(),
            mds,
            mds_eigenvalues: eigs,
            groups,
            group_column: self.group_column.clone(),
        };
        let json = serde_json::to_string_pretty(&payload)?;
        std::fs::write(path, json)?;
        Ok(())
    }
}

impl PluMAPlugin for BrayCurtisPlugin {
    fn input(&mut self, filepath: String) -> Result<(), Box<dyn Error>> {
        let param_file = PathBuf::from(&filepath);
        let params = read_parameters(&param_file)?;
        let otu_raw = params
            .get("otufile")
            .ok_or("parameter file missing 'otufile'")?;
        let otu_path = resolve_input_path(&param_file, otu_raw);
        let (samples, otus, counts) = read_otu_csv(otu_path)?;
        self.samples = samples;
        self.otus = otus;
        self.counts = counts;
        self.mapping_path = params
            .get("mapping")
            .map(|v| resolve_input_path(&param_file, v).to_string_lossy().into_owned());
        self.group_column = params.get("column").cloned();
        Ok(())
    }

    fn run(&mut self) -> Result<(), Box<dyn Error>> {
        let n = self.samples.len();
        self.dissimilarity = vec![vec![0.0; n]; n];
        if n < 2 {
            return Ok(());
        }
        if n >= 10 {
            let pairs: Vec<(usize, usize, f64)> = (0..n)
                .into_par_iter()
                .flat_map(|i| {
                    ((i + 1)..n)
                        .map(|j| (i, j, bray_curtis(&self.counts[i], &self.counts[j])))
                        .collect::<Vec<_>>()
                })
                .collect();
            for (i, j, d) in pairs {
                self.dissimilarity[i][j] = d;
                self.dissimilarity[j][i] = d;
            }
        } else {
            for i in 0..n {
                for j in (i + 1)..n {
                    let d = bray_curtis(&self.counts[i], &self.counts[j]);
                    self.dissimilarity[i][j] = d;
                    self.dissimilarity[j][i] = d;
                }
            }
        }
        Ok(())
    }

    fn output(&mut self, filepath: String) -> Result<(), Box<dyn Error>> {
        // `filepath` is treated as an output prefix to match the R plugin's
        // `<prefix>.csv`/`<prefix>.pdf` shape.
        self.write_csv(format!("{filepath}.csv"))?;
        self.write_json(format!("{filepath}.json"))?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Classical 2-D MDS (Torgerson/Gower PCoA) via power iteration
// ---------------------------------------------------------------------------

pub fn classical_mds_2d(d: &[Vec<f64>]) -> (Option<Vec<[f64; 2]>>, Option<[f64; 2]>) {
    let n = d.len();
    if n < 3 {
        return (None, None);
    }
    let mut a = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let v = d[i][j];
            a[i][j] = -0.5 * v * v;
        }
    }
    let row_means: Vec<f64> = (0..n).map(|i| a[i].iter().sum::<f64>() / n as f64).collect();
    let col_means: Vec<f64> = (0..n)
        .map(|j| (0..n).map(|i| a[i][j]).sum::<f64>() / n as f64)
        .collect();
    let grand_mean: f64 = row_means.iter().sum::<f64>() / n as f64;
    let mut b = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            b[i][j] = a[i][j] - row_means[i] - col_means[j] + grand_mean;
        }
    }
    let (l1, v1) = match power_iterate(&b, None, 256) {
        Some(p) => p,
        None => return (None, None),
    };
    let (l2, v2) = match power_iterate(&b, Some(&(l1, v1.clone())), 256) {
        Some(p) => p,
        None => return (None, None),
    };
    let s1 = l1.max(0.0).sqrt();
    let s2 = l2.max(0.0).sqrt();
    let coords: Vec<[f64; 2]> = (0..n).map(|i| [v1[i] * s1, v2[i] * s2]).collect();
    (Some(coords), Some([l1, l2]))
}

fn power_iterate(
    m: &[Vec<f64>],
    deflate: Option<&(f64, Vec<f64>)>,
    max_iter: usize,
) -> Option<(f64, Vec<f64>)> {
    let n = m.len();
    // If we're after the second eigenpair, build the deflated matrix
    // `B' = B - λ₁ v₁ v₁ᵀ` once instead of subtracting in every iteration —
    // numerically stabler than per-step deflation.
    let m_owned: Vec<Vec<f64>>;
    let m_use: &[Vec<f64>] = match deflate {
        None => m,
        Some((l, ev)) => {
            m_owned = (0..n)
                .map(|i| (0..n).map(|j| m[i][j] - l * ev[i] * ev[j]).collect())
                .collect();
            &m_owned
        }
    };

    let mut v = vec![0.0; n];
    for (i, slot) in v.iter_mut().enumerate() {
        *slot = ((i + 1) as f64).sin().abs() + 0.1;
    }
    // For second-eigenpair, project the initial vector perpendicular to v1
    // to avoid hunting back along the dominant direction.
    if let Some((_, ev)) = deflate {
        let proj: f64 = ev.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
        for i in 0..n {
            v[i] -= proj * ev[i];
        }
    }
    normalize(&mut v);

    let mut lambda = 0.0;
    for _ in 0..max_iter {
        let mut w = vec![0.0; n];
        for i in 0..n {
            let mut s = 0.0;
            for j in 0..n {
                s += m_use[i][j] * v[j];
            }
            w[i] = s;
        }
        lambda = dot(&v, &w);
        let norm = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-12 {
            // Genuinely degenerate (rank-deficient B') — return what we have
            // rather than None so callers get a usable 2-D embedding.
            return Some((lambda.max(0.0), v));
        }
        let prev = v.clone();
        let mut diff = 0.0;
        for i in 0..n {
            v[i] = w[i] / norm;
            diff += (v[i] - prev[i]).abs();
        }
        if diff < 1e-10 {
            break;
        }
    }
    Some((lambda, v))
}

fn normalize(v: &mut [f64]) {
    let n: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if n > 0.0 {
        for x in v.iter_mut() {
            *x /= n;
        }
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

// ---------------------------------------------------------------------------
// FFI exports
// ---------------------------------------------------------------------------

#[unsafe(no_mangle)]
pub extern "C" fn BrayCurtis_plugin_create() -> *mut std::ffi::c_void {
    Box::into_raw(Box::new(BrayCurtisPlugin::new())) as *mut std::ffi::c_void
}

#[unsafe(no_mangle)]
pub extern "C" fn BrayCurtis_plugin_destroy(ptr: *mut std::ffi::c_void) {
    if !ptr.is_null() {
        unsafe {
            let _ = Box::from_raw(ptr as *mut BrayCurtisPlugin);
        }
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn BrayCurtis_plugin_input(ptr: *mut std::ffi::c_void, filename: *const c_char) {
    if ptr.is_null() || filename.is_null() {
        return;
    }
    unsafe {
        let plugin = &mut *(ptr as *mut BrayCurtisPlugin);
        let s = CStr::from_ptr(filename).to_str().unwrap_or("").to_string();
        if let Err(e) = plugin.input(s) {
            eprintln!("[BrayCurtis] input error: {e}");
        }
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn BrayCurtis_plugin_run(ptr: *mut std::ffi::c_void) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        let plugin = &mut *(ptr as *mut BrayCurtisPlugin);
        if let Err(e) = plugin.run() {
            eprintln!("[BrayCurtis] run error: {e}");
        }
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn BrayCurtis_plugin_output(ptr: *mut std::ffi::c_void, filename: *const c_char) {
    if ptr.is_null() || filename.is_null() {
        return;
    }
    unsafe {
        let plugin = &mut *(ptr as *mut BrayCurtisPlugin);
        let s = CStr::from_ptr(filename).to_str().unwrap_or("").to_string();
        if let Err(e) = plugin.output(s) {
            eprintln!("[BrayCurtis] output error: {e}");
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_tmp(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        write!(f, "{content}").unwrap();
        f
    }

    #[test]
    fn identical_samples_are_zero() {
        let a = [10.0, 20.0, 30.0];
        assert!(bray_curtis(&a, &a).abs() < 1e-12);
    }

    #[test]
    fn disjoint_samples_are_one() {
        let a = [10.0, 20.0, 0.0, 0.0];
        let b = [0.0, 0.0, 5.0, 7.0];
        assert!((bray_curtis(&a, &b) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn known_textbook_value() {
        // Bray & Curtis (1957) example: a=[6,7,4], b=[10,0,6]
        // |6-10|+|7-0|+|4-6| = 13; sum = 33; d = 13/33
        let a = [6.0, 7.0, 4.0];
        let b = [10.0, 0.0, 6.0];
        assert!((bray_curtis(&a, &b) - 13.0 / 33.0).abs() < 1e-12);
    }

    #[test]
    fn matrix_is_symmetric_with_zero_diagonal() {
        let samples = vec!["A".into(), "B".into(), "C".into(), "D".into()];
        let otus = vec!["s1".into(), "s2".into(), "s3".into()];
        let counts = vec![
            vec![10.0, 0.0, 5.0],
            vec![5.0, 5.0, 0.0],
            vec![0.0, 10.0, 10.0],
            vec![3.0, 3.0, 3.0],
        ];
        let mut p = BrayCurtisPlugin::from_matrix(samples, otus, counts);
        p.run().unwrap();
        let m = p.dissimilarity();
        for i in 0..4 {
            assert!(m[i][i].abs() < 1e-12);
            for j in 0..4 {
                assert!((m[i][j] - m[j][i]).abs() < 1e-12);
                assert!(m[i][j] >= 0.0 && m[i][j] <= 1.0);
            }
        }
    }

    #[test]
    fn plugin_runs_against_csv() {
        let csv = "OTU,A,B,C\nsp1,10,0,5\nsp2,0,20,10\nsp3,5,5,0\n";
        let otu_file = write_tmp(csv);
        let params = format!("otufile {}\n", otu_file.path().display());
        let pf = write_tmp(&params);
        let mut plugin = BrayCurtisPlugin::new();
        plugin.input(pf.path().to_string_lossy().into_owned()).unwrap();
        plugin.run().unwrap();
        assert_eq!(plugin.num_samples(), 3);
        assert!(plugin.dissimilarity()[0][1] > 0.0);
    }

    #[test]
    fn mds_runs_on_small_matrix() {
        let samples = vec!["A".into(), "B".into(), "C".into(), "D".into()];
        let otus = vec!["s1".into(), "s2".into()];
        let counts = vec![
            vec![10.0, 0.0],
            vec![9.0, 1.0],
            vec![0.0, 10.0],
            vec![1.0, 9.0],
        ];
        let mut p = BrayCurtisPlugin::from_matrix(samples, otus, counts);
        p.run().unwrap();
        let (mds, eigs) = classical_mds_2d(p.dissimilarity());
        assert!(mds.is_some());
        assert!(eigs.is_some());
        assert_eq!(mds.unwrap().len(), 4);
    }

    #[test]
    fn resolve_input_path_uses_pluma_prefix_convention() {
        let tmp = tempfile::tempdir().unwrap();
        let prefix = tmp.path();
        let params_dir = prefix.join("parameters");
        let csv_dir = prefix.join("CSV");
        std::fs::create_dir_all(&params_dir).unwrap();
        std::fs::create_dir_all(&csv_dir).unwrap();
        let target = csv_dir.join("x.csv");
        std::fs::write(&target, "OTU\n").unwrap();
        let param_file = params_dir.join("foo.txt");
        std::fs::write(&param_file, "").unwrap();

        assert_eq!(resolve_input_path(&param_file, "CSV/x.csv"), target);
    }

    #[test]
    fn resolve_input_path_preserves_absolute() {
        assert_eq!(
            resolve_input_path(Path::new("/etc/foo.txt"), "/abs/path.csv"),
            PathBuf::from("/abs/path.csv")
        );
    }

    #[test]
    fn ffi_lifecycle_round_trip() {
        let ptr = BrayCurtis_plugin_create();
        assert!(!ptr.is_null());
        BrayCurtis_plugin_destroy(ptr);
    }
}
