## 2026-03-31

### Enhancements: Per-run outputs, UID/GID-safe docker, metadata logging

```diff
diff --git a/main_pipeline.py b/main_pipeline.py
+ DEFAULT_OUTPUT_ROOT = Path.home() / "Output"
+ def get_output_root() -> Path
+ def create_run_dir() -> Path  # creates run_YYYY-MM-DD_HH-MM-SS with collision suffix
+ resolve_outdir() now returns a per-run directory under the output root
+ PIPELINE_OUTPUT_DIR env is injected into all R script calls so defaults follow the run dir
+ write_run_metadata() writes run_metadata.yaml (timestamp, CLI, organism, paths, run selection)
+ Logging/outputs now target the per-run directory; example mode respects it

diff --git a/R/00_load_data.R b/R/00_load_data.R
+ pipeline_output_root now honors PIPELINE_OUTPUT_DIR or requested --outdir
+ ensure_pipeline_output_root() updates the global root and uses the caller’s outdir/run dir

diff --git a/docker-compose.yml b/docker-compose.yml
+ user: "${UID:-0}:${GID:-0}" for host-permission alignment
+ PIPELINE_OUTPUT_DIR=/app/Output and ${HOME}/Output mount retained

diff --git a/README.md b/README.md
+ Documented per-run `run_YYYY-MM-DD_HH-MM-SS` directories, PIPELINE_OUTPUT_DIR, UID/GID usage

diff --git a/config/config_template.yaml b/config/config_template.yaml
+ Output comment mentions `~/Output` or `PIPELINE_OUTPUT_DIR/run_YYYY-MM-DD_HH-MM-SS`
```

Run directory design: every execution creates a unique `run_YYYY-MM-DD_HH-MM-SS` under the resolved output root (PIPELINE_OUTPUT_DIR if set, otherwise `~/Output`). If a collision occurs, a numeric suffix is added. R scripts inherit the run directory via `PIPELINE_OUTPUT_DIR`, keeping default paths consistent across steps. Logs (`pipeline.log`), run metadata (`run_metadata.yaml`), and analysis outputs live inside the run directory, preventing overwrites of previous results.

UID/GID handling: `docker-compose.yml` now supports `user: "${UID:-0}:${GID:-0}"`. Users can export `UID=$(id -u)` and `GID=$(id -g)` before `docker compose run` to avoid root-owned files on host volumes while remaining optional for existing setups.

Verification:
- `python3 main_pipeline.py --example --run 1` (escalated to create ~/Output) → success.
- Outputs written to `/Users/yodao/Output/run_2026-03-31_16-16-49/` including `pipeline.log`, `run_metadata.yaml`, `environment_report.txt`, and DE results/plots under `01_DE/`.
- No quarantine failures or permission errors observed; DESeq2 completed with DEG数: 20.

### Code Diff

```diff
diff --git a/R/00_load_data.R b/R/00_load_data.R
index 7bbea87..a8f9ab3 100644
--- a/R/00_load_data.R
+++ b/R/00_load_data.R
@@ -44,16 +44,20 @@ pipeline_root_dir <- local({
   }
 })
 
-pipeline_output_root <- normalizePath(
-  file.path(pipeline_root_dir, "output", "result"),
-  mustWork = FALSE
-)
+pipeline_output_root <- local({
+  configured_root <- Sys.getenv("PIPELINE_OUTPUT_DIR", unset = "")
+  if (nzchar(configured_root)) {
+    normalizePath(configured_root, mustWork = FALSE)
+  } else {
+    normalizePath(path.expand("~/Output"), mustWork = FALSE)
+  }
+})
 
 #' コマンドライン引数をパースする基本オプションリスト
 base_option_list <- list(
   make_option("--counts",   type = "character", help = "Count matrix RDS path"),
   make_option("--metadata", type = "character", help = "Metadata xlsx/csv path"),
-  make_option("--outdir",   type = "character", default = "", help = "Output directory (retained for compatibility; results are written under output/result)"),
+  make_option("--outdir",   type = "character", default = "", help = "Output directory (retained for compatibility; results are written under ~/Output or PIPELINE_OUTPUT_DIR)"),
   make_option("--label_mode", type = "character", default = "sampleID_only",
               help = "Sample label display: sampleID_only | symbol_only | both"),
   make_option("--symbol_col", type = "character", default = "",
diff --git a/R/install_missing_packages.R b/R/install_missing_packages.R
index 7350b8d..0367fcc 100644
--- a/R/install_missing_packages.R
+++ b/R/install_missing_packages.R
@@ -112,13 +112,48 @@ quarantine_package <- function(pkg, reason) {
   dir.create(quarantine_root, recursive = TRUE, showWarnings = FALSE)
   stamp <- format(Sys.time(), "%Y%m%d%H%M%S")
   target <- file.path(quarantine_root, paste0(pkg, "-", stamp))
-  ok <- file.rename(pkg_dir, target)
-  if (ok) {
+  rename_ok <- suppressWarnings(file.rename(pkg_dir, target))
+  if (isTRUE(rename_ok)) {
     message("[WARN] Quarantined incompatible package from user library: ", pkg)
     message("[WARN]   Reason: ", reason)
     message("[WARN]   Moved to: ", target)
     return(invisible(TRUE))
   }
+
+  # Docker OverlayFS and bind mounts can make file.rename() fail with EXDEV.
+  copy_ok <- suppressWarnings(file.copy(pkg_dir, target, recursive = TRUE, copy.mode = TRUE, copy.date = TRUE))
+  if (isTRUE(copy_ok) && dir.exists(target)) {
+    delete_ok <- suppressWarnings(unlink(pkg_dir, recursive = TRUE, force = TRUE))
+    if (identical(delete_ok, 0L) && !dir.exists(pkg_dir)) {
+      message("[WARN] Quarantined incompatible package from user library: ", pkg)
+      message("[WARN]   Reason: ", reason)
+      message("[WARN]   file.rename() failed; used copy-and-delete fallback.")
+      message("[WARN]   Copied to: ", target)
+      return(invisible(TRUE))
+    }
+
+    message("[WARN] Package copy succeeded but source cleanup was incomplete: ", pkg)
+    message("[WARN]   Attempting direct removal of the original package directory.")
+    cleanup_ok <- suppressWarnings(unlink(pkg_dir, recursive = TRUE, force = TRUE))
+    if (identical(cleanup_ok, 0L) && !dir.exists(pkg_dir)) {
+      message("[WARN] Quarantined incompatible package from user library: ", pkg)
+      message("[WARN]   Reason: ", reason)
+      message("[WARN]   Original package removed after copy fallback.")
+      message("[WARN]   Copied to: ", target)
+      return(invisible(TRUE))
+    }
+  }
+
+  unlink(target, recursive = TRUE, force = TRUE)
+
+  final_cleanup <- suppressWarnings(unlink(pkg_dir, recursive = TRUE, force = TRUE))
+  if (identical(final_cleanup, 0L) && !dir.exists(pkg_dir)) {
+    message("[WARN] Quarantined incompatible package from user library: ", pkg)
+    message("[WARN]   Reason: ", reason)
+    message("[WARN]   Quarantine move failed; removed package directly to stop shadowing.")
+    return(invisible(TRUE))
+  }
+
   message("[WARN] Failed to quarantine package: ", pkg)
   invisible(FALSE)
 }
diff --git a/config/config_template.yaml b/config/config_template.yaml
index 6980abe..0bd95ca 100644
--- a/config/config_template.yaml
+++ b/config/config_template.yaml
@@ -15,7 +15,7 @@ metadata_path: ''        # Path to metadata (.xlsx/.csv), e.g., /data/meta.xlsx
 organism: 'org.Hs.eg.db'
 
 # --- Output ---
-outdir: ''               # Retained for compatibility; results are always written to ./output/result
+outdir: ''               # Retained for compatibility; results are always written to ~/Output or PIPELINE_OUTPUT_DIR
diff --git a/docker-compose.yml b/docker-compose.yml
index a9a167b..b167e41 100644
--- a/docker-compose.yml
+++ b/docker-compose.yml
@@ -7,8 +7,9 @@ services:
     working_dir: /app
     environment:
       R_LIBS_USER: /app/.r_libs
+      PIPELINE_OUTPUT_DIR: /app/Output
     volumes:
       - ./:/app
       - ./data:/data
-      - ./output:/output
+      - ${HOME}/Output:/app/Output
     command: ["--help"]
diff --git a/main_pipeline.py b/main_pipeline.py
index ba94f12..8d4af19 100644
--- a/main_pipeline.py
+++ b/main_pipeline.py
@@ -27,7 +27,7 @@ from pathlib import Path
 BASE_DIR = Path(__file__).resolve().parent
 R_DIR = BASE_DIR / "R"
 CONFIG_FILE = BASE_DIR / "config.yaml"
-RESULT_ROOT = BASE_DIR / "output" / "result"
+DEFAULT_OUTPUT_ROOT = Path.home() / "Output"
@@ -210,6 +210,13 @@ R_DEPS_INSTALLED = False
 R_LIBS_USER = BASE_DIR / ".r_libs"
 
 
+def get_output_root():
+    configured = os.environ.get("PIPELINE_OUTPUT_DIR", "").strip()
+    if configured:
+        return Path(configured).expanduser().resolve()
+    return DEFAULT_OUTPUT_ROOT.resolve()
+
+
 def init_rscript():
@@ -418,11 +425,12 @@ def get_param(cfg, args, key, default):
 
 
 def resolve_outdir(base_outdir):
+    result_root = get_output_root()
     if base_outdir:
         requested = Path(base_outdir).resolve()
     else:
-        requested = RESULT_ROOT
-    outdir = RESULT_ROOT
+        requested = result_root
+    outdir = result_root
     outdir.mkdir(parents=True, exist_ok=True)
     return outdir, requested
@@ -957,7 +965,7 @@ def main():
     parser.add_argument("--non-interactive", action="store_true",
                         help="Disable prompts and require CLI values")
     parser.add_argument("--run", type=str, help="Run step: all, 1, 6, or 7")
-    parser.add_argument("--outdir", type=str, help="Output directory (retained for compatibility; results are written to output/result)")
+    parser.add_argument("--outdir", type=str, help="Output directory (retained for compatibility; results are written to ~/Output or PIPELINE_OUTPUT_DIR)")
@@ -1043,7 +1051,7 @@ def main():
         if not args.run:
             args.run = "1"  # default: run DE only for quick demo
         if not args.outdir:
-            args.outdir = str(RESULT_ROOT)
+            args.outdir = str(get_output_root())
 ```

### Why `file.rename()` Failed

`file.rename()` is implemented as a filesystem move. In Docker OverlayFS, bind mounts, and host-mounted volumes, the source package directory and quarantine directory can appear under different device IDs even when they look adjacent in the container path tree. That causes the underlying rename syscall to return `EXDEV` ("Cross-device link"), which R surfaces as a failed `file.rename()`. The failure left incompatible packages in `.r_libs`, so they continued to shadow the valid system libraries. The fix treats rename as an optimization only, then falls back to `file.copy(..., recursive=TRUE)` plus `unlink(..., recursive=TRUE)` and finally to a direct `unlink()` cleanup path if quarantine relocation itself cannot be completed.

### Verification Results

Commands run:

```bash
python3 main_pipeline.py --help
python3 main_pipeline.py --example --run 1
```

Observed results:

- `python3 main_pipeline.py --help` exited successfully.
- `python3 main_pipeline.py --example --run 1` exited successfully.
- No `Failed to quarantine` warnings were present in `/Users/yodao/Output/pipeline.log`.
- DE analysis completed with `DESeq2` and reported `DEG数: 20`.
- Output artifacts were created under `/Users/yodao/Output`, including:
  - `/Users/yodao/Output/01_DE/results/de_results.csv`
  - `/Users/yodao/Output/environment_report.txt`
  - `/Users/yodao/Output/pipeline.log`

### Output Persistence Guarantee

The pipeline now resolves its canonical output root from `PIPELINE_OUTPUT_DIR` when set, otherwise from `~/Output`. For local runs, this writes directly to the host home directory. For Docker Compose, `docker-compose.yml` now mounts `${HOME}/Output:/app/Output` and sets `PIPELINE_OUTPUT_DIR=/app/Output`. That means the container writes into a bind-mounted host directory, so results persist independently of the container filesystem and remain available even if the container or image is deleted.
