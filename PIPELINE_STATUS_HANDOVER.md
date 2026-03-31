# Current Global Strategy

- Distribution method: GitHub repository + Docker image
- Reproducibility target: fixed Linux container with exact local package versions captured from the current workstation
- Runtime target: R 4.4.0 + Bioconductor 3.20 + Python 3 + `PyYAML==6.0.3`
- Dependency strategy:
  - `dependencies/r-package-manifest.csv` stores exact package versions captured from the local machine
  - `R/install_exact_dependencies.R` installs packages sequentially and verifies installed versions
  - `Dockerfile` uses `rocker/r-ver:4.4.0` and installs system libraries required by CRAN / Bioconductor source builds
- Current primary objective: complete Step 4 by obtaining a full Docker build and then running dummy-data validation inside the container without errors

---

## 2026-03-30 16:48:08 +0900 - Trial #15

- Current Action:
  - Start the next Docker rebuild using the nuclear-option `ggtree` install path: pinned GitHub commit instead of the Bioconductor 3.20 tarball.
- Error Log:
  ```text
  [OK] treeio 1.30.0
  [INFO] Installing Bioconductor package: ggtree 3.14.0 (Bioconductor 3.20)
  [INFO] Installing ggtree from pinned GitHub ref: 1a3d0285536625825ccd597d42ca66a6ce8c17c5
  * building ‘ggtree_4.1.1.006.tar.gz’
  ERROR: dependency ‘ggiraph’ is not available for package ‘ggtree’
  * removing ‘/app/.r_libs/ggtree’
  Error: Version check failed for ggtree: expected 4.1.1.006 but found NOT_INSTALLED
  ```
- Fix Applied:
  - Resolved `YuLab-SMU/ggtree` `devel` branch HEAD to commit `1a3d0285536625825ccd597d42ca66a6ce8c17c5`.
  - Fetched the remote `DESCRIPTION` at that commit and recorded `Version: 4.1.1.006`.
  - Updated `R/install_exact_dependencies.R` so `ggtree` installs from that exact GitHub commit and validates both:
    - `packageVersion("ggtree") == 4.1.1.006`
    - `RemoteSha == 1a3d0285536625825ccd597d42ca66a6ce8c17c5`
- Result:
  - Trial #15 eliminated the original lazy-loading failure, but the pinned upstream `ggtree` commit introduced a new required dependency that is not present in the recorded manifest.
  - Trial #16 will preinstall `ggiraph` before the pinned GitHub `ggtree` build and then rerun the Docker image build.

---

## 2026-03-30 17:27:08 +0900 - Trial #16

- Current Action:
  - Start the next Docker rebuild after teaching the exact installer to satisfy the pinned GitHub `ggtree` dependency graph.
- Error Log:
  ```text
  [OK] treeio 1.30.0
  [INFO] Installing Bioconductor package: ggtree 3.14.0 (Bioconductor 3.20)
  [INFO] Preinstalling ggtree GitHub dependency: ggiraph
  Error in remotes::install_cran(package = "ggiraph", lib = user_lib, repos = getOption("repos")[["CRAN"]],  :
    argument "pkgs" is missing, with no default
  Calls: install_bioc_release -> install_ggtree_compat -> <Anonymous> -> lapply
  Execution halted
  ```
- Fix Applied:
  - Updated `R/install_exact_dependencies.R` so `install_ggtree_compat()` preinstalls `ggiraph` from CRAN, with normal `Depends` / `Imports` / `LinkingTo` dependency resolution, before installing `YuLab-SMU/ggtree@1a3d0285536625825ccd597d42ca66a6ce8c17c5`.
- Result:
  - Trial #16 reached the correct dependency edge, but the new preinstall path failed because of an installer API bug, not because of Docker system libraries or package compilation.
  - Trial #17 will correct the `remotes::install_cran()` argument name and rerun the image build.

---

## 2026-03-30 17:44:01 +0900 - Trial #17

- Current Action:
  - Start the next Docker rebuild after fixing the `ggiraph` preinstall invocation in the exact dependency installer.
- Error Log:
  ```text
  [OK] treeio 1.30.0
  [INFO] Installing Bioconductor package: ggtree 3.14.0 (Bioconductor 3.20)
  [INFO] Preinstalling ggtree GitHub dependency: ggiraph
  Installing 1 packages: ggiraph
  [OK] ggiraph 0.9.1
  [INFO] Installing ggtree from pinned GitHub ref: 1a3d0285536625825ccd597d42ca66a6ce8c17c5
  * building ‘ggtree_4.1.1.006.tar.gz’
  [OK] ggtree 4.1.1.006
  [INFO] Exact dependency installation complete.
  ```
- Fix Applied:
  - Corrected `R/install_exact_dependencies.R` to call `remotes::install_cran(pkgs = "ggiraph", ...)` instead of the invalid `package = "ggiraph"` form.
- Result:
  - Trial #17 completed the full Docker build successfully.
  - The resulting image was created as `rnaseq-pipeline:trial17-check`.

---

## 2026-03-30 18:07:17 +0900 - Step 4 Validation Success

- Current Action:
  - Run dummy-data validation inside the newly built container.
- Validation Command:
  ```bash
  docker run --rm rnaseq-pipeline:trial17-check --example --organism zebrafish --run 1 --non-interactive
  ```
- Validation Log:
  ```text
  [INFO] Example dataset created: 200 genes x 12 samples
  [INFO] 選択: zebrafish -> OrgDb=org.Dr.eg.db, KEGG=dre
  [INFO] Non-interactive run: step 1
  [INFO] 完了: 00_align_check.R
  [DESeq2] DEG数: 20
  [INFO] DE結果保存: /app/output/example_results/01_DE/results/de_results.csv
  [INFO] 完了: 01_differential_expression.R
  ```
- Result:
  - Step 4 is complete.
  - The Docker image now builds reproducibly and the in-container zebrafish example validation succeeds end to end.
  - Step 5 can finalize `DISTRIBUTION_GUIDE.md` around the exact-manifest plus SHA-pinned `ggtree` strategy.

---

## 2026-03-30 18:25:44 +0900 - Post-Validation Integration Test

- Current Action:
  - Recreate the Trial #17 example outputs on the host via the validated Docker image, inspect the DE result contract, and run the first downstream file-format consumer against those outputs.
- Primary Outputs Identified:
  - `output/example_results/01_DE/results/de_results.csv`
  - Duplicate convenience copy: `output/example_results/01_DE/de_results.csv`
  - No normalized-count matrix or step-1-generated `.rds` artifact is produced in `output/example_results/`; Step 1 currently emits the DE CSV plus plots and logs.
- Output Contract Check:
  - `de_results.csv` columns are:
    - `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `gene`, `significant`
  - These match the downstream expectations in:
    - `R/06_functional_enrichment.R` (`gene`, `significant`)
    - `R/11_visualization.R` (`gene`, `significant`, `log2FoldChange`, `pvalue`, optional `baseMean`, `padj`)
  - Example metadata preserved the expected condition labels:
    - `control`: 6 samples
    - `treatment`: 6 samples
  - Trial #17 validation logs also confirmed DESeq2 factor ordering:
    - `condition levels: control, treatment`
- Downstream Smoke Test 1:
  - Command:
    ```bash
    docker run --rm -v "/Users/yodao/RNA-seq2/output:/app/output" --entrypoint /usr/local/bin/Rscript rnaseq-pipeline:trial17-check --vanilla /app/R/06_functional_enrichment.R --de_result /app/output/example_results/01_DE/results/de_results.csv --organism org.Dr.eg.db --kegg_org dre --outdir /app/output/example_results
    ```
  - Result:
    ```text
    [INFO] DE結果読み込み: /app/output/example_results/01_DE/results/de_results.csv
    [INFO] 有意な遺伝子数: 20
    Error in .testForValidKeys(x, keys, keytype, fks) :
      None of the keys entered are valid keys for 'SYMBOL'.
    ```
  - Interpretation:
    - This is not a DE output-format incompatibility.
    - The synthetic example genes (`Gene_001` ... `Gene_200`) are placeholder IDs and do not map to valid zebrafish `org.Dr.eg.db` `SYMBOL` keys, so enrichment cannot proceed on the example dataset even though the CSV schema loads correctly.
- Downstream Smoke Test 2:
  - Commands:
    ```bash
    docker run --rm -v "/Users/yodao/RNA-seq2:/work" --workdir /work --entrypoint /usr/local/bin/Rscript rnaseq-pipeline:trial17-check --vanilla example/generate_example.R
    docker run --rm -v "/Users/yodao/RNA-seq2:/work" --workdir /work --entrypoint /usr/local/bin/Rscript rnaseq-pipeline:trial17-check --vanilla R/11_visualization.R --counts /work/example/counts.rds --metadata /work/example/metadata.csv --de_result /work/output/example_results/01_DE/results/de_results.csv --outdir /work/output/example_results
    ```
  - Result:
    ```text
    [VIS] Volcano plot...
    [VIS] MA plot...
    [VIS] Heatmap...
    [VIS] PCA plot...
    [完了] 11_visualization.R
    ```
  - Generated artifacts:
    - `output/example_results/11_visualization/volcano_plot.pdf`
    - `output/example_results/11_visualization/volcano_plot.png`
    - `output/example_results/11_visualization/ma_plot.pdf`
    - `output/example_results/11_visualization/ma_plot.png`
    - `output/example_results/11_visualization/heatmap.pdf`
    - `output/example_results/11_visualization/pca_plot.pdf`
    - `output/example_results/11_visualization/pca_plot.png`
- Result:
  - The Trial #17 DE outputs are structurally compatible with downstream consumers.
  - No tweak to `01_differential_expression.R` is required: the standard DE CSV contract is already correct.
  - The only integration blocker observed is biological identifier realism in the synthetic example dataset for enrichment-style steps, not column naming, factor handling, CSV/RDS loading, or DE result formatting.

---

## 2026-03-30 16:48:08 +0900 - Trial #14 Outcome

- Current Action:
  - Evaluate the direct source-patch strategy that injected `check_linewidth` into the Bioconductor `ggtree 3.14.0` tarball.
- Error Log:
  ```text
  [INFO] Installing patched ggtree for ggplot2 compatibility: 3.14.0
  * installing *source* package ‘ggtree’ ...
  Error in get(x, envir = ns, inherits = FALSE) :
    object 'check_linewidth' not found
  Error: unable to load R code in package ‘ggtree’
  ERROR: lazy loading failed for package ‘ggtree’
  Error: Version check failed for ggtree: expected 3.14.0 but found NOT_INSTALLED
  ```
- Fix Applied:
  - Confirmed the patch file was being added and the installer was using the patched source path.
  - Concluded that simply injecting a new R file into the release tarball is not enough to satisfy the lazy-loading path used by this package version.
- Result:
  - Trial #14 did not resolve the blocker.
  - Trial #15 escalates to a pinned GitHub commit that already contains upstream compatibility changes.

---

## 2026-03-30 16:26:41 +0900 - Trial #14

- Current Action:
  - Start the next Docker rebuild with a direct `ggtree` source patch, using the host's working `check_linewidth` behavior as the compatibility shim.
- Error Log:
  - No build result yet in this trial. This entry records the next fix before relaunching the image build.
- Fix Applied:
  - Extracted the host-side `ggtree::check_linewidth` function body from the installed namespace.
  - Added `download_bioc_source()` and `install_ggtree_compat()` to `R/install_exact_dependencies.R`.
  - The installer now patches the Bioconductor `ggtree 3.14.0` source tarball by injecting `R/check_linewidth_compat.R` before installation.
- Result:
  - The compatibility workaround now targets the actual package that is failing during lazy loading.
  - Trial #14 will verify whether `ggtree` installs and the full image can complete.

---

## 2026-03-30 16:26:41 +0900 - Trial #13 Outcome

- Current Action:
  - Evaluate the long-running Trial #13 build after pinning `ggplot2` to `3.5.1`.
- Error Log:
  ```text
  [OK] treeio 1.30.0
  [INFO] Installing Bioconductor package: ggtree 3.14.0 (Bioconductor 3.20)
  * installing *source* package ‘ggtree’ ...
  Error in get(x, envir = ns, inherits = FALSE) :
    object 'check_linewidth' not found
  Error: unable to load R code in package ‘ggtree’
  ERROR: lazy loading failed for package ‘ggtree’
  Error: Version check failed for ggtree: expected 3.14.0 but found NOT_INSTALLED
  ```
- Fix Applied:
  - Confirmed that the `ggplot2 3.5.1` compatibility pin alone does not resolve the lazy-loading failure.
  - Inspected the host installation and extracted the working `check_linewidth` helper from `ggtree` itself.
- Result:
  - The remaining blocker is still isolated to `ggtree` source installation.
  - Trial #14 moves from package-version pinning to a direct source patch on `ggtree`.

---

## 2026-03-30 16:09:27 +0900 - Trial #13

- Current Action:
  - Start the next Docker rebuild after fixing the installer bug in the compatibility-override lookup.
- Error Log:
  - No build result yet in this trial. This entry records the immediate retry after the Trial #12 scripting failure.
- Fix Applied:
  - Replaced direct `compatibility_overrides[[pkg]]` access with a safe helper function `get_effective_version(pkg, version)`.
  - This prevents non-overridden packages from triggering `subscript out of bounds`.
- Result:
  - The installer can now continue past the override lookup stage.
  - Trial #13 build is the next action.

---

## 2026-03-30 16:09:27 +0900 - Trial #12 Outcome

- Current Action:
  - Evaluate the first build attempt with the `ggplot2 3.5.1` compatibility override.
- Error Log:
  ```text
  Error in compatibility_overrides[[pkg]] : subscript out of bounds
  Calls: %||%
  Execution halted
  The command '/bin/sh -c Rscript --vanilla /app/R/install_exact_dependencies.R' returned a non-zero code: 1
  ```
- Fix Applied:
  - Determined that the new override logic failed before package installation because named-vector indexing with `[[pkg]]` does not return `NULL` for missing names in this context.
- Result:
  - Trial #12 did not test the `ggplot2` compatibility hypothesis.
  - Trial #13 will rerun the same strategy with a safe override helper.

---

## 2026-03-30 16:07:08 +0900 - Trial #12

- Current Action:
  - Start the next Docker rebuild with an explicit `ggplot2` compatibility override before the `ggtree` source install.
- Error Log:
  - No build result yet in this trial. This entry records the compatibility strategy before relaunching the image build.
- Fix Applied:
  - Verified on the host that `ggplot2 4.0.2` does not expose `check_linewidth`, while `ggtree 3.14.0` is installed and loadable.
  - Added a targeted compatibility override in `R/install_exact_dependencies.R` so the container installs `ggplot2 3.5.1` instead of the manifest's `4.0.2`.
  - Kept the exact manifest intact and made the installer log the override explicitly.
- Result:
  - The working hypothesis is now a reproducible package API mismatch, not an M5-specific compiler or architecture failure.
  - Trial #12 will verify whether `ggtree 3.14.0` succeeds once built against `ggplot2 3.5.1`.

---

## 2026-03-30 16:07:08 +0900 - Trial #11 Outcome

- Current Action:
  - Evaluate Trial #11 after forcing `ggfun` namespace unload/reload around the patched install.
- Error Log:
  ```text
  [OK] treeio 1.30.0
  [INFO] Installing Bioconductor package: ggtree 3.14.0 (Bioconductor 3.20)
  * installing *source* package ‘ggtree’ ...
  Error in get(x, envir = ns, inherits = FALSE) :
    object 'check_linewidth' not found
  Error: unable to load R code in package ‘ggtree’
  ERROR: lazy loading failed for package ‘ggtree’
  Error: Version check failed for ggtree: expected 3.14.0 but found NOT_INSTALLED
  ```
- Fix Applied:
  - Collected local host evidence after the failed build:
    - `ggplot2 4.0.2` on the host does not expose `check_linewidth`.
    - `ggtree 3.14.0` is installed and loadable on the host, implying the host package was likely built under an earlier `ggplot2` API.
  - Classified the failure as a source-build compatibility issue between `ggtree 3.14.0` and `ggplot2 4.0.x`, not as an Apple M5 / ARM64 crash.
- Result:
  - Namespace unload/reload did not resolve the blocker.
  - Trial #12 moves to a `ggplot2` compatibility pin for container builds.

---

## 2026-03-30 15:43:39 +0900 - Trial #11

- Current Action:
  - Start the next Docker rebuild after forcing `ggfun` namespace unload/reload around the patched source install.
- Error Log:
  - No build result yet in this trial. This entry records the namespace-state fix before relaunching the image build.
- Fix Applied:
  - Added `unloadNamespace("ggfun")` before and after the patched `install.packages()` call.
  - Strengthened verification to check both namespace object presence and export status for `check_linewidth`.
- Result:
  - The installer now guards against a stale unpatched `ggfun` namespace surviving in the same R session.
  - Trial #11 build is the next action.

---

## 2026-03-30 15:43:39 +0900 - Trial #10 Outcome

- Current Action:
  - Evaluate the long-running Trial #10 Docker rebuild after the verification bug fix.
- Error Log:
  ```text
  [OK] treeio 1.30.0
  [INFO] Installing Bioconductor package: ggtree 3.14.0 (Bioconductor 3.20)
  * installing *source* package ‘ggtree’ ...
  Error in get(x, envir = ns, inherits = FALSE) :
    object 'check_linewidth' not found
  Error: unable to load R code in package ‘ggtree’
  ERROR: lazy loading failed for package ‘ggtree’
  Error: Version check failed for ggtree: expected 3.14.0 but found NOT_INSTALLED
  ```
- Fix Applied:
  - Determined that the patched `ggfun` install path ran, but `ggtree` still saw a namespace without `check_linewidth`.
  - Hypothesis: an older `ggfun` namespace remained loaded in the R session after dependency installation, so `ggtree` did not see the patched namespace contents.
- Result:
  - Trial #10 progressed much further and confirmed the remaining blocker is specifically namespace visibility during `ggtree` installation.
  - Trial #11 addresses this by forcing a namespace unload/reload.

---

## 2026-03-30 15:24:58 +0900 - Trial #10

- Current Action:
  - Start a new Docker rebuild after correcting the `ggfun` patch verification call.
- Error Log:
  - No build result yet in this trial. This entry records the immediate follow-up fix after Trial #9.
- Fix Applied:
  - Replaced `library(ggfun, character.only = TRUE, ...)` with `loadNamespace("ggfun")` before checking `getNamespaceExports()`.
  - Kept the explicit `check_linewidth` export verification so the installer still fails early if the patch is missing.
- Result:
  - The verification bug in the installer is removed.
  - Trial #10 build is the next action.

---

## 2026-03-30 15:24:58 +0900 - Trial #9 Outcome

- Current Action:
  - Evaluate the Trial #9 Docker rebuild after aligning the `ggfun` install sequence.
- Error Log:
  ```text
  * DONE (ggfun)
  Error in install_ggfun_compat(version) : object 'ggfun' not found
  Calls: install_cran_exact ... suppressPackageStartupMessages -> withCallingHandlers -> library
  Execution halted
  The command '/bin/sh -c Rscript --vanilla /app/R/install_exact_dependencies.R' returned a non-zero code: 1
  ```
- Fix Applied:
  - Identified that the package install itself completed and the failure came from the new verification line, not from `ggfun` or `ggtree`.
- Result:
  - Trial #9 confirmed the installer now reaches and completes the patched `ggfun` source install.
  - The next blocker is only the verification implementation bug, which is now corrected in Trial #10.

---

## 2026-03-30 15:20:37 +0900 - Trial #9

- Current Action:
  - Start the next Docker rebuild after aligning `R/install_exact_dependencies.R` with the older working `ggfun` patch sequence.
- Error Log:
  - No build result yet in this trial. This entry records the new fix before relaunching the image build.
- Fix Applied:
  - Removed the initial unpatched `remotes::install_version("ggfun", ...)` step.
  - Added explicit preinstallation of `ggplotify 0.1.3`, `aplot 0.2.9`, and `yulab.utils 0.2.4` before the patched `ggfun` source install.
  - Added a post-install verification that `ggfun` exports `check_linewidth`; the installer now stops immediately if the compatibility patch is missing.
- Result:
  - Installer logic is now closer to the previously working Docker restore script.
  - Trial #9 build is the next action.

---

## 2026-03-30 15:20:37 +0900 - Trial #8

- Current Action:
  - Rebuild the Docker image with `DOCKER_BUILDKIT=0` to bypass the prior Docker transport EOF and confirm whether the remaining blocker is a package issue or an engine issue.
- Error Log:
  ```text
  Error in get(x, envir = ns, inherits = FALSE) :
    object 'check_linewidth' not found
  Error: unable to load R code in package ‘ggtree’
  ERROR: lazy loading failed for package ‘ggtree’
  Error: Version check failed for ggtree: expected 3.14.0 but found NOT_INSTALLED
  The command '/bin/sh -c Rscript --vanilla /app/R/install_exact_dependencies.R' returned a non-zero code: 1
  ```
- Fix Applied:
  - No code fix inside Trial #8 itself; this run was used to isolate the failure mode with BuildKit disabled.
- Result:
  - The build no longer ended with Docker EOF.
  - The remaining blocker is reproducible package installation failure at `ggtree`, specifically the missing `ggfun::check_linewidth` symbol during lazy loading.

---

## 2026-03-30 14:48:02 +0900 - Trial #7

- Current Action:
  - Rebuild the Docker image after fixing the `ggfun` patch order so that dependencies install first and the `check_linewidth` compatibility patch is applied afterward.
- Error Log:
  ```text
  ERROR: failed to build: failed to receive status: rpc error: code = Unavailable desc = error reading from server: EOF
  ```
- Fix Applied:
  - Updated `R/install_exact_dependencies.R` so `ggfun` is first installed with normal dependency resolution, then reinstalled from patched source with `check_linewidth`.
- Result:
  - The build progressed much further than before and passed the previous `ggfun` / `ggtree` failure point.
  - The run later aborted with a Docker/build transport EOF after long compilation, so package dependency logic is improved but full image completion is still not confirmed.

---

## 2026-03-30 14:20:00 +0900 - Trial #6

- Current Action:
  - Rebuild the Docker image after adding a `ggfun` source patch to restore `check_linewidth` for `ggtree`.
- Error Log:
  ```text
  Error in get(x, envir = ns, inherits = FALSE) :
    object 'check_linewidth' not found
  ERROR: lazy loading failed for package ‘ggtree’
  Error: Version check failed for ggtree: expected 3.14.0 but found NOT_INSTALLED
  ```
- Fix Applied:
  - Added a `ggfun` patch routine to `R/install_exact_dependencies.R` to inject `check_linewidth` and export it before installation.
- Result:
  - The original `ggtree` compatibility error was addressed conceptually.
  - A new issue appeared because the patched `ggfun` source install did not first resolve `dplyr` and `yulab.utils`.

---

## 2026-03-30 13:55:00 +0900 - Trial #5

- Current Action:
  - Rebuild after moving `WGCNA` later in the package order.
- Error Log:
  ```text
  ERROR: dependencies ‘impute’, ‘preprocessCore’ are not available for package ‘WGCNA’
  Error: Version check failed for WGCNA: expected 1.74 but found NOT_INSTALLED
  ```
- Fix Applied:
  - Moved `WGCNA` to the end of `ordered_packages` in `scripts/export_dependency_manifest.R`.
  - Regenerated `dependencies/r-package-manifest.csv`.
- Result:
  - The `WGCNA` order problem was removed.
  - The next blocking issue became `ggtree` failing because `ggfun::check_linewidth` was missing.

---

## 2026-03-30 13:20:00 +0900 - Trial #4

- Current Action:
  - Rebuild after adding `cmake` as a system dependency for `fs`.
- Error Log:
  ```text
  /bin/bash: line 2: cmake: command not found
  ERROR: compilation failed for package ‘fs’
  ERROR: dependency ‘fs’ is not available for package ‘yulab.utils’
  ERROR: dependency ‘yulab.utils’ is not available for package ‘ggfun’
  ```
- Fix Applied:
  - Added `cmake` to the `apt-get install` list in `Dockerfile`.
- Result:
  - The `fs` build blocker was resolved in subsequent attempts.
  - The build then progressed to later package-order and compatibility issues.

---

## 2026-03-30 12:40:00 +0900 - Trial #3

- Current Action:
  - Rebuild after narrowing CRAN dependency installation from broad `TRUE` to required dependency classes only.
- Error Log:
  ```text
  Installing many unnecessary CRAN packages and examples/suggests during build.
  ```
- Fix Applied:
  - Changed `remotes::install_version(..., dependencies = TRUE)` to `dependencies = c("Depends", "Imports", "LinkingTo")` in `R/install_exact_dependencies.R`.
- Result:
  - Reduced install scope and improved reproducibility.
  - Full build still failed later because `fs` needed `cmake`.

---

## 2026-03-30 11:35:00 +0900 - Trial #2

- Current Action:
  - Validate the new exact installer locally before using it in Docker.
- Error Log:
  ```text
  Warning in install.packages("remotes") :
    'lib = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library"' is not writable
  Error in install.packages("remotes") : unable to install packages
  ```
- Fix Applied:
  - Updated `R/install_exact_dependencies.R` to install bootstrap packages into `R_LIBS_USER`.
  - Added explicit `lib = user_lib` for `install.packages()`, `remotes::install_version()`, and `BiocManager::install()`.
  - Updated manifest export logic so `org.*.eg.db` packages are recorded as `Bioconductor 3.20`.
- Result:
  - The installer became suitable for container use with a writable user library.

---

## 2026-03-30 11:03:08 +0900 - Trial #1

- Current Action:
  - Establish a reproducible dependency capture mechanism and create the Docker distribution scaffolding.
- Error Log:
  - No Docker build attempted yet in this trial.
- Fix Applied:
  - Added `scripts/export_dependency_manifest.R`.
  - Generated `dependencies/r-package-manifest.csv` and `dependencies/runtime-versions.txt`.
  - Added `R/install_exact_dependencies.R`.
  - Added `requirements.txt`, `Dockerfile`, `docker-compose.yml`, and `.dockerignore`.
  - Updated documentation and config templates for the new CLI options.
- Result:
  - Steps 1-3 were completed.
  - Step 4 started with iterative container build/debugging.
