#!/usr/bin/env python3
import argparse
import os
import platform
import re
import subprocess
from pathlib import Path
from datetime import datetime
from importlib import metadata


OPENMP_STEPS = {5, 6, 7, 13, 14, 15}


def run(cmd, log_path):
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(proc.stdout or "")
        if proc.stderr:
            f.write("\n[STDERR]\n")
            f.write(proc.stderr)
    return proc.returncode, proc.stdout + "\n" + proc.stderr


def ensure_r_dependencies(log_path):
    cmd = ["Rscript", "--vanilla", "R/install_missing_packages.R"]
    return run(cmd, log_path)


def check_openmp():
    script = (
        "res <- tryCatch({"
        "suppressMessages(library(WGCNA)); TRUE"
        "}, error = function(e) {"
        "if (grepl('OMP: Error #179', e$message)) FALSE else stop(e)"
        "});"
        "if (isTRUE(res)) {cat('OPENMP_OK\\n')} else {cat('OPENMP_FAIL: OMP: Error #179\\n')}"
    )
    proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True)
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    out = out.strip()
    if proc.returncode == 0 and "OPENMP_OK" in proc.stdout:
        return True, ""
    if "OMP: Error #179" in out:
        return False, out
    raise RuntimeError(out or f"OpenMP check failed (exit {proc.returncode})")


def detect_r_packages(r_dir):
    pkgs = set()
    pattern = re.compile(r"(?:library|require)\(([^)]+)\)")
    for path in Path(r_dir).glob("*.R"):
        for line in path.read_text(encoding="utf-8").splitlines():
            m = pattern.search(line)
            if not m:
                continue
            raw = m.group(1).strip()
            raw = raw.strip("\"'")
            if raw.startswith("opt$"):
                continue
            if raw:
                pkgs.add(raw)
    return sorted(pkgs)


def detect_python_packages():
    return ["pyyaml"]


def r_pkg_versions(pkgs):
    if not pkgs:
        return {}
    pkg_vec = "c(" + ",".join([f'"{p}"' for p in pkgs]) + ")"
    script = (
        f"pkgs <- {pkg_vec}; "
        "ip <- installed.packages(); "
        "for (p in pkgs) {"
        "if (p %in% rownames(ip)) {"
        "cat(p, ':', ip[p, 'Version'], '\\n', sep='')"
        "} else {"
        "cat(p, ':', 'NOT_INSTALLED', '\\n', sep='')"
        "}"
        "}"
    )
    proc = subprocess.run(["Rscript", "-e", script], capture_output=True, text=True)
    versions = {}
    for line in (proc.stdout or "").splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            versions[k.strip()] = v.strip()
    return versions


def python_pkg_versions(pkgs):
    versions = {}
    for p in pkgs:
        try:
            versions[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            versions[p] = "NOT_INSTALLED"
    return versions


def build_steps(args):
    outdir = Path(args.outdir)
    steps = [
        (0, ["Rscript", "R/00_align_check.R", "-c", args.counts, "-m", args.metadata, "-o", str(outdir)]),
        (1, ["Rscript", "R/01_differential_expression.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
             "--tool", "deseq2", "--condition_col", args.condition_col,
             "--ref_level", args.ref_level, "--treat_level", args.treat_level]),
        (2, ["Rscript", "R/02_sample_clustering.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir)]),
        (3, ["Rscript", "R/03_dimension_reduction.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
             "--color_by", args.color_by]),
        (4, ["Rscript", "R/04_normalization.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir)]),
        ("05_clustering", ["Rscript", "R/05_clustering.R",
                           "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
                           "--normalization_method", "vst", "--cluster_method", "both"]),
        (5, ["Rscript", "R/05_expression_pattern.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
             "--trait_col", args.condition_col]),
        (6, ["Rscript", "R/06_functional_enrichment.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
             "--organism", args.organism, "--kegg_org", args.kegg_org]),
        (7, ["Rscript", "R/07_gsea.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
             "--organism", args.organism, "--kegg_org", args.kegg_org,
             "--geneset", args.gsea_geneset, "--min_gs", str(args.gsea_min_gs),
             "--max_gs", str(args.gsea_max_gs)]),
        (8, ["Rscript", "R/08_batch_correction.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
             "--batch_col", args.batch_col, "--condition_col", args.condition_col]),
        (9, ["Rscript", "R/09_regression.R",
             "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
             "--formula_str", args.formula_str]),
        (10, ["Rscript", "R/10_sample_qc.R",
              "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir)]),
        (11, ["Rscript", "R/11_visualization.R",
              "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
              "--color_by", args.color_by]),
        (12, ["Rscript", "R/12_phenotype_association.R",
              "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
              "--phenotype_col", args.phenotype_col]),
        (13, ["Rscript", "R/13_machine_learning.R",
              "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir),
              "--target_col", args.target_col]),
        (14, ["Rscript", "R/14_network_analysis.R",
              "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir)]),
        (15, ["Rscript", "R/15_single_cell_integration.R",
              "--counts", args.counts, "--metadata", args.metadata, "--outdir", str(outdir)]),
    ]
    return steps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--condition-col", dest="condition_col", default="condition")
    ap.add_argument("--ref-level", dest="ref_level", default="control")
    ap.add_argument("--treat-level", dest="treat_level", default="treatment")
    ap.add_argument("--batch-col", dest="batch_col", default="batch")
    ap.add_argument("--color-by", dest="color_by", default="condition")
    ap.add_argument("--phenotype-col", dest="phenotype_col", default="condition")
    ap.add_argument("--target-col", dest="target_col", default="condition")
    ap.add_argument("--formula-str", dest="formula_str", default="~ condition")
    ap.add_argument("--organism", default="org.Hs.eg.db")
    ap.add_argument("--kegg-org", dest="kegg_org", default="hsa")
    ap.add_argument("--gsea-geneset", dest="gsea_geneset", default="all")
    ap.add_argument("--gsea-min-gs", dest="gsea_min_gs", type=int, default=15)
    ap.add_argument("--gsea-max-gs", dest="gsea_max_gs", type=int, default=500)
    args = ap.parse_args()

    user_lib = Path(__file__).resolve().parent.parent / ".r_libs"
    user_lib.mkdir(parents=True, exist_ok=True)
    os.environ["R_LIBS_USER"] = str(user_lib)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_dir = outdir / "validation_logs"
    log_dir.mkdir(exist_ok=True)

    install_rc, _ = ensure_r_dependencies(log_dir / "install_r_packages.log")

    openmp_available = None
    openmp_error = ""
    try:
        openmp_available, openmp_error = check_openmp()
    except Exception as exc:
        openmp_available, openmp_error = None, str(exc)

    results = []
    for step, cmd in build_steps(args):
        step_slug = f"{step:02d}" if isinstance(step, int) else str(step)
        log_path = log_dir / f"step_{step_slug}.log"
        if openmp_available is False and step in OPENMP_STEPS:
            results.append((step, "SKIPPED", "OpenMP unavailable"))
            continue
        rc, _ = run(cmd, log_path)
        status = "OK" if rc == 0 else f"FAIL (exit {rc})"
        results.append((step, status, str(log_path)))

    # Dependency audit
    r_pkgs = detect_r_packages("R")
    py_pkgs = detect_python_packages()
    r_versions = r_pkg_versions(r_pkgs)
    py_versions = python_pkg_versions(py_pkgs)

    # R version
    r_ver = subprocess.run(["Rscript", "-e", "cat(R.version.string)"], capture_output=True, text=True).stdout.strip()

    report_path = outdir / "final_validation_report.txt"
    with report_path.open("w", encoding="utf-8") as f:
        f.write("Final Validation Report\n")
        f.write("========================\n")
        f.write(f"Timestamp: {datetime.now().isoformat()}\n")
        f.write(f"OS: {platform.platform()}\n")
        f.write(f"Python: {platform.python_version()}\n")
        f.write(f"Python executable: {os.path.realpath(os.sys.executable)}\n")
        f.write(f"R: {r_ver}\n")
        f.write(f"R_LIBS_USER: {os.environ.get('R_LIBS_USER', '')}\n")
        f.write(f"R dependency install status: {'OK' if install_rc == 0 else 'FAIL'}\n")
        f.write("R dependency install log: results/validation_logs/install_r_packages.log\n\n")

        if openmp_available is True:
            f.write("OpenMP availability: TRUE\n")
        elif openmp_available is False:
            f.write("OpenMP availability: FALSE\n")
        else:
            f.write("OpenMP availability: UNKNOWN (check failed)\n")

        if openmp_available is False or openmp_available is None:
            f.write("OpenMP error: " + openmp_error + "\n")
            if openmp_available is False:
                f.write("Skipped steps: " + ", ".join(f"{s:02d}" for s in sorted(OPENMP_STEPS)) + "\n\n")

        f.write("Detected R packages:\n")
        for p in r_pkgs:
            f.write(f"- {p}: {r_versions.get(p, 'UNKNOWN')}\n")
        f.write("\nDetected Python packages:\n")
        for p in py_pkgs:
            f.write(f"- {p}: {py_versions.get(p, 'UNKNOWN')}\n")

        f.write("\nExternal tools:\n")
        f.write("- Rscript (required)\n")
        f.write("- conda (recommended for reproducibility)\n")

        f.write("\nR script validation results:\n")
        for step, status, info in results:
            if isinstance(step, int):
                step_label = f"{step:02d}"
            else:
                step_label = str(step)
            f.write(f"- Step {step_label}: {status} ({info})\n")

    print(f"[INFO] Wrote {report_path}")


if __name__ == "__main__":
    main()
