#!/usr/bin/env python3
"""
RNA-seq 統合解析パイプライン
メインUI および R スクリプト起動制御
"""

import os
import sys
import subprocess
try:
    import yaml
except ModuleNotFoundError:
    print("[ERROR] Python module 'pyyaml' is missing. Install it or use the conda environment.")
    print("        Example: conda env create -f environment.yml && conda activate rnaseq_env")
    sys.exit(1)
import shutil
import argparse
import logging
from datetime import datetime
from pathlib import Path


# ─────────────────────────────────────────────
# ユーティリティ
# ─────────────────────────────────────────────

BASE_DIR = Path(__file__).resolve().parent
R_DIR = BASE_DIR / "R"
CONFIG_FILE = BASE_DIR / "config.yaml"
DEFAULT_OUTPUT_ROOT = Path.home() / "Output"

LOG = logging.getLogger("pipeline")

STEP_DESCRIPTIONS = {
    1: {
        "title": "Differential Expression",
        "en": "Identifies genes with significant expression differences between conditions.",
        "jp": "条件間で発現量が有意に異なる遺伝子を同定する解析です。",
    },
    6: {
        "title": "Functional Enrichment",
        "en": "Tests over-representation of functional categories such as GO and KEGG among DE genes.",
        "jp": "差次的発現遺伝子に対してGOやKEGGなどの機能カテゴリの富化を解析します。",
    },
    7: {
        "title": "GSEA",
        "en": "Ranks genes by expression change and tests enrichment of predefined gene sets.",
        "jp": "遺伝子を発現変化でランキングし、既知の遺伝子セットの富化を評価します。",
    },
}

STEP_NAMES = {
    1: "Differential Expression",
    2: "Sample Clustering",
    3: "Dimension Reduction",
    4: "Normalization",
    5: "Expression Pattern (WGCNA)",
    6: "Functional Enrichment",
    7: "GSEA",
    8: "Batch Correction",
    9: "Regression",
    10: "Sample QC",
    11: "Visualization",
    12: "Phenotype Association",
    13: "Machine Learning",
    14: "Network Analysis",
    15: "Single-cell Integration",
}

STEP_ALIASES = {
    "1": 1, "de": 1, "diffexp": 1, "differential_expression": 1,
    "6": 6, "enrichment": 6, "functional_enrichment": 6,
    "7": 7, "gsea": 7,
    "all": "all",
}

OPENMP_STEPS = {5, 6, 7, 13, 14, 15}
OPENMP_AVAILABLE = True
OPENMP_ERROR = ""
def run_r_script(script_path, args):
    # .Rprofile が存在する /app を作業ディレクトリにする
    result = subprocess.run(
        ["Rscript", script_path] + args,
        cwd="/app", 
        capture_output=True,
        text=True
    )
    return result

def setup_logging(outdir=None):
    LOG.setLevel(logging.INFO)
    LOG.handlers.clear()
    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    LOG.addHandler(ch)

    if outdir:
        log_path = Path(outdir) / "pipeline.log"
        fh = logging.FileHandler(log_path, encoding="utf-8")
        fh.setFormatter(formatter)
        LOG.addHandler(fh)

# 生物種の候補と派生情報
ORG_OPTIONS = {
    "human": {"orgdb": "org.Hs.eg.db", "kegg": "hsa", "label": "Human"},
    "mouse": {"orgdb": "org.Mm.eg.db", "kegg": "mmu", "label": "Mouse"},
    "zebrafish": {"orgdb": "org.Dr.eg.db", "kegg": "dre", "label": "Zebrafish"},
}

# グローバル保持（GSEA/Enrichment 用）
ORG_ARGS = {}


def ensure_dirs(outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    for sub in [
        "01_DE", "02_clustering", "03_dimreduc", "04_normalization",
        "05_expression_pattern", "06_enrichment", "07_gsea",
        "08_batch_correction", "09_regression", "10_qc",
        "11_visualization", "12_phenotype", "13_ml",
        "14_network", "15_singlecell"
    ]:
        (outdir / sub).mkdir(exist_ok=True)


def find_rscript():
    """Rscript のパスを検索"""
    rscript = shutil.which("Rscript")
    if rscript is None:
        # Windows 典型パス
        candidates = [
            r"C:\Program Files\R\R-4.3.1\bin\Rscript.exe",
            r"C:\Program Files\R\R-4.4.0\bin\Rscript.exe",
            r"C:\Program Files\R\R-4.2.3\bin\Rscript.exe",
            "/usr/bin/Rscript",
            "/usr/local/bin/Rscript",
        ]
        for c in candidates:
            if os.path.isfile(c):
                return c
        LOG.error("Rscript が見つかりません。PATHを確認してください。")
        sys.exit(1)
    return rscript


def check_openmp_runtime():
    rscript = find_rscript()
    script = (
        "res <- tryCatch({"
        "suppressMessages(library(WGCNA)); TRUE"
        "}, error = function(e) {"
        "if (grepl('OMP: Error #179', e$message)) FALSE else stop(e)"
        "});"
        "if (isTRUE(res)) {cat('OPENMP_OK\\n')} else {cat('OPENMP_FAIL: OMP: Error #179\\n')}"
    )
    try:
        env = os.environ.copy()
        env["R_LIBS_USER"] = str(R_LIBS_USER)
        proc = subprocess.run(
            [rscript, "-e", script],
            capture_output=True,
            text=True,
            check=False,
            env=env,
        )
    except Exception as exc:
        return False, f"OpenMP check failed: {exc}"

    output = (proc.stdout or "") + "\n" + (proc.stderr or "")
    output = output.strip()
    if proc.returncode == 0 and "OPENMP_OK" in proc.stdout:
        return True, ""
    if "OMP: Error #179" in output:
        return False, output
    if "OPENMP_FAIL" in output:
        return False, output
    raise RuntimeError(output or f"OpenMP check failed (exit {proc.returncode})")


def write_environment_report(outdir, openmp_available, openmp_error):
    report_path = Path(outdir) / "environment_report.txt"
    skipped = sorted(OPENMP_STEPS)
    with report_path.open("w", encoding="utf-8") as f:
        f.write("OpenMP availability: " + ("TRUE" if openmp_available else "FALSE") + "\n")
        if not openmp_available:
            f.write("Skipped steps: " + ", ".join(f"{n:02d}" for n in skipped) + "\n")
            f.write("Error message:\n")
            f.write(openmp_error.strip() + "\n")
            f.write("\nSuggested fixes:\n")
            f.write("- Use conda environment (recommended)\n")
            f.write("- Install libomp (Mac)\n")
            f.write("- Use non-OpenMP builds\n")


def write_run_metadata(outdir, args, cfg, organism_args, params, run_selection):
    meta = {
        "timestamp": datetime.now().isoformat(),
        "run_dir": str(Path(outdir).resolve()),
        "output_root": str(get_output_root()),
        "counts": cfg.get("counts_path"),
        "metadata": cfg.get("metadata_path"),
        "organism": organism_args.get("organism"),
        "kegg_org": organism_args.get("kegg_org"),
        "analysis_name": params.get("analysis_name"),
        "run_selection": run_selection,
        "cli": sys.argv,
    }
    path = Path(outdir) / "run_metadata.yaml"
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(meta, f, allow_unicode=True, default_flow_style=False)


def filter_openmp_steps(run_list):
    if OPENMP_AVAILABLE:
        return run_list
    skipped = [n for n in run_list if n in OPENMP_STEPS]
    if skipped:
        LOG.warning("OpenMP runtime failure detected (OMP Error #179)")
        LOG.warning("Skipping steps: %s", ", ".join(f"{n:02d}" for n in skipped))
    return [n for n in run_list if n not in OPENMP_STEPS]


RSCRIPT = None
R_DEPS_INSTALLED = False
R_LIBS_USER = BASE_DIR / ".r_libs"


def get_output_root():
    configured = os.environ.get("PIPELINE_OUTPUT_DIR", "").strip()
    if configured:
        return Path(configured).expanduser().resolve()
    return DEFAULT_OUTPUT_ROOT.resolve()


def create_run_dir(base_root: Path) -> Path:
    base_root.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("run_%Y-%m-%d_%H-%M-%S")
    candidate = base_root / ts
    suffix = 1
    while candidate.exists():
        candidate = base_root / f"{ts}_{suffix}"
        suffix += 1
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate


def init_rscript():
    global RSCRIPT
    if RSCRIPT is None:
        RSCRIPT = find_rscript()
        if RSCRIPT is None:
            LOG.error("Rscript が見つかりません。R をインストールしてください。")
            sys.exit(1)


def ensure_r_dependencies():
    """R依存関係をインストール/確認"""
    global R_DEPS_INSTALLED
    if R_DEPS_INSTALLED:
        return
    init_rscript()
    install_script = R_DIR / "install_missing_packages.R"
    if not install_script.exists():
        LOG.error("依存関係インストールスクリプトが見つかりません: %s", install_script)
        sys.exit(1)
    R_LIBS_USER.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["R_LIBS_USER"] = str(R_LIBS_USER)
    cmd = [RSCRIPT, "--vanilla", str(install_script)]
    LOG.info("R依存関係を確認/インストールします: %s", " ".join(cmd))
    result = subprocess.run(cmd, cwd=str(BASE_DIR), capture_output=True, text=True, env=env)
    if result.stdout:
        for line in result.stdout.strip().splitlines():
            LOG.info("%s", line)
    if result.stderr:
        for line in result.stderr.strip().splitlines():
            LOG.error("%s", line)
    if result.returncode != 0:
        LOG.error("R依存関係のインストールに失敗しました (code=%s)", result.returncode)
        sys.exit(1)
    R_DEPS_INSTALLED = True


def load_config(path):
    path = Path(path)
    if path.exists():
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f) or {}
    return {}


def save_config(cfg, path):
    path = Path(path)
    with open(path, "w", encoding="utf-8") as f:
        yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True)


def infer_choice_from_orgdb(orgdb):
    """org.*.eg.db から ORG_OPTIONS のキーを推定"""
    if not orgdb:
        return None
    orgdb_lower = orgdb.lower()
    if "org.hs.eg.db" in orgdb_lower:
        return "human"
    if "org.mm.eg.db" in orgdb_lower:
        return "mouse"
    if "org.dr.eg.db" in orgdb_lower:
        return "zebrafish"
    return None


def run_r_script(script_name, args_dict, dry_run=False, precheck=True):
    """R スクリプトを実行する。引数は --key value 形式で渡す"""
    init_rscript()
    if not dry_run:
        ensure_r_dependencies()
    script_path = R_DIR / script_name
    if not script_path.exists():
        LOG.error("スクリプトが見つかりません: %s", script_path)
        return False

    cmd = [RSCRIPT, "--vanilla", str(script_path)]
    for k, v in args_dict.items():
        cmd.extend([f"--{k}", str(v)])

    if precheck and script_name != "00_align_check.R":
        LOG.info("サンプル整合性チェックを実行します: 00_align_check.R")
        align_args = {
            "counts": args_dict.get("counts"),
            "metadata": args_dict.get("metadata"),
            "outdir": args_dict.get("outdir"),
        }
        ok = run_r_script("00_align_check.R", align_args, dry_run=False, precheck=False)
        if not ok:
            LOG.error("サンプル整合性チェックに失敗しました。")
            return False

    LOG.info("============================================================")
    LOG.info("実行: %s", script_name)
    LOG.info("コマンド: %s", " ".join(cmd))
    LOG.info("============================================================")

    if dry_run and script_name != "00_align_check.R":
        LOG.info("dry-run: 実行をスキップしました。")
        return True

    env = os.environ.copy()
    env["R_LIBS_USER"] = str(R_LIBS_USER)
    if "outdir" in args_dict and args_dict["outdir"]:
        env["PIPELINE_OUTPUT_DIR"] = str(args_dict["outdir"])
    else:
        env["PIPELINE_OUTPUT_DIR"] = str(get_output_root())
    result = subprocess.run(cmd, cwd=str(BASE_DIR), capture_output=True, text=True, env=env)
    if result.stdout:
        for line in result.stdout.strip().splitlines():
            LOG.info("%s", line)
    if result.stderr:
        for line in result.stderr.strip().splitlines():
            LOG.error("%s", line)
    if result.returncode != 0:
        LOG.error("%s がエラーで終了しました (code=%s)", script_name, result.returncode)
        return False
    LOG.info("完了: %s", script_name)
    return True


# ─────────────────────────────────────────────
# 入力ヘルパー
# ─────────────────────────────────────────────

def input_path(prompt, must_exist=True):
    while True:
        p = input(prompt).strip().strip('"').strip("'")
        if not p:
            return None
        path = Path(p)
        if must_exist and not path.exists():
            LOG.warning("ファイルが見つかりません: %s", path)
            continue
        return str(path.resolve())


def input_choice(prompt, choices):
    cstr = "/".join(choices)
    while True:
        v = input(f"{prompt} [{cstr}]: ").strip().lower()
        if v in [c.lower() for c in choices]:
            return v
        LOG.warning("次のいずれかを入力してください: %s", cstr)


def input_int(prompt, default=None):
    while True:
        v = input(prompt).strip()
        if not v and default is not None:
            return default
        try:
            return int(v)
        except ValueError:
            LOG.warning("整数を入力してください。")


def input_float(prompt, default=None):
    while True:
        v = input(prompt).strip()
        if not v and default is not None:
            return default
        try:
            return float(v)
        except ValueError:
            LOG.warning("数値を入力してください。")


def input_text(prompt, default=""):
    v = input(prompt).strip()
    return v if v else default


def value_choice(prompt, choices, default, non_interactive):
    if non_interactive:
        return default if default in choices else choices[0]
    return input_choice(prompt, choices)


def value_int(prompt, default, non_interactive):
    if non_interactive:
        return default
    return input_int(prompt, default)


def value_float(prompt, default, non_interactive):
    if non_interactive:
        return default
    return input_float(prompt, default)


def value_text(prompt, default, non_interactive):
    if non_interactive:
        return default
    return input_text(prompt, default)


def value_path(prompt, default, non_interactive, must_exist=False):
    if non_interactive:
        return default
    return input_path(prompt, must_exist=must_exist)


def get_param(cfg, args, key, default):
    val = getattr(args, key, None)
    if val is not None:
        return val
    if key in cfg and cfg[key] not in (None, ""):
        return cfg[key]
    return default


def resolve_outdir(base_outdir):
    result_root = get_output_root()
    requested = Path(base_outdir).resolve() if base_outdir else result_root
    run_dir = create_run_dir(result_root)
    return run_dir, requested


def get_step_output_dir(base_outdir, step_dir, analysis_name=None):
    outdir = Path(base_outdir) / step_dir
    if analysis_name and analysis_name != step_dir:
        outdir = outdir / analysis_name
    return outdir


def parse_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    v = str(value).strip().lower()
    if v in {"true", "1", "yes", "y"}:
        return True
    if v in {"false", "0", "no", "n"}:
        return False
    return default


# ─────────────────────────────────────────────
# 共通引数の取得 / 保存
# ─────────────────────────────────────────────

def get_common_args(cfg, args, interactive=True, config_path=CONFIG_FILE):
    """count matrix と metadata のパスを取得・保持（CLI > config > defaults）"""
    counts = args.counts or cfg.get("counts_path")
    meta = args.metadata or cfg.get("metadata_path")

    if interactive:
        if counts and Path(counts).exists():
            LOG.info("Count matrix (保存済み): %s", counts)
            ch = input_choice("  このパスを使用しますか？", ["y", "n"])
            if ch == "n":
                counts = None

        if not counts:
            counts = input_path("  Count matrix (.rds) のパス: ")
            if not counts:
                LOG.error("Count matrix は必須です。")
                return None

        if meta and Path(meta).exists():
            LOG.info("Metadata (保存済み): %s", meta)
            ch = input_choice("  このパスを使用しますか？", ["y", "n"])
            if ch == "n":
                meta = None

        if not meta:
            meta = input_path("  Metadata (.xlsx / .csv) のパス: ")
            if not meta:
                LOG.error("Metadata は必須です。")
                return None
    else:
        if not counts or not meta:
            LOG.error("--non-interactive では --counts と --metadata が必須です。")
            return None

    counts = str(Path(counts).resolve())
    meta = str(Path(meta).resolve())

    if not Path(counts).exists():
        LOG.error("Count matrix が見つかりません: %s", counts)
        return None
    if not Path(meta).exists():
        LOG.error("Metadata が見つかりません: %s", meta)
        return None

    cfg["counts_path"] = counts
    cfg["metadata_path"] = meta
    if interactive:
        save_config(cfg, config_path)
    return {"counts": counts, "metadata": meta}


def get_organism_args(cfg, args, interactive=True, config_path=CONFIG_FILE):
    """生物種選択を一元管理し、OrgDb/KEGGコードを返す"""
    global ORG_ARGS
    if args.organism:
        choice = args.organism
    else:
        saved_orgdb = cfg.get("organism")
        if saved_orgdb in ORG_OPTIONS:
            inferred_choice = saved_orgdb
        else:
            inferred_choice = infer_choice_from_orgdb(saved_orgdb)

        if interactive and saved_orgdb:
            LOG.info("生物種 (保存済み): %s", saved_orgdb)
            ch = input_choice("  この生物種を使用しますか？", ["y", "n"])
            if ch == "n":
                inferred_choice = None

        if inferred_choice is None:
            if interactive:
                LOG.info("生物種を選択してください (human / mouse / zebrafish)")
                choice = input_choice("  Organism", list(ORG_OPTIONS.keys()))
            else:
                choice = "human"
        else:
            choice = inferred_choice

    orgdb = ORG_OPTIONS[choice]["orgdb"]
    kegg = ORG_OPTIONS[choice]["kegg"]
    cfg["organism"] = orgdb
    if interactive:
        save_config(cfg, config_path)

    ORG_ARGS = {"organism": orgdb, "kegg_org": kegg, "organism_label": choice}
    LOG.info("選択: %s -> OrgDb=%s, KEGG=%s", choice, orgdb, kegg)
    return ORG_ARGS


# ─────────────────────────────────────────────
# 各解析のパラメータ取得 & 実行
# ─────────────────────────────────────────────

def run_01_de(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 1. 差次発現解析 (Differential Expression) ---")
    while True:
        tool = value_choice("  ツール", ["deseq2", "edger", "limma"], "deseq2", non_interactive)
        LOG.info("  条件 → サンプルをグループ分けする列")
        LOG.info("  参照 → 比較の基準となるグループ")
        LOG.info("  処理 → 比較したいグループ")
        condition_col = value_text("  条件カラム名 (例: condition): ", "condition", non_interactive)
        ref_level = value_text("  参照レベル (例: control): ", "control", non_interactive)
        treat_level = value_text("  処理レベル (例: treatment): ", "treatment", non_interactive)
        fdr_cutoff = value_float("  FDR閾値 [0.05]: ", params["fdr_cutoff"], non_interactive)
        lfc_cutoff = value_float("  log2FC閾値 [1.0]: ", params["lfc_cutoff"], non_interactive)

        design_mode = value_choice("  デザインモード [simple/additive/interaction]: ",
                                   ["simple", "additive", "interaction"], params["design_mode"], non_interactive)
        factors = value_text("  因子 (カンマ区切り。空=condition): ", params["factors"], non_interactive)
        factors = factors.strip("'\"")  # strip accidental surrounding quotes
        contrast_mode = value_choice("  コントラストモード [pairwise/auto/manual]: ",
                                     ["pairwise", "auto", "manual"], params["contrast_mode"], non_interactive)
        contrast = value_text("  manual用コントラスト名 (空=未指定): ", params["contrast"], non_interactive)
        test = value_choice("  DESeq2 検定 [Wald/LRT]: ",
                            ["Wald", "LRT"], params["de_test"], non_interactive)
        reduced_design = value_text("  LRT の reduced design (空=未指定): ",
                                    params["reduced_design"], non_interactive)
        subset_col = value_text("  subset metadata column (空=なし): ",
                                params["subset_col"], non_interactive)
        subset_val = value_text("  subset values (カンマ区切り、空=なし): ",
                                params["subset_val"], non_interactive)
        analysis_name = value_text("  解析名 (出力フォルダ名) [01_DE]: ",
                                   params["analysis_name"], non_interactive)
        explain_results = value_choice("  結果の解釈レポートを作成しますか？ [true/false]: ",
                                       ["true", "false"],
                                       "true" if params["explain_results"] else "false",
                                       non_interactive)

        if not non_interactive:
            LOG.info("  実行設定: tool=%s, condition=%s, %s vs %s, design=%s, contrast=%s",
                     tool, condition_col, treat_level, ref_level, design_mode, contrast_mode)
            ch = value_choice("  この設定で実行しますか？ [y/n]: ", ["y", "n"], "y", False)
            if ch == "n":
                LOG.info("  再入力します。")
                continue

        args = {**common,
                "tool": tool,
                "condition_col": condition_col,
                "ref_level": ref_level,
                "treat_level": treat_level,
                "fdr_cutoff": fdr_cutoff,
                "lfc_cutoff": lfc_cutoff,
                "design_mode": design_mode,
                "factors": factors,
                "contrast_mode": contrast_mode,
                "contrast": contrast,
                "test": test,
                "reduced_design": reduced_design,
                "subset_col": subset_col,
                "subset_val": subset_val,
                "analysis_name": analysis_name,
                "explain_results": explain_results}
        return run_r_script("01_differential_expression.R", args, dry_run=dry_run)


def run_02_clustering(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 2. サンプルクラスタリング ---")
    method = value_choice("  手法", ["hierarchical", "kmeans", "both"], "both", non_interactive)
    k = value_int("  クラスタ数 (k-means用) [3]: ", 3, non_interactive)
    dist_method = value_choice("  距離指標", ["euclidean", "pearson", "spearman"], "euclidean", non_interactive)
    top_n = value_int("  使用する上位変動遺伝子数 [1000]: ", 1000, non_interactive)
    target_genes = value_path("  対象遺伝子CSV (空=top_nを使用): ",
                              params["target_genes"], non_interactive, must_exist=False)
    trend_x_col = value_text("  Trend plot の x 軸 metadata column (空=スキップ): ",
                             params["trend_x_col"], non_interactive)
    label_mode = value_choice("  サンプルラベル表示 [sampleID_only/symbol_only/both]: ",
                              ["sampleID_only", "symbol_only", "both"],
                              params["label_mode"], non_interactive)
    symbol_col = value_text("  ラベルに使う metadata column (必要時のみ): ",
                            params["symbol_col"], non_interactive)

    args = {**common, "method": method, "k": k,
            "dist_method": dist_method, "top_n": top_n,
            "target_genes": target_genes, "trend_x_col": trend_x_col,
            "label_mode": label_mode, "symbol_col": symbol_col}
    return run_r_script("02_sample_clustering.R", args, dry_run=dry_run)


def run_03_dimreduc(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 3. 次元削減解析 (Diagnostic) ---")
    methods = value_text("  手法 (カンマ区切り: pca,tsne,umap) [pca,umap]: ", "pca,umap", non_interactive)
    color_by = value_text("  色分けカラム (metadata) [condition]: ", "condition", non_interactive)
    top_n = value_int("  使用する上位変動遺伝子数 [2000]: ", 2000, non_interactive)
    perplexity = value_int("  t-SNE perplexity [30]: ", 30, non_interactive)
    n_neighbors = value_int("  UMAP n_neighbors [15]: ", 15, non_interactive)

    # New diagnostic options
    transform = value_choice("  データ変換 [log2/vst/rlog]: ",
                             ["log2", "vst", "rlog"], "log2", non_interactive)
    split_by = value_text("  因子で分割解析 (空=なし, 例: dpf): ", "", non_interactive)
    color_by_all = value_choice("  全メタデータカラムでPCA色分け [true/false]: ",
                                ["true", "false"], "false", non_interactive)
    feature_method = value_choice("  特徴量選択 [variance/mad]: ",
                                  ["variance", "mad"], "variance", non_interactive)
    label_samples = value_choice("  サンプルラベル表示 [true/false]: ",
                                 ["true", "false"], "true", non_interactive)
    label_mode = value_choice("  サンプルラベル形式 [sampleID_only/symbol_only/both]: ",
                              ["sampleID_only", "symbol_only", "both"],
                              params["label_mode"], non_interactive)
    symbol_col = value_text("  ラベルに使う metadata column (必要時のみ): ",
                            params["symbol_col"], non_interactive)
    perplexity_list = value_text("  t-SNE perplexityリスト (カンマ区切り, 空=単一値): ", "", non_interactive)
    n_neighbors_list = value_text("  UMAP n_neighborsリスト (カンマ区切り, 空=単一値): ", "", non_interactive)

    args = {**common, "methods": methods, "color_by": color_by,
            "top_n": top_n, "perplexity": perplexity, "n_neighbors": n_neighbors,
            "transform": transform, "split_by": split_by,
            "color_by_all": color_by_all, "feature_method": feature_method,
            "label_samples": label_samples, "label_mode": label_mode,
            "symbol_col": symbol_col, "perplexity_list": perplexity_list,
            "n_neighbors_list": n_neighbors_list}
    return run_r_script("03_dimension_reduction.R", args, dry_run=dry_run)


def run_04_normalization(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 4. 発現量の正規化 ---")
    method = value_choice("  正規化手法", ["deseq2", "tmm", "log2", "all"], "all", non_interactive)
    label_mode = value_choice("  サンプルラベル形式 [sampleID_only/symbol_only/both]: ",
                              ["sampleID_only", "symbol_only", "both"],
                              params["label_mode"], non_interactive)
    symbol_col = value_text("  ラベルに使う metadata column (必要時のみ): ",
                            params["symbol_col"], non_interactive)
    args = {**common, "method": method, "label_mode": label_mode, "symbol_col": symbol_col}
    return run_r_script("04_normalization.R", args, dry_run=dry_run)


def run_05_expression_pattern(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 5. 発現パターン解析 (WGCNA) ---")
    power_range = value_text("  ソフト閾値の探索範囲 [1:20]: ", "1:20", non_interactive)
    min_module_size = value_int("  最小モジュールサイズ [30]: ", 30, non_interactive)
    merge_cut = value_float("  モジュール統合閾値 [0.25]: ", 0.25, non_interactive)
    trait_col = value_text("  形質カラム (trait association用) [condition]: ", "condition", non_interactive)

    args = {**common, "power_range": power_range,
            "min_module_size": min_module_size,
            "merge_cut": merge_cut, "trait_col": trait_col}
    return run_r_script("05_expression_pattern.R", args, dry_run=dry_run)


def run_06_enrichment(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 6. 機能解析 (GO/KEGG/Reactome) ---")
    de_default = str(get_step_output_dir(common["outdir"], "01_DE", params["analysis_name"]) / "de_results.csv")
    de_result = value_path("  DE結果ファイル (CSVパス、空欄で自動検索): ",
                           de_default, non_interactive, must_exist=False)
    if not de_result:
        de_result = de_default
    pval_cutoff = value_float("  p値閾値 [0.05]: ", params["enrich_pval_cutoff"], non_interactive)

    if not ORG_ARGS:
        LOG.warning("生物種情報が初期化されていません。デフォルト(human)を使用します。")
        orgdb = ORG_OPTIONS["human"]["orgdb"]
        kegg_org = ORG_OPTIONS["human"]["kegg"]
    else:
        orgdb = ORG_ARGS["organism"]
        kegg_org = ORG_ARGS["kegg_org"]

    args = {**common, "de_result": de_result, "organism": orgdb,
            "kegg_org": kegg_org, "pval_cutoff": pval_cutoff}
    return run_r_script("06_functional_enrichment.R", args, dry_run=dry_run)


def run_07_gsea(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 7. GSEA ---")
    de_default = str(get_step_output_dir(common["outdir"], "01_DE", params["analysis_name"]) / "de_results.csv")
    de_result = value_path("  DE結果ファイル (CSVパス、空欄で自動検索): ",
                           de_default, non_interactive, must_exist=False)
    if not de_result:
        de_result = de_default
    geneset = value_choice("  Gene set", ["GO", "KEGG", "Reactome", "all"], params["gsea_geneset"], non_interactive)
    ranking_metric = value_choice("  Ranking metric", ["log2FoldChange", "stat", "custom"],
                                  params["gsea_ranking_metric"], non_interactive)
    min_gs = value_int("  最小遺伝子セットサイズ [15]: ", params["gsea_min_gs"], non_interactive)
    max_gs = value_int("  最大遺伝子セットサイズ [500]: ", params["gsea_max_gs"], non_interactive)

    if not ORG_ARGS:
        LOG.warning("生物種情報が初期化されていません。デフォルト(human)を使用します。")
        orgdb = ORG_OPTIONS["human"]["orgdb"]
        kegg_org = ORG_OPTIONS["human"]["kegg"]
    else:
        orgdb = ORG_ARGS["organism"]
        kegg_org = ORG_ARGS["kegg_org"]

    args = {**common, "de_result": de_result, "organism": orgdb,
            "kegg_org": kegg_org, "geneset": geneset,
            "ranking_metric": ranking_metric, "min_gs": min_gs, "max_gs": max_gs}
    return run_r_script("07_gsea.R", args, dry_run=dry_run)


def run_08_batch(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 8. バッチ効果補正 ---")
    batch_col = value_text("  Metadataのバッチカラム名 [batch]: ", "batch", non_interactive)
    method = value_choice("  補正手法", ["combat", "limma", "both"], "both", non_interactive)
    condition_col = value_text("  保持する生物学的条件カラム [condition]: ", "condition", non_interactive)
    label_mode = value_choice("  サンプルラベル形式 [sampleID_only/symbol_only/both]: ",
                              ["sampleID_only", "symbol_only", "both"],
                              params["label_mode"], non_interactive)
    symbol_col = value_text("  ラベルに使う metadata column (必要時のみ): ",
                            params["symbol_col"], non_interactive)

    args = {**common, "batch_col": batch_col, "method": method,
            "condition_col": condition_col, "label_mode": label_mode,
            "symbol_col": symbol_col}
    return run_r_script("08_batch_correction.R", args, dry_run=dry_run)


def run_09_regression(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 9. 回帰モデル解析 ---")
    formula_str = value_text(
        "  モデル式 (例: ~ treatment + sex + age) [~ condition]: ", "~ condition", non_interactive)
    gene_list = value_text("  対象遺伝子リスト (カンマ区切り、空=全遺伝子) []: ", "", non_interactive)
    method = value_choice("  手法", ["lm", "glm", "limma"], "limma", non_interactive)

    args = {**common, "formula_str": formula_str,
            "gene_list": gene_list, "method": method}
    return run_r_script("09_regression.R", args, dry_run=dry_run)


def run_10_qc(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 10. サンプル品質評価 ---")
    color_by = value_text("  色分けカラム [condition]: ", "condition", non_interactive)
    args = {**common, "color_by": color_by}
    return run_r_script("10_sample_qc.R", args, dry_run=dry_run)


def run_11_vis(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 11. 可視化 ---")
    de_default = str(get_step_output_dir(common["outdir"], "01_DE", params["analysis_name"]) / "de_results.csv")
    de_result = value_path("  DE結果CSVパス (空欄で自動検索): ",
                           de_default, non_interactive, must_exist=False)
    if not de_result:
        de_result = de_default
    plots = value_text(
        "  作成する図 (カンマ区切り: heatmap,volcano,ma,pca) [heatmap,volcano,ma,pca]: ",
        "heatmap,volcano,ma,pca", non_interactive)
    top_n = value_int("  ヒートマップ上位遺伝子数 [50]: ", 50, non_interactive)
    color_by = value_text("  PCA色分けカラム [condition]: ", "condition", non_interactive)
    label_mode = value_choice("  サンプルラベル形式 [sampleID_only/symbol_only/both]: ",
                              ["sampleID_only", "symbol_only", "both"],
                              params["label_mode"], non_interactive)
    symbol_col = value_text("  ラベルに使う metadata column (必要時のみ): ",
                            params["symbol_col"], non_interactive)

    args = {**common, "de_result": de_result, "plots": plots,
            "top_n": top_n, "color_by": color_by,
            "label_mode": label_mode, "symbol_col": symbol_col}
    return run_r_script("11_visualization.R", args, dry_run=dry_run)


def run_12_phenotype(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 12. 疾患・表現型関連解析 ---")
    phenotype_col = value_text("  表現型カラム名 [condition]: ", "condition", non_interactive)
    method = value_choice("  解析手法", ["correlation", "association", "signature", "all"], "all", non_interactive)
    top_n = value_int("  上位遺伝子数 [100]: ", 100, non_interactive)

    args = {**common, "phenotype_col": phenotype_col,
            "method": method, "top_n": top_n}
    return run_r_script("12_phenotype_association.R", args, dry_run=dry_run)


def run_13_ml(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 13. 機械学習解析 ---")
    target_col = value_text("  予測対象カラム [condition]: ", "condition", non_interactive)
    method = value_choice("  手法", ["rf", "svm", "both"], "both", non_interactive)
    top_n = value_int("  特徴量とする上位変動遺伝子数 [500]: ", 500, non_interactive)
    cv_folds = value_int("  交差検証 fold数 [5]: ", 5, non_interactive)
    test_ratio = value_float("  テストデータ比率 [0.3]: ", 0.3, non_interactive)

    args = {**common, "target_col": target_col, "method": method,
            "top_n": top_n, "cv_folds": cv_folds, "test_ratio": test_ratio}
    return run_r_script("13_machine_learning.R", args, dry_run=dry_run)


def run_14_network(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 14. ネットワーク解析 ---")
    method = value_choice("  手法", ["coexpression", "wgcna", "both"], "both", non_interactive)
    top_n = value_int("  使用する上位変動遺伝子数 [500]: ", 500, non_interactive)
    cor_threshold = value_float("  相関閾値 [0.8]: ", 0.8, non_interactive)
    hub_top = value_int("  ハブ遺伝子上位数 [20]: ", 20, non_interactive)
    detect_modules = value_choice("  WGCNA モジュール検出 [true/false]: ",
                                  ["true", "false"],
                                  "true" if params["detect_modules"] else "false",
                                  non_interactive)
    trait_cols = value_text("  trait columns (カンマ区切り、空=なし): ",
                            params["trait_cols"], non_interactive)

    args = {**common, "method": method, "top_n": top_n,
            "cor_threshold": cor_threshold, "hub_top": hub_top,
            "detect_modules": detect_modules, "trait_cols": trait_cols}
    return run_r_script("14_network_analysis.R", args, dry_run=dry_run)


def run_15_sc(common, params, non_interactive=False, dry_run=False):
    LOG.info("--- 15. 疑似bulk / single-cell 統合 ---")
    sc_data = value_path("  Single-cell データパス (.rds): ", None, non_interactive, must_exist=False)
    method = value_choice("  手法", ["umap", "clustering", "marker", "all"], "all", non_interactive)
    resolution = value_float("  クラスタリング resolution [0.5]: ", 0.5, non_interactive)
    n_neighbors = value_int("  UMAP n_neighbors [30]: ", 30, non_interactive)

    args = {**common, "method": method,
            "resolution": resolution, "n_neighbors": n_neighbors}
    if sc_data:
        args["sc_data"] = sc_data
    return run_r_script("15_single_cell_integration.R", args, dry_run=dry_run)


# ─────────────────────────────────────────────
# バッチ実行
# ─────────────────────────────────────────────

def run_batch(common, params, dry_run=False):
    """複数解析を連続実行"""
    LOG.info("--- バッチ実行モード ---")
    LOG.info("実行する解析番号をカンマ区切りで入力 (例: 10,4,1,3)")
    LOG.info("'all' で全解析実行")
    sel = input_text("  選択: ", "")
    if not sel:
        return

    if sel.lower() == "all":
        nums = list(range(1, 16))
    else:
        try:
            nums = [int(x.strip()) for x in sel.split(",")]
        except ValueError:
            LOG.error("無効な入力です。")
            return

    nums = filter_openmp_steps(nums)
    if not nums:
        LOG.warning("実行可能な解析がありません。")
        return

    dispatch = {
        1: run_01_de, 2: run_02_clustering, 3: run_03_dimreduc,
        4: run_04_normalization, 5: run_05_expression_pattern,
        6: run_06_enrichment, 7: run_07_gsea, 8: run_08_batch,
        9: run_09_regression, 10: run_10_qc, 11: run_11_vis,
        12: run_12_phenotype, 13: run_13_ml, 14: run_14_network,
        15: run_15_sc,
    }

    for n in nums:
        if n in dispatch:
            ok = dispatch[n](common, params, dry_run=dry_run)
            if not ok:
                LOG.error("step %s failed. Exiting batch.", n)
                return
        else:
            LOG.warning("解析番号 %s は無効です。スキップします。", n)


# ─────────────────────────────────────────────
# メインメニュー
# ─────────────────────────────────────────────

MENU = """
╔══════════════════════════════════════════════════════════╗
║           RNA-seq 統合解析パイプライン                    ║
╠══════════════════════════════════════════════════════════╣
║  [ 統計解析 ]                                            ║
║   1.  差次発現解析 (DESeq2 / edgeR / limma)              ║
║   9.  回帰モデル解析                                     ║
║                                                          ║
║  [ サンプル構造解析 ]                                    ║
║   2.  サンプルクラスタリング                             ║
║   3.  次元削減解析 (PCA / t-SNE / UMAP)                  ║
║   4.  発現量の正規化                                     ║
║  10.  サンプル品質評価 (QC)                              ║
║   8.  バッチ効果補正                                     ║
║                                                          ║
║  [ 遺伝子機能解析 ]                                      ║
║   5.  発現パターン解析 (WGCNA)                           ║
║   6.  機能解析 (GO / KEGG / Reactome)                    ║
║   7.  GSEA                                               ║
║                                                          ║
║  [ ネットワーク解析 ]                                    ║
║  14.  ネットワーク解析                                   ║
║                                                          ║
║  [ 予測・機械学習 ]                                      ║
║  13.  機械学習解析 (Random Forest / SVM)                 ║
║                                                          ║
║  [ 関連解析 ]                                            ║
║  12.  疾患・表現型関連解析                               ║
║  15.  疑似bulk / single-cell 統合                        ║
║                                                          ║
║  [ その他 ]                                              ║
║  11.  可視化                                             ║
║   B.  バッチ実行 (複数解析を連続実行)                    ║
║   S.  設定確認・変更                                     ║
║   Q.  終了                                               ║
╚══════════════════════════════════════════════════════════╝
"""


def main():
    parser = argparse.ArgumentParser(description="RNA-seq 統合解析パイプライン")
    parser.add_argument("--counts", type=str, help="Count matrix (.rds) path")
    parser.add_argument("--metadata", type=str, help="Metadata (.xlsx/.csv) path")
    parser.add_argument("--organism", type=str, choices=list(ORG_OPTIONS.keys()),
                        help="Organism: human, mouse, zebrafish")
    parser.add_argument("--config", type=str, default=str(CONFIG_FILE),
                        help="Config YAML path")
    parser.add_argument("--non-interactive", action="store_true",
                        help="Disable prompts and require CLI values")
    parser.add_argument("--run", type=str, help="Run step: all, 1, 6, or 7")
    parser.add_argument("--outdir", type=str, help="Output directory (retained for compatibility; results are written to ~/Output or PIPELINE_OUTPUT_DIR)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate and print commands without executing")
    parser.add_argument("--describe", action="store_true",
                        help="Describe analysis steps in English and Japanese")
    parser.add_argument("--list-steps", action="store_true",
                        help="List available analysis step IDs")
    parser.add_argument("--fdr-cutoff", type=float, help="FDR cutoff for DE")
    parser.add_argument("--lfc-cutoff", type=float, help="log2FC cutoff for DE")
    parser.add_argument("--enrich-pval-cutoff", type=float, help="p-value cutoff for enrichment")
    parser.add_argument("--gsea-geneset", type=str, help="GSEA geneset: GO/KEGG/Reactome/all")
    parser.add_argument("--gsea-ranking-metric", type=str, choices=["log2FoldChange", "stat", "custom"],
                        help="GSEA ranking metric")
    parser.add_argument("--gsea-min-gs", type=int, help="GSEA min gene set size")
    parser.add_argument("--gsea-max-gs", type=int, help="GSEA max gene set size")
    parser.add_argument("--analysis-name", type=str, help="DE analysis output folder name")
    parser.add_argument("--design-mode", type=str, choices=["simple", "additive", "interaction"],
                        help="Design mode for DE")
    parser.add_argument("--factors", type=str, help="Comma-separated metadata columns for design")
    parser.add_argument("--contrast-mode", type=str, choices=["pairwise", "auto", "manual"],
                        help="Contrast mode for DE")
    parser.add_argument("--contrast", type=str, help="Manual contrast name (advanced)")
    parser.add_argument("--de-test", type=str, choices=["Wald", "LRT"], help="DESeq2 test")
    parser.add_argument("--reduced-design", type=str, help="Reduced design for DESeq2 LRT")
    parser.add_argument("--subset-col", type=str, help="Metadata column used for DE subsetting")
    parser.add_argument("--subset-val", type=str, help="Comma-separated metadata values used for DE subsetting")
    parser.add_argument("--target-genes", type=str, help="Target gene CSV for clustering")
    parser.add_argument("--trend-x-col", type=str, help="Metadata column used on the clustering trend-plot x-axis")
    parser.add_argument("--label-mode", type=str, choices=["sampleID_only", "symbol_only", "both"],
                        help="Sample label format for plots")
    parser.add_argument("--symbol-col", type=str, help="Metadata column used as the sample symbol label")
    parser.add_argument("--detect-modules", type=str, help="TRUE/FALSE to run WGCNA module detection")
    parser.add_argument("--trait-cols", type=str, help="Comma-separated metadata columns for module-trait correlation")
    parser.add_argument("--explain-results", type=str, help="TRUE/FALSE to generate interpretation report")
    parser.add_argument("--example", action="store_true",
                        help="Run pipeline with built-in example dataset (no input files needed)")
    args = parser.parse_args()

    setup_logging()

    if args.list_steps:
        for k in sorted(STEP_NAMES.keys()):
            LOG.info("%s: %s", k, STEP_NAMES[k])
        return
    if args.describe:
        for step_id, desc in STEP_DESCRIPTIONS.items():
            LOG.info("[Step %s: %s]", step_id, desc["title"])
            LOG.info("EN: %s", desc["en"])
            LOG.info("JP: %s", desc["jp"])
        return

    # ---- Example mode ----
    if args.example:
        LOG.info("============================================================")
        LOG.info("EXAMPLE MODE: Running pipeline with synthetic dataset")
        LOG.info("============================================================")
        example_dir = BASE_DIR / "example"
        example_counts = example_dir / "counts.rds"
        example_meta = example_dir / "metadata.csv"

        # Generate example data if not present
        if not example_counts.exists() or not example_meta.exists():
            gen_script = example_dir / "generate_example.R"
            if not gen_script.exists():
                LOG.error("Example generator not found: %s", gen_script)
                sys.exit(1)
            init_rscript()
            LOG.info("Generating example dataset...")
            result = subprocess.run(
                [RSCRIPT, "--vanilla", str(gen_script)],
                cwd=str(BASE_DIR), capture_output=True, text=True
            )
            if result.returncode != 0:
                LOG.error("Example generation failed: %s", result.stderr)
                sys.exit(1)
            for line in (result.stderr or "").strip().splitlines():
                LOG.info("%s", line)

        # Override args for example mode
        args.counts = str(example_counts)
        args.metadata = str(example_meta)
        args.organism = args.organism or "human"
        args.non_interactive = True
        if not args.run:
            args.run = "1"  # default: run DE only for quick demo
        if not args.outdir:
            args.outdir = str(get_output_root())
        LOG.info("Example counts: %s", args.counts)
        LOG.info("Example metadata: %s", args.metadata)
        LOG.info("Example output: %s", args.outdir)

    LOG.info("============================================================")
    LOG.info("RNA-seq 統合解析パイプライン v1.0")
    LOG.info("============================================================")

    cfg = load_config(args.config)
    outdir_value = args.outdir or cfg.get("outdir")
    outdir, outdir_base = resolve_outdir(outdir_value)
    setup_logging(outdir)
    if outdir_base != outdir:
        LOG.info("出力先は固定です。指定値 %s ではなく %s を使用します。", outdir_base, outdir)

    ensure_dirs(outdir)

    ensure_r_dependencies()

    global OPENMP_AVAILABLE, OPENMP_ERROR
    try:
        OPENMP_AVAILABLE, OPENMP_ERROR = check_openmp_runtime()
    except Exception as exc:
        LOG.error("OpenMP pre-flight check failed: %s", exc)
        sys.exit(1)
    if not OPENMP_AVAILABLE:
        LOG.warning("OpenMP runtime failure detected (OMP Error #179)")
        LOG.warning("Skipping steps: 05, 06, 07, 13, 14, 15")
    write_environment_report(outdir, OPENMP_AVAILABLE, OPENMP_ERROR)

    # 共通引数
    common = get_common_args(cfg, args, interactive=not args.non_interactive, config_path=args.config)
    if common is None:
        sys.exit(1)
    common["outdir"] = str(outdir)
    if not args.non_interactive:
        cfg["outdir"] = str(outdir)
        save_config(cfg, args.config)
    params = {
        "fdr_cutoff": get_param(cfg, args, "fdr_cutoff", 0.05),
        "lfc_cutoff": get_param(cfg, args, "lfc_cutoff", 1.0),
        "enrich_pval_cutoff": get_param(cfg, args, "enrich_pval_cutoff", 0.05),
        "gsea_geneset": get_param(cfg, args, "gsea_geneset", "all"),
        "gsea_ranking_metric": get_param(cfg, args, "gsea_ranking_metric", "log2FoldChange"),
        "gsea_min_gs": get_param(cfg, args, "gsea_min_gs", 15),
        "gsea_max_gs": get_param(cfg, args, "gsea_max_gs", 500),
        "analysis_name": get_param(cfg, args, "analysis_name", "01_DE"),
        "design_mode": get_param(cfg, args, "design_mode", "simple"),
        "factors": get_param(cfg, args, "factors", ""),
        "contrast_mode": get_param(cfg, args, "contrast_mode", "pairwise"),
        "contrast": get_param(cfg, args, "contrast", ""),
        "de_test": get_param(cfg, args, "de_test", "Wald"),
        "reduced_design": get_param(cfg, args, "reduced_design", ""),
        "subset_col": get_param(cfg, args, "subset_col", ""),
        "subset_val": get_param(cfg, args, "subset_val", ""),
        "target_genes": get_param(cfg, args, "target_genes", ""),
        "trend_x_col": get_param(cfg, args, "trend_x_col", ""),
        "label_mode": get_param(cfg, args, "label_mode", "sampleID_only"),
        "symbol_col": get_param(cfg, args, "symbol_col", ""),
        "detect_modules": parse_bool(get_param(cfg, args, "detect_modules", False), default=False),
        "trait_cols": get_param(cfg, args, "trait_cols", ""),
        "explain_results": parse_bool(get_param(cfg, args, "explain_results", False), default=False),
    }
    if params["gsea_geneset"] not in {"GO", "KEGG", "Reactome", "all"}:
        LOG.warning("gsea_geneset が不正なため 'all' に変更します: %s", params["gsea_geneset"])
        params["gsea_geneset"] = "all"
    if params["design_mode"] not in {"simple", "additive", "interaction"}:
        LOG.warning("design_mode が不正なため 'simple' に変更します: %s", params["design_mode"])
        params["design_mode"] = "simple"
    if params["contrast_mode"] not in {"pairwise", "auto", "manual"}:
        LOG.warning("contrast_mode が不正なため 'pairwise' に変更します: %s", params["contrast_mode"])
        params["contrast_mode"] = "pairwise"
    if params["de_test"] not in {"Wald", "LRT"}:
        LOG.warning("de_test が不正なため 'Wald' に変更します: %s", params["de_test"])
        params["de_test"] = "Wald"
    if params["gsea_ranking_metric"] not in {"log2FoldChange", "stat", "custom"}:
        LOG.warning("gsea_ranking_metric が不正なため 'log2FoldChange' に変更します: %s",
                    params["gsea_ranking_metric"])
        params["gsea_ranking_metric"] = "log2FoldChange"
    if params["label_mode"] not in {"sampleID_only", "symbol_only", "both"}:
        LOG.warning("label_mode が不正なため 'sampleID_only' に変更します: %s", params["label_mode"])
        params["label_mode"] = "sampleID_only"

    get_organism_args(cfg, args, interactive=not args.non_interactive, config_path=args.config)
    run_selection = args.run if args.run else ("interactive" if not args.non_interactive else "none")

    if args.non_interactive:
        if args.run:
            run_sel = args.run.strip().lower()
            run_sel = STEP_ALIASES.get(run_sel, run_sel)
            if run_sel == "all":
                run_list = list(range(1, 16))
            elif run_sel in {1, 6, 7}:
                run_list = [run_sel]
            else:
                LOG.error("--run は all / 1 / 6 / 7 のみ対応しています。")
                sys.exit(1)

            run_list = filter_openmp_steps(run_list)
            if not run_list:
                LOG.warning("実行可能な解析がありません。")
                return

            dispatch = {
                1: run_01_de, 2: run_02_clustering, 3: run_03_dimreduc,
                4: run_04_normalization, 5: run_05_expression_pattern,
                6: run_06_enrichment, 7: run_07_gsea, 8: run_08_batch,
                9: run_09_regression, 10: run_10_qc, 11: run_11_vis,
                12: run_12_phenotype, 13: run_13_ml, 14: run_14_network,
                15: run_15_sc,
            }

            write_run_metadata(outdir, args, cfg, ORG_ARGS, params, run_list)

            for n in run_list:
                LOG.info("Non-interactive run: step %s", n)
                ok = dispatch[n](common, params, non_interactive=True, dry_run=args.dry_run)
                if not ok:
                    LOG.error("step %s failed. Exiting.", n)
                    sys.exit(1)
            return
        else:
            write_run_metadata(outdir, args, cfg, ORG_ARGS, params, run_selection)
            LOG.info("Non-interactive mode: configuration resolved. Exiting without prompts.")
            return
    else:
        write_run_metadata(outdir, args, cfg, ORG_ARGS, params, run_selection)

    dispatch = {
        "1": run_01_de, "2": run_02_clustering, "3": run_03_dimreduc,
        "4": run_04_normalization, "5": run_05_expression_pattern,
        "6": run_06_enrichment, "7": run_07_gsea, "8": run_08_batch,
        "9": run_09_regression, "10": run_10_qc, "11": run_11_vis,
        "12": run_12_phenotype, "13": run_13_ml, "14": run_14_network,
        "15": run_15_sc,
    }

    while True:
        LOG.info("\n%s", MENU)
        choice = input("  解析番号を選択 > ").strip().lower()

        if choice == "q":
            LOG.info("終了します。お疲れ様でした。")
            break
        elif choice == "b":
            run_batch(common, params, dry_run=args.dry_run)
        elif choice == "s":
            LOG.info("Count matrix: %s", common["counts"])
            LOG.info("Metadata:     %s", common["metadata"])
            LOG.info("Output dir:   %s", common["outdir"])
            if ORG_ARGS:
                LOG.info("Organism:     %s (KEGG %s)", ORG_ARGS["organism"], ORG_ARGS["kegg_org"])
            ch = input_choice("  パスを変更しますか？", ["y", "n"])
            if ch == "y":
                cfg.pop("counts_path", None)
                cfg.pop("metadata_path", None)
                common = get_common_args(cfg, args, interactive=True, config_path=args.config)
                if common is None:
                    sys.exit(1)
                common["outdir"] = str(outdir)
            ch_org = input_choice("  生物種を変更しますか？", ["y", "n"])
            if ch_org == "y":
                get_organism_args(cfg, args, interactive=True, config_path=args.config)
        elif choice in dispatch:
            step_num = int(choice)
            if not OPENMP_AVAILABLE and step_num in OPENMP_STEPS:
                LOG.warning("OpenMP runtime failure detected (OMP Error #179)")
                LOG.warning("Skipping steps: 05, 06, 07, 13, 14, 15")
                continue
            try:
                dispatch[choice](common, params, dry_run=args.dry_run)
            except KeyboardInterrupt:
                LOG.warning("中断されました")
        else:
            LOG.warning("無効な選択です。")


if __name__ == "__main__":
    main()
