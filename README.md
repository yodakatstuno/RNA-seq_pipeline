# RNA-seq統合解析パイプライン

Docker と厳密依存管理による再現可能な RNA-seq パイプラインです。差次的発現解析、機能富化解析、GSEA、次元削減、可視化などを一貫して実行できます。

**対応生物種:** ヒト・マウス・ゼブラフィッシュ

## クイックスタート（Docker）

検証済みの `v1.0.0` イメージを使う最短手順です。

```bash
git clone https://github.com/yodakatstuno/RNA-seq2.git
cd RNA-seq2

docker pull ghcr.io/yodakatstuno/rnaseq-pipeline:v1.1.0
docker tag ghcr.io/yodakatstuno/rnaseq-pipeline:v1.1.0 rnaseq-pipeline:v1.1.0

docker run --rm rnaseq-pipeline:v1.1.0 --help

docker run --rm rnaseq-pipeline:v1.1.0 \
  --example \
  --organism zebrafish \
  --run 1 \
  --non-interactive
```

## クイックスタート（ローカル実行）

```bash
git clone https://github.com/yodakatstuno/RNA-seq2.git
cd RNA-seq2
conda env create -f environment.yml
conda activate rnaseq_env
Rscript R/install_missing_packages.R
python main_pipeline.py --example
```

---

## 目次

- [インストール](#インストール)
- [入力形式](#入力形式)
- [使い方](#使い方)
- [解析ステップ](#解析ステップ)
- [CLI リファレンス](#cli-リファレンス)
- [出力構成](#出力構成)
- [Docker と再現性](#docker-と再現性)
- [設定ファイル](#設定ファイル)
- [トラブルシューティング](#トラブルシューティング)
- [ライセンス](#ライセンス)

---

## インストール

### 方法1: Docker（推奨）

[Docker](https://docs.docker.com/get-docker/) をインストール済みであることを前提とします。

配布済みイメージを使う場合:

```bash
docker pull ghcr.io/yodakatstuno/rnaseq-pipeline:v1.1.0
docker tag ghcr.io/yodakatstuno/rnaseq-pipeline:v1.1.0 rnaseq-pipeline:v1.1.0
```

ソースからビルドする場合:

```bash
docker compose build
```

これにより、R / Python / Bioconductor の依存関係がコンテナ内へ事前導入されます。

### 方法2: Conda

```bash
conda env create -f environment.yml
conda activate rnaseq_env
Rscript R/install_missing_packages.R
```

### 方法3: 手動セットアップ

1. R 4.3 以上と Python 3.9 以上をインストールします。
2. `pip install pyyaml openpyxl` を実行します。
3. `R/install_missing_packages.R` に記載された R パッケージを導入します。

---

## 入力形式

### カウント行列（`.rds`）

| | sample_01 | sample_02 | sample_03 | ... |
|---|---|---|---|---|
| Gene_001 | 1523 | 892 | 2301 | ... |
| Gene_002 | 45 | 38 | 52 | ... |

- **行:** 遺伝子 ID（ENSEMBL、SYMBOL、ENTREZ など）
- **列:** サンプル ID（メタデータと一致する必要があります）
- **値:** 正規化前の raw integer count

### メタデータ（`.xlsx` または `.csv`）

| sample_id | condition | genotype | dpf | batch |
|---|---|---|---|---|
| sample_01 | control | WT | 3 | A |
| sample_02 | control | WT | 5 | A |
| sample_03 | treatment | Ho | 3 | B |

#### 列の説明

| 列名 | 必須 | 説明 |
|---|---|---|
| `sample_id` | ✅ | カウント行列の列名と一致する必要があります。自動検出にも対応します。 |
| `condition` | ✅ | 差次的発現解析で使う主比較条件です。 |
| `genotype` | 任意 | 生物学的群情報です（例: WT, Ho, Het）。 |
| `dpf` | 任意 | 発生ステージです。 |
| `batch` | 任意 | バッチや反復情報です。 |

> **注意:** サンプル ID 列は `ID`、`sample`、`sample_id`、`sampleID`、`index`、`name` などの一般的な列名から自動検出します。該当列がない場合は、カウント行列の列名との一致率が最も高い列を使用します。

---

## 使い方

### 対話モード

```bash
python main_pipeline.py
```

メニュー形式で順に解析条件を設定できます。

### 非対話モード

```bash
python main_pipeline.py \
  --counts /path/to/counts.rds \
  --metadata /path/to/meta.xlsx \
  --organism zebrafish \
  --run all \
  --non-interactive
```

高度な指定を含む例:

```bash
python main_pipeline.py \
  --counts /path/to/counts.rds \
  --metadata /path/to/meta.xlsx \
  --organism zebrafish \
  --run 1 \
  --non-interactive \
  --design-mode additive \
  --factors dpf,genotype \
  --de-test LRT \
  --reduced-design "~ dpf" \
  --subset-col genotype \
  --subset-val WT,Het \
  --analysis-name 01_DE_genotype_subset
```

### サンプルデータ実行モード

合成データを自動生成してパイプラインを実行します。

```bash
python main_pipeline.py --example
python main_pipeline.py --example --run all
```

### Docker での実行

自分のデータを使う場合:

```bash
docker compose run pipeline \
  --counts /data/counts.rds \
  --metadata /data/meta.xlsx \
  --organism zebrafish \
  --run all \
  --non-interactive
```

ホスト側では `./data/` に入力ファイルを置き、結果はリポジトリ内の `./output/result/` に出力されます。

### 検証済みサンプル実行

`v1.0.0` イメージで検証済みのコマンドです。

```bash
docker run --rm rnaseq-pipeline:v1.0.0 \
  --example \
  --organism zebrafish \
  --run 1 \
  --non-interactive
```

この実行では zebrafish 設定で `00_align_check.R` と `01_differential_expression.R` が完了し、synthetic data から DESeq2 により 20 個の DEGs が得られることを確認済みです。

---

## 解析ステップ

| ステップ | 名称 | 説明 |
|---|---|---|
| 1 | **差次的発現解析** | DESeq2 / edgeR / limma による DE 遺伝子検出 |
| 2 | **サンプルクラスタリング** | 階層クラスタリングと k-means |
| 3 | **次元削減** | PCA / t-SNE / UMAP と診断解析 |
| 4 | **正規化** | DESeq2 / TMM / log2 正規化比較 |
| 5 | **発現パターン解析** | 遺伝子発現パターン解析 |
| 6 | **機能富化解析** | GO / KEGG / Reactome の ORA |
| 7 | **GSEA** | Gene Set Enrichment Analysis |
| 8 | **バッチ補正** | ComBat / ComBat-seq による補正 |
| 9 | **回帰解析** | 発現量の回帰解析 |
| 10 | **サンプル QC** | 品質評価指標の集計 |
| 11 | **可視化** | 発現データの総合可視化 |
| 12 | **表現型関連解析** | 表現型と発現の関連評価 |
| 13 | **機械学習** | RF / SVM 分類 |
| 14 | **ネットワーク解析** | 共発現ネットワーク解析 |
| 15 | **単一細胞統合** | bulk RNA-seq と scRNA-seq の統合 |

---

## CLI リファレンス

### 基本オプション

| オプション | 説明 | 既定値 |
|---|---|---|
| `--counts <path>` | カウント行列（`.rds`） | なし |
| `--metadata <path>` | メタデータ（`.xlsx` / `.csv`） | なし |
| `--organism <name>` | `human`、`mouse`、`zebrafish` | なし |
| `--config <path>` | 設定 YAML ファイル | `config.yaml` |
| `--run <step>` | `all`、`1`–`15`、`de`、`enrichment`、`gsea` | なし |
| `--outdir <path>` | 互換性維持のため受理。実際の出力先は常に `output/result` | `output/result` |
| `--non-interactive` | 対話なしで実行 | off |
| `--dry-run` | コマンドのみ表示 | off |
| `--example` | 内蔵 example データで実行 | off |
| `--describe` | ステップ説明を表示 | off |
| `--list-steps` | ステップ一覧を表示 | off |

### 差次的発現解析オプション（Step 1）

| オプション | 説明 | 既定値 |
|---|---|---|
| `--fdr-cutoff` | FDR 閾値 | `0.05` |
| `--lfc-cutoff` | log2FC 閾値 | `1.0` |
| `--design-mode` | `simple`、`additive`、`interaction` | `simple` |
| `--factors` | 因子列をカンマ区切りで指定 | なし |
| `--contrast-mode` | `pairwise`、`auto`、`manual` | `pairwise` |
| `--analysis-name` | 出力サブフォルダ名 | `01_DE` |
| `--de-test` | DESeq2 検定法: `Wald` または `LRT` | `Wald` |
| `--reduced-design` | `--de-test LRT` 時の縮約モデル | 空 |
| `--subset-col` | DE 前に絞り込みに使うメタデータ列 | 空 |
| `--subset-val` | `--subset-col` で残す値をカンマ区切りで指定 | 空 |

#### デザイン式の例

```text
simple:       ~ condition
additive:     ~ dpf + genotype
interaction:  ~ dpf + genotype + dpf:genotype
```

#### 高度な DE 解析の補足

- `--de-test LRT` を指定すると、DESeq2 は係数検定（Wald）ではなくモデル比較（LRT）を実行します。
- `--reduced-design` は `--de-test LRT` と組み合わせて必須です。`~ dpf` や `~ genotype` のような R の式文字列を指定してください。
- `--subset-col` と `--subset-val` を使うと、元の入力ファイルを編集せずに、特定サブセットだけで DE 解析を実行できます。
- `--subset-val` に複数値を渡す場合は `WT,Het` のようにカンマ区切りで指定します。

### 富化解析 / GSEA オプション（Step 6-7）

| オプション | 説明 | 既定値 |
|---|---|---|
| `--enrich-pval-cutoff` | 富化解析の p 値閾値 | `0.05` |
| `--gsea-geneset` | `GO`、`KEGG`、`Reactome`、`all` | `all` |
| `--gsea-ranking-metric` | ランキング指標: `log2FoldChange`、`stat`、`custom` | `log2FoldChange` |
| `--gsea-min-gs` | 遺伝子セットの最小サイズ | `15` |
| `--gsea-max-gs` | 遺伝子セットの最大サイズ | `500` |

`--gsea-ranking-metric` は、GSEA 実行前にどの値で遺伝子を順位付けするかを制御します。`stat` は DESeq2 の検定統計量を使いたい場合に適しており、`custom` は上流の DE 結果ファイルに必要な独自ランキング列が含まれている場合にのみ使用してください。

### クラスタリング / ネットワーク解析オプション（Step 2 と Step 14）

| オプション | 対象ステップ | 説明 | 既定値 |
|---|---|---|---|
| `--target-genes` | 2 | `gene` 列を持つ CSV を指定し、その遺伝子のみクラスタリング | 空 |
| `--trend-x-col` | 2 | クラスタ trend plot の x 軸に使う metadata 列 | 空 |
| `--detect-modules` | 14 | ネットワーク解析で WGCNA モジュール検出を有効化する `TRUE/FALSE` | `FALSE` |
| `--trait-cols` | 14 | モジュール相関に使うメタデータ列をカンマ区切りで指定 | 空 |

- `--target-genes` は `gene` 列を持つ CSV を受け取り、既定では `padj < 0.05` の遺伝子だけをクラスタリングに使います。
- `--trend-x-col` を指定すると、k-means クラスタごとの平均発現トレンドを描画します。列が存在しない場合はエラーになります。
- `--detect-modules TRUE` は Step 14 の WGCNA 系処理を有効化します。
- `--trait-cols` には、入力メタデータ内に存在する列名を指定してください。例: `condition,dpf,batch`

### サンプルラベルオプション（Steps 2 / 3 / 4 / 8 / 11）

| オプション | 説明 | 既定値 |
|---|---|---|
| `--label-mode` | `sampleID_only`、`symbol_only`、`both` | `sampleID_only` |
| `--symbol-col` | `symbol_only` または `both` で使う metadata 列 | 空 |

- `sampleID_only` は従来どおりサンプル ID を表示します。
- `symbol_only` は `--symbol-col` の値だけを表示します。
- `both` は `sampleID | symbol` の形式で表示します。
- PCA、k-means、バッチ補正、各種 heatmap ではラベル数に応じて文字サイズを下げ、散布図では `ggrepel` により重なりを抑えます。

### 次元削減オプション（Step 3）

対話モードでは以下を順に選択できます。

- **変換法:** `log2` / `vst` / `rlog`
- **特徴量選択:** `variance` / `mad`
- **多因子 PCA:** メタデータ列ごとの色分け
- **分割解析:** 因子レベルごとの PCA
- **寄与率評価:** メタデータ因子と PCA 成分の R²
- **パラメータ探索:** 複数の perplexity / n_neighbors で比較

---

## 出力構成

```text
output/
└── result/
    ├── 01_DE/
    │   ├── de_results.csv
    │   ├── results/de_results.csv
    │   ├── plots/volcano_plot.pdf
    │   └── logs/sessionInfo.txt
    ├── 02_clustering/
    │   ├── dendrogram.pdf
    │   ├── distance_heatmap.pdf
    │   └── cluster_expression_trend.pdf
    ├── 03_dimreduc/
    │   ├── pca_plot.pdf
    │   ├── pca_variance_explained.csv
    │   ├── sample_distance_heatmap.pdf
    │   └── pca_by_factor/
    ├── 06_enrichment/
    │   ├── GO_BP_results.csv
    │   └── KEGG_dotplot.pdf
    ├── 07_gsea/
    │   ├── GSEA_GO_results.csv
    │   └── GSEA_GO_dotplot.pdf
    └── pipeline.log
```

---

## Docker と再現性

### `docker-compose.yml`

ホスト側のデータと出力を以下のようにマウントして利用します。

```yaml
volumes:
  - ./data:/data
  - ./output:/output
```

解析結果はコンテナ内外ともに `output/result/` 配下へ集約されます。

### 変更後に再ビルドする

```bash
docker compose build --no-cache
```

### 再現性の考え方

このプロジェクトでは、「動いたときの依存状態」を Docker イメージに固定することを重視しています。

- `dependencies/r-package-manifest.csv` にホスト環境から採取した厳密な R パッケージ版を保存します。
- `R/install_exact_dependencies.R` がその順序で逐次インストールし、導入後に版番号を検証します。
- `ggtree` は moving target を避けるため、Bioconductor の通常リリースではなく GitHub 固定コミット `1a3d0285536625825ccd597d42ca66a6ce8c17c5` を使用します。
- `ggtree` は `Version == 4.1.1.006` と `RemoteSha == 1a3d0285536625825ccd597d42ca66a6ce8c17c5` の両方を検証します。

このため、`v1.0.0` の Docker イメージは「検証済みの実行環境」そのものとして扱えます。

---

## 設定ファイル

テンプレートをコピーして利用してください。

```bash
cp config/config_template.yaml config.yaml
```

その後、`config.yaml` に入力ファイルパスや解析条件を記述します。利用可能な項目は `config/config_template.yaml` を参照してください。

`outdir` は互換性維持のため設定ファイルに残っていますが、実際の出力先は常にリポジトリ直下の `output/result` です。

---

## トラブルシューティング

| エラー | 対処法 |
|---|---|
| `Count matrix が見つかりません` | `--counts` のパスを確認してください |
| `Metadata が見つかりません` | `--metadata` のパスを確認してください |
| `some values in assay are not integers` | カウント値は自動的に丸められます。入力データの元形式を確認してください |
| R パッケージが不足している | `Rscript R/install_missing_packages.R` を実行してください |
| `python not found` | `python3` を使うか conda 環境を有効化してください |
| `OMP: Error #179`（macOS） | `brew install libomp` を実行するか Docker を使ってください |
| Reactome がスキップされる | Reactome はヒト専用です |
| 遺伝子 ID 変換率が低い | 遺伝子 ID 形式（ENSEMBL / SYMBOL など）を確認してください |

---

## ライセンス

MIT License
