# RNA-seq Pipeline 配布ガイド

## 配布者（開発者）向け：リポジトリとイメージの更新手順

このセクションは、配布者がホスト側の開発環境（M5 Mac など）で最終成果物を公開するための運用手順です。  
Trial #17 で確立した `ggtree` の SHA 固定戦略と検証結果を維持したまま、Git リポジトリと Docker イメージを同じバージョンでそろえて公開してください。

### Step A: Git Push

リリース前に、更新済みのコード・マニフェスト・インストーラ・ドキュメントをリモートへ反映します。  
最低限、`dependencies/r-package-manifest.csv`、`R/install_exact_dependencies.R`、`Dockerfile`、配布ドキュメント類がコミットに含まれていることを確認してください。

```bash
git add .
git commit -m "Release v1.0.0: SHA-pinned ggtree and validated pipeline"
git push origin main
```

必要に応じて、あわせて Git タグ `v1.0.0` も push します。

```bash
git tag v1.0.0
git push origin v1.0.0
```

### Step B: Docker Image Tagging & Pushing

ホスト上で検証済みのローカルイメージ `rnaseq-pipeline:v1.0.0` を、そのままコンテナレジストリへ公開します。  
以下は GitHub Container Registry（GHCR）の例です。Docker Hub を使う場合も、`docker tag` と `docker push` の流れは同じです。

```bash
docker login
docker tag rnaseq-pipeline:v1.0.0 ghcr.io/<USERNAME>/rnaseq-pipeline:v1.0.0
docker push ghcr.io/<USERNAME>/rnaseq-pipeline:v1.0.0
```

公開後、利用者向け README にはこのレジストリタグを記載し、`docker pull ghcr.io/<USERNAME>/rnaseq-pipeline:v1.0.0` で取得できる状態にします。

### Step C: Versioning

Docker タグと Git タグ / GitHub Release の版番号は必ず一致させてください。  
たとえばソースコードが `v1.0.0` なのに、配布イメージが `latest` や `r4.4.0` だけだと、利用者はどのコードに対応したコンテナか判断できません。

推奨運用:

- Git タグ: `v1.0.0`
- GitHub Release: `v1.0.0`
- Docker タグ: `rnaseq-pipeline:v1.0.0`
- GHCR タグ: `ghcr.io/<USERNAME>/rnaseq-pipeline:v1.0.0`

この整合を保つことで、Trial #17 で検証済みの内容と配布物を一対一で追跡できます。

## 1. 前提

このパイプラインは Docker を前提に配布します。  
推奨構成は以下です。

- Windows 11 + Docker Desktop
- macOS 13 以降 + Docker Desktop
- Ubuntu 22.04 / 24.04 + Docker Engine + Docker Compose Plugin

再現性の基準環境は以下です。

- R 4.4.0
- Bioconductor 3.20
- Python 3
- 固定 Python 依存: `PyYAML==6.0.3`

## 2. リポジトリ取得

新しい PC で任意の作業ディレクトリを開き、以下を実行します。

```bash
git clone https://github.com/<your-org>/RNA-seq2.git
cd RNA-seq2
```

GitHub から ZIP をダウンロードして展開しても構いませんが、更新追従のため `git clone` を推奨します。

## 3. Docker のインストール

### Windows / macOS

Docker Desktop をインストールします。

- https://www.docker.com/products/docker-desktop/

インストール後、以下で動作確認します。

```bash
docker --version
docker compose version
```

### Linux

Docker Engine と Compose Plugin を導入します。Ubuntu の例:

```bash
sudo apt-get update
sudo apt-get install -y docker.io docker-compose-plugin
sudo systemctl enable --now docker
docker --version
docker compose version
```

必要に応じて現在ユーザーを `docker` グループへ追加してください。

```bash
sudo usermod -aG docker "$USER"
newgrp docker
```

## 4. 入力データの配置

ホスト側で以下のようなディレクトリを用意します。

```bash
mkdir -p data output
```

以下を `data/` に配置します。

- カウント行列 `.rds`
- メタデータ `.csv` または `.xlsx`

例:

```text
RNA-seq2/
├── data/
│   ├── counts.rds
│   └── metadata.csv
└── output/
```

## 5. イメージをビルドする方法

ソースからビルドする場合:

```bash
docker build -t rnaseq-pipeline:v1.0.0 .
```

または Compose を使う場合:

```bash
docker compose build
```

初回ビルドは時間がかかります。R / Bioconductor の固定依存を順番に導入するため、数十分かかる場合があります。

このリポジトリでは再現性のため、通常の `install.packages()` 連打ではなく以下の固定戦略を使っています。

- `dependencies/r-package-manifest.csv` にローカル環境から採取した R パッケージ版を保存
- `R/install_exact_dependencies.R` がその順序で厳密インストールとバージョン検証を実施
- `ggtree` だけは Bioconductor 3.20 のリリース版ではなく、互換性修正を含む GitHub 固定コミット `1a3d0285536625825ccd597d42ca66a6ce8c17c5` を使用
- `ggtree` は `Version == 4.1.1.006` と `RemoteSha == 1a3d0285536625825ccd597d42ca66a6ce8c17c5` の両方を検証

つまり、ビルドが通った Docker イメージ自体が配布用の再現基準になります。

## 6. イメージを pull する方法

GitHub Container Registry などに配布済みイメージがある場合は、ビルドせず pull できます。

```bash
docker pull ghcr.io/<your-org>/rnaseq-pipeline:v1.0.0
```

pull 後の実行例:

```bash
docker tag ghcr.io/<your-org>/rnaseq-pipeline:v1.0.0 rnaseq-pipeline:v1.0.0
```

## 7. ヘルプ表示

```bash
docker run --rm rnaseq-pipeline:v1.0.0 --help
```

または Compose:

```bash
docker compose run --rm pipeline --help
```

## 7.1 ビルド後の推奨検証

ビルド直後は、最低でもヘルプ表示と example データ実行を確認します。

```bash
docker run --rm rnaseq-pipeline:v1.0.0 --help
```

```bash
docker run --rm rnaseq-pipeline:v1.0.0 \
  --example \
  --organism zebrafish \
  --run 1 \
  --non-interactive
```

期待される確認ポイント:

- example データ生成が完了する
- `zebrafish -> OrgDb=org.Dr.eg.db, KEGG=dre` が表示される
- `00_align_check.R` が完了する
- `01_differential_expression.R` が完了し、synthetic data に対して DESeq2 の結果が出力される

## 8. 自分のデータで実行する方法

### `docker run` の例

```bash
docker run --rm \
  -v "$(pwd)/data:/data" \
  -v "$(pwd)/output:/output" \
  rnaseq-pipeline:v1.0.0 \
  --counts /data/counts.rds \
  --metadata /data/metadata.csv \
  --organism zebrafish \
  --run all \
  --non-interactive \
  --outdir /output/run_01
```

### `docker compose` の例

`docker-compose.yml` を使う場合:

```bash
docker compose run --rm pipeline \
  --counts /data/counts.rds \
  --metadata /data/metadata.csv \
  --organism zebrafish \
  --run all \
  --non-interactive \
  --outdir /output/run_01
```

## 9. 追加オプションの実行例

### LRT を使った差次的発現解析

注意: この CLI では `--test` ではなく `--de-test` を使います。

```bash
docker compose run --rm pipeline \
  --counts /data/counts.rds \
  --metadata /data/metadata.csv \
  --organism zebrafish \
  --run 1 \
  --non-interactive \
  --design-mode additive \
  --factors dpf,genotype \
  --de-test LRT \
  --reduced-design "~ dpf" \
  --outdir /output/lrt_run
```

### サブセット解析

```bash
docker compose run --rm pipeline \
  --counts /data/counts.rds \
  --metadata /data/metadata.csv \
  --organism zebrafish \
  --run 1 \
  --non-interactive \
  --subset-col genotype \
  --subset-val WT,Het \
  --outdir /output/subset_run
```

### GSEA のランキング指標指定

```bash
docker compose run --rm pipeline \
  --counts /data/counts.rds \
  --metadata /data/metadata.csv \
  --organism zebrafish \
  --run 7 \
  --non-interactive \
  --gsea-geneset all \
  --gsea-ranking-metric stat \
  --outdir /output/gsea_run
```

### クラスタリング対象遺伝子を指定

`target_genes.csv` は `gene` 列を必須とします。

```bash
docker compose run --rm pipeline \
  --counts /data/counts.rds \
  --metadata /data/metadata.csv \
  --organism zebrafish \
  --run 2 \
  --non-interactive \
  --target-genes /data/target_genes.csv \
  --outdir /output/clustering_run
```

### ネットワーク解析でモジュール検出を有効化

```bash
docker compose run --rm pipeline \
  --counts /data/counts.rds \
  --metadata /data/metadata.csv \
  --organism zebrafish \
  --run 14 \
  --non-interactive \
  --detect-modules TRUE \
  --trait-cols condition,dpf,batch \
  --outdir /output/network_run
```

## 10. 設定ファイルを使う方法

テンプレートをコピーして編集します。

```bash
cp config/config_template.yaml config.yaml
```

設定ファイルを使う実行例:

```bash
docker compose run --rm pipeline \
  --config /app/config.yaml \
  --run all \
  --non-interactive
```

Compose の既定ボリューム設定ではリポジトリ全体が `/app` にマウントされるため、`config.yaml` はそのまま利用できます。

## 11. 出力の確認

実行後の結果はホスト側 `output/` に保存されます。

例:

```text
output/
└── run_01/
    ├── 01_DE/
    ├── 06_enrichment/
    ├── 07_gsea/
    ├── pipeline.log
    └── environment_report.txt
```

## 12. よくあるトラブル

### Docker が起動していない

```bash
docker ps
```

失敗する場合は Docker Desktop または Docker Engine を起動してください。

### 入力ファイルが見つからない

コンテナ内ではホストの `./data` が `/data` にマウントされます。  
そのため、CLI 引数はホストパスではなくコンテナパスを指定します。

正しい例:

```bash
--counts /data/counts.rds --metadata /data/metadata.csv
```

### `.xlsx` を使う場合

メタデータが Excel の場合はそのまま `/data/metadata.xlsx` を指定できます。

### 古いビルドキャッシュを消したい

```bash
docker builder prune -f
docker image rm rnaseq-pipeline:v1.0.0
```

## 13. 推奨運用

- 本番解析前に `--help` と小規模データで動作確認する
- 出力先ごとに `--outdir` を分ける
- 解析条件を残すため `config.yaml` をプロジェクトごとに保存する
- zebrafish 解析では `--organism zebrafish` を明示する
