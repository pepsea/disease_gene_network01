# NW重複評価アプリ (Gene–Disease Network Overlap Evaluator)

疾患と遺伝子リストを入力すると、遺伝子ごとに PPI ネットワークを構築し、
その近傍が**疾患ネットワーク（Open Targets 上位遺伝子群）を何％カバーするか**を
評価して表にします。

1. **疾患** — Open Targets を検索し、候補から選んで EFO/MONDO ID を確定
2. **遺伝子** — HGNC (HUGO) の公式シンボルと照合（別名・旧シンボルは自動変換）
3. **PPI** — SIGNOR / STRING / BioGRID の組み合わせとスコア条件を選択
4. **解析** — 2つの表で重複を評価
   - **表1 遺伝子重複** — PPIパートナーと疾患遺伝子の重複
   - **表2 パスウェイ重複** — 疾患のエンリッチメント解析結果との重複

PPI の収集・統合・ランキングは
[aiagent_hypothesis_generator002](https://github.com/pepsea/aiagent_hypothesis_generator002)
の解析法をそのまま移植しています。

## Docker で常駐させる（推奨）

```bash
cp .env.example .env      # 必要に応じて編集（BioGRID キーなど）
docker compose up -d --build
```

http://127.0.0.1:5005 を開きます。

`restart: unless-stopped` を設定しているため、コンテナが落ちた場合もホストを
再起動した場合も自動で復帰します（`docker compose down` で明示的に停止した
場合は復帰しません）。

```bash
docker compose ps                  # 稼働状態と healthcheck の結果
docker compose logs -f             # ログ追尾
docker compose up -d --build       # コード更新後の再デプロイ
docker compose down                # 停止・削除（常駐解除）
```

ヘルスチェックは 30 秒間隔で `/healthz` を叩き、`docker compose ps` の
`STATUS` 列に `healthy` / `unhealthy` として出ます。ログは 10MB × 3 世代で
ローテートします。PPI キャッシュは名前付きボリューム `ppi-cache` に保存され、
再起動・再ビルドをまたいで保持されます。

ポートを変えたい場合は `.env` の `HOST_PORT` を変更します
（コンテナ内は常に 5005 を listen します）。

### ローカルで直接動かす場合

```bash
pip install -r requirements.txt
python app.py                                  # 開発用サーバ
gunicorn -c gunicorn.conf.py app:app           # 本番同等
```

## 使い方

### STEP 1 — 疾患を検索して指定

疾患名（英語）を入力して「Open Targets で検索」を押すと候補が出ます。
完全一致がある場合は自動で選択されます（`Alzheimer disease` が
`Alzheimer disease 2` に負けないようにするため）。
`EFO_0000249` のようなオントロジーIDを直接入力すると検索を省略します。

### STEP 2 — 遺伝子を HGNC 照合

「HGNCシンボルを確認」で各遺伝子の状態が表示されます。

| 状態 | 意味 | 解析での扱い |
|---|---|---|
| 承認シンボル | HGNC の approved symbol | そのまま解析 |
| 別名 / 旧シンボル | alias_symbol / prev_symbol に一致 | 承認シンボルに変換して解析 |
| 廃止 | Entry Withdrawn | 警告のうえ解析 |
| HGNC未登録 | どれにも一致しない | 警告のうえ解析（新しいシンボルの可能性があるためブロックしない） |
| 照合できず | HGNC に到達できなかった | 入力のまま解析 |

別名が承認シンボルに変換された結果、同じ遺伝子が重複した場合は1件にまとめます。

### STEP 3 — PPI ソースとスコア条件

| 設定 | 既定値 | 説明 |
|---|---|---|
| ソース | SIGNOR + STRING | SIGNOR / STRING / BioGRID を自由に組み合わせ |
| STRING 信頼度閾値 | 700 | STRING 側の `required_score`（0–1000、400=中, 700=高） |
| 共通スコア下限 | なし | 全ソース共通のエッジスコア下限（0–1） |
| ハブ判定の相互作用数 | 1000 | IntAct のグローバル相互作用数がこれを超えると非特異ハブとして除外 |
| ハブ遺伝子を除外する | ON | 上記に加え、ユビキチン・アクチン・チューブリン等の非特異ハブ族も除外 |
| PPIパートナー上限 | 100 | ランキング上位から採用する件数 |

BioGRID は API キーが必要です。未設定の場合はチェックボックスが無効になり、
リクエストで指定されても自動的に外されます（`biogrid_skipped: true` で通知）。

## PPI 解析法（移植元の仕様）

各ソースの相互作用を1つのグラフに統合し、パートナーを次の優先順で並べます。

1. **裏付けDB数** — 複数のDBが支持する相互作用ほど再現性が高い
2. **代表スコア** — `SIGNOR > BioGRID > STRING` の順で最初に見つかったDB固有スコア
3. **エッジの重み** — 観測された相互作用の回数

スコアを持たないエッジ（BioGRID の `SCORE` は多くが空）には、
パートナーの**グローバル相互作用数の逆数**を代用スコアとして与えます。
相互作用数が少ない＝無差別なハブでないパートナーほど、その相互作用は特異的で
意味がある、という考え方です。

ソース別の取得仕様:

| ソース | 取得方法 | 除外・整形 |
|---|---|---|
| SIGNOR | ヒト全データをTSVで一括取得（7日キャッシュ） | パートナーが protein / complex 以外（化合物・フェノタイプ等）を除外 |
| STRING | `required_score` 付きで遺伝子ごとに取得（3日キャッシュ） | パートナー同士のエッジを除外、パートナーごとに最高スコアの1件のみ |
| BioGRID | 遺伝子ごとに取得（`selfInteractionsExcluded`） | 実験手法・文献違いの重複をパートナーごとに1件へ集約 |

## 2つの表

### 表1 — 疾患ネットワーク重複（遺伝子レベル）

遺伝子のPPIパートナーが、疾患の Open Targets 上位遺伝子を何％カバーするか。

### 表2 — 疾患エンリッチメント重複（パスウェイレベル）

疾患遺伝子を g:Profiler でエンリッチメント解析して「疾患パスウェイ signature」を
作り、各遺伝子の**PPIネットワーク（自身＋パートナー）のエンリッチメント結果**が
それを何％再現するかを見ます。

遺伝子レベルでは重複しなくても、同じパスウェイ上にいる場合があります。
表2はその機能的な近さを拾います。

| 指標 | 定義 |
|---|---|
| 加重パスウェイ重複率 | 重複パスウェイの `-log10(p)` 合計 ÷ 疾患パスウェイ全体の `-log10(p)` 合計。有意なパスウェイほど重く数えます |
| 単純パスウェイ重複率 | 重複パスウェイ数 ÷ 疾患パスウェイ数 |
| ターゲット自身の該当 | 疾患パスウェイのうち、ターゲット遺伝子自身が構成遺伝子として含まれる数（元プロジェクトの `pathway_fit` に相当） |

エンリッチメント解析は STEP 3 のチェックボックスで ON/OFF できます（既定 ON）。
OFF にすると表2は表示されず、g:Profiler も呼びません。g:Profiler に到達できない
場合も表1はそのまま表示されます。

疾患パスウェイの算出には、**表1と同じ疾患遺伝子リスト**を使います
（`NW_ENRICH_GENE_N` で件数を変更可能）。2つの表が同じ疾患を指すようにするため
です。移植元は 20 遺伝子でエンリッチメントしていました。

## ターゲット遺伝子自身の扱い

両方の表で、ターゲット遺伝子自身の寄与を重複に含めます。

| | ターゲット自身の数え方 |
|---|---|
| 表1 | ターゲットが疾患遺伝子リストに含まれる場合、そのOTスコアを分子に加算し `matched_count` にも計上 |
| 表2 | ターゲットが疾患パスウェイの構成遺伝子である場合、そのパスウェイを重複として計上 |

いずれも**二重計上はしません**。表1ではターゲット自身を重複遺伝子リストから除外し、
表2ではネットワーク経由と自身の所属が同じパスウェイを指す場合も1件として数えます
（`via` が `both`）。

表2で自身を含めるのは重要です。PPIパートナーが少ない遺伝子はエンリッチメント
クエリが実質1遺伝子になり、g:Profiler は単一遺伝子ではほぼ有意なタームを返しません。
自身の所属を合算しないと、疾患パスウェイに含まれる遺伝子でもスコアが 0% になって
しまいます。

`overlapping_pathways` の各要素には、どの経路で重複したかを示す `via` が付きます。

| `via` | 意味 |
|---|---|
| `network` | PPIネットワークのエンリッチメント経由 |
| `target` | ターゲット遺伝子自身が構成遺伝子 |
| `both` | 両方（1件として計上） |

## スコアの定義

```
加重重複率 = ( Σ 重複遺伝子のOTスコア + ターゲット自身のOTスコア )
             / Σ 疾患ネットワーク全体のOTスコア

単純重複率 = 重複した疾患遺伝子数 / 疾患ネットワークの遺伝子数
```

どちらも**疾患ネットワークに対する割合（％）**です。ターゲット自身が疾患遺伝子
リストに含まれる場合、そのスコアが分子に加算されます（重複遺伝子側では二重
計上しません）。表2も同様にターゲット自身の所属パスウェイを含めます
（「ターゲット遺伝子自身の扱い」を参照）。

| 加重重複率 | 解釈 |
|---|---|
| ≥ 30% | ネットワーク的に疾患と強く関連 |
| 10〜30% | 中程度の関連 |
| < 10% | 関連が薄い（または PPI 情報不足） |

## API

### `GET /api/diseases?q=<疾患名>&limit=10`

```jsonc
{ "query": "alzheimer disease",
  "results": [
    { "id": "EFO_0000249", "name": "Alzheimer disease",
      "description": "...", "exact": true }
  ] }
```

### `POST /api/genes/validate`

```jsonc
// リクエスト: { "genes": ["APP", "PS1", "NOTAGENE"] }
{ "genes": [
    { "input": "PS1", "symbol": "PSEN1", "status": "alias",
      "hgnc_id": "HGNC:8828", "name": "presenilin 1" }
  ],
  "summary": { "approved": 1, "alias": 1, "unknown": 1 } }
```

### `POST /api/analyze`

```jsonc
{
  "disease_id": "EFO_0000249",        // STEP 1 で選んだID（`disease` 名でも可）
  "genes": ["APP", "PS1", "APOE"],    // 配列でも改行・カンマ区切り文字列でも可
  "top_n": 100,                        // 疾患上位遺伝子数（任意、上限500）
  "enrichment": true,                  // 表2を計算するか（既定 true）
  "ppi": {
    "sources": ["signor", "string", "biogrid"],
    "string_score": 700,
    "min_score": null,
    "hub_threshold": 1000,
    "exclude_hubs": true,
    "max_nodes": 100
  }
}
```

```jsonc
{
  "disease_id": "EFO_0000249",
  "disease_label": "Alzheimer disease",
  "disease_gene_count": 100,
  "ppi": { "...実際に使われた設定..." },
  "genes": { "resolved": [ /* HGNC照合結果 */ ], "summary": { "approved": 3 } },
  "enrichment": {
    "enabled": true, "disease_pathway_count": 4, "disease_gene_n": 100,
    "top_pathways": [ { "term_id": "R-HSA-977225", "name": "Amyloid fiber formation",
                        "source": "REAC", "p_value": 1e-12 } ]
  },
  "results": [
    {
      "gene": "PSEN1", "input_gene": "PS1", "symbol_status": "alias",
      "hgnc_id": "HGNC:8828", "gene_name": "presenilin 1",
      "weighted_percent": 39.7, "weighted_score": 0.397,
      "overlap_percent": 25.0, "simple_ratio": 0.25,
      "matched_count": 2, "overlap_count": 1,
      "disease_gene_count": 8, "ppi_partner_count": 2,
      "excluded_hub_count": 0, "source_counts": { "SIGNOR": 2 },
      "target_self": "PSEN1", "target_self_score": 0.71,
      "interpretation": "strong",
      "overlapping_genes": [ { "symbol": "APP", "score": 0.92 } ],

      // 表2（enrichment 有効時のみ）
      "pathway_weighted_percent": 58.6, "pathway_overlap_percent": 75.0,
      "pathway_matched_count": 3,     // 重複パスウェイ総数（自身を含む）
      "pathway_overlap_count": 2,     // うちネットワーク経由
      "disease_pathway_count": 4,
      "gene_pathway_count": 3,
      "target_in_pathway_count": 2, "target_fit_percent": 50.0,
      "pathway_interpretation": "strong",
      "overlapping_pathways": [
        { "term_id": "R-HSA-977225", "name": "Amyloid fiber formation",
          "source": "REAC", "p_value": 1e-12, "via": "both" }
      ]
    }
  ]
}
```

`results` は加重重複率の降順です。個別の遺伝子で PPI 取得に失敗した場合、
その要素は `{"gene": ..., "error": ...}` になり、他の遺伝子の結果は返ります。

ステータスコード: `400`（入力不備・遺伝子数超過）、`404`（疾患が見つからない・
関連遺伝子ゼロ）、`502`（Open Targets 側の障害）。

### `GET /healthz`

現在の設定値と稼働状態を返します。

## 設定（環境変数）

| 変数 | 既定値 | 説明 |
|---|---|---|
| `PORT` / `HOST` | `5005` / `127.0.0.1` | 待ち受け |
| `HOST_PORT` | `5005` | Docker でホスト側に公開するポート |
| `BIOGRID_KEY` | （空） | BioGRID API キー（`BIOGRID_API_KEY` でも可） |
| `NW_MAX_GENES` | `30` | 1リクエストあたりの遺伝子数上限 |
| `NW_DISEASE_TOP_N` | `100` | 疾患上位遺伝子の取得件数 |
| `NW_MAX_WORKERS` | `5` | 遺伝子の並列処理数 |
| `NW_STRING_SCORE` | `700` | STRING 信頼度閾値の既定値 |
| `NW_HUB_THRESHOLD` | `1000` | ハブ判定の既定値 |
| `NW_MAX_NODES` | `100` | PPIパートナー上限の既定値 |
| `NW_ENRICH_GENE_N` | `100` | エンリッチメントに使う疾患遺伝子の件数 |
| `NW_ENRICH_MAX_TERMS` | `50` | 取得するパスウェイ数の上限 |
| `PPI_CACHE_DIR` | `./ppi_cache` | PPI キャッシュの保存先 |
| `PPI_CACHE_DISABLED` | 無効 | `1` でディスクキャッシュを無効化 |
| `GUNICORN_WORKERS` / `GUNICORN_THREADS` | `2` / `8` | 常駐時のサーバ規模 |
| `GUNICORN_TIMEOUT` | `300` | リクエストタイムアウト（秒） |
| `LOG_LEVEL` | `info` | ログレベル |

外部APIのURL（`OT_API_URL` / `HGNC_API_URL` / `SIGNOR_TSV_URL` / `STRING_URL` /
`BIOGRID_URL` / `INTACT_URL` / `GPROFILER_URL`）も上書きできます。ミラーの利用やテストに使います。

## データソース

| ソース | 用途 | APIキー | ライセンス |
|---|---|---|---|
| [Open Targets](https://platform.opentargets.org/) | 疾患検索・疾患上位遺伝子 | 不要 | CC0 |
| [HGNC](https://www.genenames.org/) | 遺伝子シンボル照合 | 不要 | CC0 |
| [SIGNOR](https://signor.uniroma2.it/) | 制御性シグナル相互作用 | 不要 | CC BY 4.0 |
| [STRING](https://string-db.org/) | 機能的・物理的関連 | 不要 | CC BY 4.0 |
| [BioGRID](https://thebiogrid.org/) | キュレーション済み相互作用 | **必要** | 非商用・学術利用限定 |
| [IntAct](https://www.ebi.ac.uk/intact/) | ハブ判定用の相互作用数 | 不要 | CC BY 4.0 |
| [g:Profiler](https://biit.cs.ut.ee/gprofiler/) | パスウェイ/GO エンリッチメント | 不要 | BSD 2-Clause |

各ソースは best-effort です。1つが落ちても残りで解析を続行し、
HGNC に到達できない場合も入力シンボルのまま解析を進めます。

## 構成

```
app.py                    # Flask ルーティング・入力検証・並列実行
nw_overlap.py             # 表1: 遺伝子レベルの重複スコア
enrichment_overlap.py     # 表2: パスウェイレベルの重複スコア
ppi_network.py            # PPIグラフ構築・パートナーランキング・ハブ判定
collectors/
├── _http.py              # リトライ付き共有 HTTP セッション
├── _cache.py             # ディスクキャッシュの場所
├── opentargets.py        # 疾患検索・疾患上位遺伝子
├── hgnc.py               # HGNC シンボル照合
├── gprofiler.py          # パスウェイエンリッチメント
├── signor.py             # SIGNOR
├── string_db.py          # STRING
└── biogrid.py            # BioGRID
templates/index.html      # フロントエンド（単一ファイル）
gunicorn.conf.py          # 常駐時のサーバ設定
Dockerfile / docker-compose.yml
tests/                    # pytest（外部APIはすべてモック）
```

## 常駐時のチューニング

処理は外部APIへの待ち時間が支配的な IO バウンドのため、並行数はプロセスでは
なくスレッドで確保しています（`gthread` ワーカー）。

SIGNOR の索引と HGNC の照合結果は**プロセスごと**にメモリ保持されるため、
`GUNICORN_WORKERS` を増やすとその分だけ複製されます。まず `GUNICORN_THREADS`
を増やし、ワーカー数は必要な場合だけ増やしてください。ディスクキャッシュ
（SIGNOR TSV・STRING・IntAct 相互作用数）はワーカー間で共有されます。

`GUNICORN_TIMEOUT` は既定で 300 秒です。初回リクエストは SIGNOR の一括取得
（最大 120 秒）を含むため、これを大きく下回る値にしないでください。

## テスト

```bash
pip install -r requirements-dev.txt
pytest
```

外部APIはすべてモックしているため、ネットワークなしで実行できます。

## 既知の制限

- 遺伝子シンボルは HGNC で照合しますが、別名が複数の遺伝子を指す場合は
  先頭の候補を採用し `candidates` に候補一覧を返します。
- 結果はキャッシュされないため、同じ問い合わせでも毎回 Open Targets を叩きます
  （PPI・エンリッチメント側はキャッシュされます）。
- 表2は遺伝子ごとに g:Profiler を1回呼ぶため、遺伝子数に比例して解析時間が
  伸びます。不要な場合は STEP 3 のチェックを外してください。
- IntAct の相互作用数が取得できない遺伝子はハブ判定の対象外です
  （取得失敗を「ハブでない」とも「ハブである」とも扱いません）。
