# NW重複評価アプリ (Gene–Disease Network Overlap Evaluator)

疾患名と遺伝子リストを入力すると、遺伝子ごとに PPI パートナーを取得し、
その近傍が疾患の Open Targets 上位遺伝子群をどれだけカバーするかを
加重スコアで評価して表にします。

- **インプット**: 疾患名（または EFO/MONDO ID） + 遺伝子リスト
- **解析**: 遺伝子ごとに PPI パートナーを収集 → 疾患上位遺伝子との重複を加重スコア化
- **アウトプット**: 遺伝子ごとのスコア・重複遺伝子リストのテーブル（CSV 出力可）

## セットアップ

```bash
pip install -r requirements.txt
python app.py
```

http://127.0.0.1:5005 を開きます。

## 設定（環境変数）

| 変数 | 既定値 | 説明 |
|---|---|---|
| `PORT` | `5005` | 待ち受けポート |
| `HOST` | `127.0.0.1` | 待ち受けアドレス |
| `BIOGRID_KEY` | （空） | BioGRID API キー。未設定なら SIGNOR + STRING のみで動作 |
| `NW_MAX_GENES` | `30` | 1リクエストあたりの遺伝子数上限 |
| `NW_DISEASE_TOP_N` | `100` | 疾患上位遺伝子の取得件数 |
| `NW_PPI_TOP_N` | `30` | 遺伝子あたりの PPI パートナー上限 |
| `NW_MAX_WORKERS` | `5` | 並列ワーカー数 |
| `NW_ENABLE_SIGNOR` / `NW_ENABLE_STRING` | 有効 | 個別のデータソースを無効化する場合に `0` |
| `FLASK_DEBUG` | 無効 | `1` でデバッグモード |

## スコアの定義

```
weighted_score = ( Σ 重複遺伝子のOTスコア + ターゲット自身のOTスコア )
                 / Σ 全疾患遺伝子のOTスコア
```

ターゲット自身が疾患遺伝子リストに含まれる場合、そのスコアが分子に加算されます
（重複遺伝子側では二重計上されません）。

| 加重スコア | 解釈 |
|---|---|
| ≥ 0.3 | ネットワーク的に疾患と強く関連 |
| 0.1〜0.3 | 中程度の関連 |
| < 0.1 | 関連が薄い（または PPI 情報不足） |

`simple_ratio` は重複遺伝子数を疾患遺伝子数で割った、スコアで重み付けしない指標です。

## API

### `POST /api/analyze`

```jsonc
// リクエスト
{ "disease": "Alzheimer disease", "genes": ["APP", "PSEN1", "APOE"], "top_n": 100 }
```

`genes` は配列でも、改行・カンマ区切りの文字列でも受け付けます。`top_n` は任意
（省略時は `NW_DISEASE_TOP_N`、上限 500）。

```jsonc
// レスポンス
{
  "disease_id": "EFO_0000249",
  "disease_label": "Alzheimer disease",
  "disease_gene_count": 100,
  "ppi_top_n": 30,
  "results": [
    {
      "gene": "APP",
      "weighted_score": 0.81,
      "simple_ratio": 0.625,
      "overlap_count": 4,
      "disease_gene_count": 100,
      "ppi_partner_count": 6,
      "target_self": "APP",
      "target_self_score": 0.92,
      "interpretation": "strong",
      "overlapping_genes": [{ "symbol": "PSEN1", "score": 0.71, "target_id": "ENSG..." }]
    }
  ]
}
```

`results` は `weighted_score` の降順です。個別の遺伝子で PPI 取得に失敗した場合、
その要素は `{"gene": ..., "error": ...}` になり、他の遺伝子の結果は返ります。

ステータスコード: `400`（入力不備・遺伝子数超過）、`404`（疾患が見つからない・
関連遺伝子ゼロ）、`502`（Open Targets 側の障害）。

### `GET /healthz`

現在の設定値と稼働状態を返します。

## データソース

| ソース | 用途 | APIキー |
|---|---|---|
| [Open Targets Platform](https://platform.opentargets.org/) (GraphQL v4) | 疾患ID解決・疾患上位遺伝子 | 不要 |
| [SIGNOR](https://signor.uniroma2.it/) | 制御性シグナル相互作用 | 不要 |
| [STRING](https://string-db.org/) | 機能的関連ネットワーク | 不要 |
| [BioGRID](https://thebiogrid.org/) | キュレーション済み相互作用 | 必要 |

各ソースは best-effort です。1つが落ちても残りで解析は続行し、
どのソースからもパートナーが得られなかった遺伝子はスコア 0 になります。

SIGNOR は遺伝子単位のエンドポイントを持たないため、ヒト相互作用データを
プロセスごとに1回だけ一括取得してメモリ上に索引化します（毎リクエスト・
毎遺伝子での再取得を避けるため）。

PPI パートナーが `NW_PPI_TOP_N` を超える場合は、**何ソースに支持されているか** →
**STRING の combined score** → アルファベット順、の優先度で上位を採用します。

## 構成

```
app.py                    # Flask アプリ本体・入力バリデーション・並列実行
nw_overlap.py             # NW重複スコア計算
collectors/
├── _http.py              # リトライ付き共有 HTTP セッション
├── opentargets.py        # 疾患ID解決・疾患上位遺伝子取得
└── ppi.py                # SIGNOR / STRING / BioGRID からの PPI 収集
templates/index.html      # フロントエンド（単一ファイル）
tests/                    # pytest（外部APIはモック）
```

## テスト

```bash
pip install -r requirements-dev.txt
pytest
```

外部APIは全てモックしているため、ネットワークなしで実行できます。

## 既知の制限

- ハブ遺伝子の除外（degree 閾値によるフィルタ）は行っていません。相互作用数が
  極端に多い遺伝子はスコアが高めに出ることがあります。
- 遺伝子シンボルはエイリアス解決を行わず、大文字化して照合します。
  Open Targets の `approvedSymbol` 表記に揃えるのが確実です。
- 結果はキャッシュされないため、同じ問い合わせでも毎回外部APIを叩きます。
