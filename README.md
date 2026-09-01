# NW重複評価アプリ (Gene–Disease Network Overlap Evaluator)

疾患と遺伝子リストを入力すると、遺伝子ごとに PPI ネットワークを構築し、
その近傍が**疾患ネットワーク（Open Targets 上位遺伝子群）を何％カバーするか**を
評価して表にします。

1. **疾患** — Open Targets を検索し、候補から選んで EFO/MONDO ID を確定
2. **遺伝子** — HGNC (HUGO) の公式シンボルと照合（別名・旧シンボルは自動変換）
3. **PPI** — SIGNOR / STRING / BioGRID の組み合わせとスコア条件を選択
4. **解析** — 3つの表で重複を評価
   - **表1 遺伝子重複** — PPIパートナーと疾患遺伝子の重複
   - **表2 パスウェイ重複** — 疾患のエンリッチメント解析結果との重複
   - **表3 症状別マトリックス** — 症状ごとにスコア化し、平均でランキング

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

候補は**最大200件まで取得し、リストはスクロール**できます。ヒット件数と表示件数
を候補リストの下に出すので、絞り込みが必要かどうかが分かります。サブタイプの
多い疾患でも、目的の候補が切り捨てられません。

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
| ハブ判定の相互作用数 | 500 | IntAct のグローバル相互作用数がこれを超えると非特異ハブとして除外 |
| ハブ遺伝子を除外する | ON | 上記に加え、ユビキチン・アクチン・チューブリン等の非特異ハブ族も除外 |
| PPIパートナー上限 | 1000 | ランキング上位から採用する件数。各ソースへの取得件数もこれに追従します |

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

## 3つの表

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

### 表3 — 症状由来遺伝子との重複

疾患の**症状(HPO表現型)**を Open Targets から取得し、各症状に関連する遺伝子を
集めて1つの「症状由来遺伝子セット」を作り、各遺伝子のPPIネットワークがそれを
何％カバーするかを見ます。

```
疾患 → 症状(HPO)リスト → 症状ごとの関連遺伝子 → 統合 → PPIネットワークとの重複
```

HPO の各表現型は Open Targets の疾患インデックスに `HP_0002354` の形で収載されて
いるため、疾患遺伝子と**同じ associations クエリ**でそのまま遺伝子を引けます。
追加のAPIキーやID変換は不要です。

### 症状データの取得元（Open Targets / HPO 直接）

Open Targets 経由には2つの弱点があります。`phenotypes` フィールド自体が返らない
疾患があること、症状→遺伝子の展開が「その HPO term が OT の疾患インデックスにも
収載されていること」に依存することです。そこで **HPO の公開ファイルから直接
取得する経路**を用意し、STEP 3 で切り替えられるようにしています。

| 設定 | 動作 |
|---|---|
| `auto`（既定） | Open Targets を優先し、症状が取れなければ HPO に切り替え。症状ごとの遺伝子も OT で引けなければ HPO の curated リンクを使用 |
| `opentargets` | Open Targets のみ。取れなければ表3は出ません |
| `hpo` | HPO を直接参照 |

HPO 直接経路が使うファイル:

| ファイル | 用途 |
|---|---|
| `phenotype.hpoa` | 疾患 → 症状（頻度・aspect・NOT修飾つき）。OT が取り込んでいる元データそのもの |
| `phenotype_to_genes.txt` | 症状 → 遺伝子（キュレーション済み） |

いずれも初回に一度ダウンロードしてディスク（`PPI_CACHE_DIR`）と メモリに索引化
します（SIGNOR と同じ方式・7日キャッシュ）。列名はヘッダ行から解決するため、
将来のリリースで列順が変わっても値がずれません。

HPO は疾患を **OMIM / ORPHA / DECIPHER** の ID で管理しているため、EFO/MONDO からの
変換が必要です。Open Targets の `dbXRefs` から取得し、取れない場合は疾患名の
完全一致で照合します（自動処理では選択者がいないため、部分一致は行いません）。

### 表3の症状を手動で指定する

STEP 1 の「表3の症状を手動で指定する」から、2通りの指定ができます。

#### A. 症状 (HP ID) を直接指定 — HPO 本来のID

`HP:0002354` のような **HP ID は HPO 自身の表現型ID**です。これを指定すると
疾患を介さずに解析対象の症状が決まり、マトリックスの列が指定したものだけに
なります。名前（`Memory impairment`）でも ID（`HP:0002354` / `HP_0002354`）でも
検索できます。

疾患の症状セットに縛られず、着目したい症状だけで遺伝子を比較したい場合に使います。
疾患指定より優先され、指定した場合は疾患からの症状取得を行いません。

検索対象は**遺伝子が紐づく HP 用語のみ**です（遺伝子のない症状はマトリックスの
列になり得ないため）。

#### B. 疾患 (OMIM / Orphanet ID) を指定

**HPO の疾患注釈は OMIM / Orphanet / DECIPHER の ID で管理されており、
HPO 独自の疾患ID体系はありません。** `phenotype.hpoa` はそれらの疾患を固有の
名前で持っているため、Open Targets を経由せずに疾患を特定できます。

自動変換に頼らず、**HPO 自身の疾患レジストリから疾患を選ぶ**こともできます。
`phenotype.hpoa` は OMIM / Orphanet / DECIPHER の全疾患を固有の名前で持っているため、
Open Targets を経由せずに疾患を特定できます。

**部分一致**で候補を出し、完全一致 → 前方一致 → 部分一致の順、同じ段では注釈の
多い疾患を上位に並べます（選ぶのは人なので、候補を出すのが正解）。

使いどころ:

- Open Targets の `dbXRefs` に OMIM / Orphanet ID がない疾患
- OMIM / Orphanet での登録名が Open Targets と異なる疾患
- より限定的なサブタイプで症状を見たい場合
  （例: `Alzheimer disease` ではなく `Early-onset autosomal dominant Alzheimer disease`）

同じ疾患が OMIM と Orphanet の両方に登録されていることはよくあります。
**両方選べば注釈が統合され**、Orphanet の頻度情報と OMIM の症状の両方が入ります。

指定した場合、表1・表2 は従来どおり STEP 1 で選んだ Open Targets の疾患を使い、
**表3 だけが指定した HPO 疾患のもの**になります。

なお HPO のパスウェイ-遺伝子リンクは**キュレーション済みで重み付けがない**ため、
HPO 由来の症状ではスコアを一律 1.0 とし、加重重複率は単純重複率と一致します。
マトリックスの列見出しに `遺伝子: HPO` と表示されます。

### 症状の出典（Orphanet を含む）

症状注釈は HPO が **OMIM・Orphanet・DECIPHER** を統合したもので、**Orphanet 由来の
症状も含まれます**。複数の情報源が同じ症状を注釈することがあるため出典はすべて
保持し（`resources`）、マトリックスの列見出しと症状リストの下に表示します。

Orphanet の ID 表記は `Orphanet:1020` / `ORPHA:1020` / `ORPHAcode` と揺れるため、
すべて HPO の `ORPHA:` 形式に正規化して照合します。

**頻度注釈（必発／高頻度／頻発／時折／稀）の多くは Orphanet 由来**です。OMIM の
行には頻度が入らないことが多いため、同じ症状に複数の出典がある場合は
Orphanet の頻度を優先します。

`phenotype.hpoa` の aspect が `P`（表現型異常）以外の行 — 遺伝形式(`I`)や
臨床経過(`C`) — は症状ではないので除外します。

結果は**遺伝子 × 症状の数値マトリックス**で表示します。

| | Memory impairment | Dementia | Parkinsonism | … | 総合(平均) |
|---|---|---|---|---|---|
| APP | 100.0 ★ | 88.8 ★ | 41.7 | … | **73.0%** |
| MAPT | 73.7 ★ | 33.5 | 41.7 ★ | … | **49.8%** |
| SNCA | 0.0 | 0.0 | 58.3 ★ | … | **11.7%** |

各セルは、その遺伝子のPPIネットワークが**その症状に関連する遺伝子群**を何％
カバーしているか（加重重複率）です。**総合は全症状の平均値**で、この順に
ランキングします。★ はターゲット遺伝子自身がその症状の関連遺伝子である印です。

症状ごとに分けることで、単一のスコアでは埋もれる**症状特異性**が見えます。
上の例では SNCA は総合 11.7% と低いものの Parkinsonism だけ 58.3% で、
特定の症状に限って強く効く遺伝子だと分かります。

| 指標 | 定義 |
|---|---|
| セル（症状ごと） | その症状の関連遺伝子群に対する加重重複率 |
| 総合 | 全症状のセル値の平均（ランキングに使用） |
| 加重重複率（統合） | 全症状の遺伝子を統合した集合に対する重複率（API・CSVのみ） |
| 単純重複率（統合） | 統合集合に対する遺伝子数ベースの重複率（API・CSVのみ） |

遺伝子のスコアは、**その遺伝子が関与する各症状での関連スコアの合計**です。
強さと広さの両方を反映するので、5つの症状に関わる遺伝子は同じ強さで1症状のみの
遺伝子より重くなります。重複率は合計どうしの比なので [0, 1] に収まります。

**「この疾患には伴わない」と注釈された症状**（HPO の `qualifierNot`）は遺伝子展開の
対象外です。その疾患が起こさない症状の遺伝子を数えることになるためです。全ての
情報源が「伴わない」としている症状のみを除外し、情報源によって判断が割れる場合は
残します。

表1と同様、ターゲット遺伝子自身が症状由来遺伝子である場合はそれを含めて数えます。

症状注釈は疾患によって粗密があります。一般に希少・単一遺伝子疾患は充実し、
common disease では疎です。注釈がない疾患では表3は表示されず、表1・表2は
通常どおり動作します。

### 幅広く症状に関連する遺伝子の有意性

症状ごとのスコアだけでは「1症状に強く当たる遺伝子」と「多くの症状に一貫して
当たる遺伝子」を区別できません。後者を統計的に評価します。

1. **各症状で超幾何検定** — その遺伝子のPPIネットワークとその症状の遺伝子群の
   重複が、同じ大きさの近傍を無作為に取った場合より多いか
2. **Fisher の方法で統合** — 全症状の p 値をまとめる。中程度の重複が多くの症状に
   わたる方が、1症状だけ強く当たるより有意になる
3. **Benjamini-Hochberg 補正** — 入力遺伝子すべてが同じ症状集合に対して検定される
   ため、遺伝子間で FDR 補正

| 列 | 意味 |
|---|---|
| カバー症状数 | 重複が1件以上あった症状の数と割合 |
| q値 (FDR) | 上記の統合 p 値を BH 補正したもの。`< 0.05` で「有意」表示 |

実例（平均がほぼ同じでも、広さで差がつきます）:

| 遺伝子 | カバー | q値 | 総合(平均) | |
|---|---|---|---|---|
| APOE | 2/5 | 1.5e-03 | 11.8% | **有意** |
| SNCA | 1/5 | 6.0e-02 | 11.7% | 非有意 |

背景遺伝子数は既定 20,000（ヒトのタンパク質コード遺伝子数）で、
`NW_SYMPTOM_BACKGROUND` で変更できます。表は「総合（平均）／カバー症状数／
有意性(q値)」で並び替えられます。

SciPy は使わず標準ライブラリのみで実装しています（超幾何分布は対数空間、
Fisher の統合は自由度が常に偶数なので閉形式）。コンテナを軽く保つためです。

## 表のコピー（PowerPoint / Excel）

各表の上に2つのボタンがあります。

| ボタン | 用途 |
|---|---|
| 表をコピー（PowerPoint / Word） | クリップボードの `text/html` に書き出します。PowerPoint に貼ると**色つきの編集可能な表**になります |
| TSVでコピー（Excel） | タブ区切りテキスト。Excel にそのまま貼れます |

PowerPoint と Word は `<style>` ブロックを無視するため、**計算済みスタイルを
各セルの `style` 属性に展開**しています。あわせて:

- **クリップボードへの書き込みは `copy` イベント経由**です。
  `navigator.clipboard` は secure context（HTTPS または localhost）でしか
  使えず、Docker のアプリを LAN の IP やホスト名で開くと存在しません。
  `copy` イベントなら平文 HTTP でも動き、こちらで生成した HTML が
  そのままクリップボードに載ります
- セルの塗りつぶしは inline CSS に加えて **`bgcolor` 属性でも出力**します。
  **貼り付け先は `style` 属性を落とすことがあり、その場合 `bgcolor` だけが
  残ります** — inline CSS だけでは色が完全に消えます
- `color(srgb ...)`（CSS Color 4）は Office が解釈できないため **16進に変換**
- ヒートマップの半透明色は**白に合成して単色化**（Office はセルのアルファを扱えない）
- クラス名は除去したうえで、`display` などレイアウトに効くプロパティも
  インライン化します（除去しただけだと見出しの改行が潰れるため）
- ダークモードで表示していても、**コピーは常にライトテーマ**で生成されます
  （スライドに貼る前提のため）
- 進捗バーなどの装飾要素は除去

**貼り付け時の注意:** PowerPoint で貼り付けた直後に出る
［貼り付けのオプション］で **［元の書式を保持］** を選んでください。
［貼り付け先のテーマを使用］だとセルの塗りつぶしが外れます。
この案内はコピーボタンの下にも表示しています。

画像として貼りたい場合は、PowerPoint 側で「形式を選択して貼り付け → 図」を
選んでください。

## ターゲット遺伝子自身の扱い

両方の表で、ターゲット遺伝子自身の寄与を重複に含めます。

| | ターゲット自身の数え方 |
|---|---|
| 表1 | ターゲットが疾患遺伝子リストに含まれる場合、そのOTスコアを分子に加算し `matched_count` にも計上 |
| 表2 | ターゲットが疾患パスウェイの構成遺伝子である場合、そのパスウェイを重複として計上 |
| 表3 | 症状ごとのセルでも統合スコアでも、ターゲットが関連遺伝子である場合はそれを含めて計上（マトリックスでは ★ 表示） |

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

### `GET /api/diseases?q=<疾患名>&limit=50`

`limit` は既定 50、最大 200。Open Targets へのページサイズとしても渡すので、
大きくすれば実際に候補が増えます。`total` はヒット総数です。

```jsonc
{ "query": "alzheimer disease", "total": 42,
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
  "symptoms": true,                    // 表3を計算するか（既定 true）
  "symptom_source": "auto",            // auto / opentargets / hpo
  "hpo_phenotype_ids": ["HP:0002354"], // 症状を直接指定（HP ID、最優先）
  "hpo_disease_ids": ["ORPHA:1848"],   // 疾患を指定（OMIM/Orphanet ID、最大10件）
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
  "symptoms": {
    "enabled": true,
    "requested_source": "auto", "phenotype_source": "hpo",
    "gene_sources": ["hpo"], "xrefs": ["OMIM:104300", "Orphanet:1020"],
    "xref_origin": "dbxrefs",          // phenotypes / selected / dbxrefs / name
    "phenotype_count": 5, "expanded_count": 3,
    "gene_count": 6, "excluded_count": 1, "unindexed_count": 1,
    "phenotypes": [ { "hpo_id": "HP:0002354", "name": "Memory impairment",
                      "frequency": "HP:0040281",
                      "resources": ["ORPHANET", "OMIM"], "excluded": false } ],
    "matrix_symptoms": [ { "hpo_id": "HP:0002354", "name": "Memory impairment",
                           "frequency": "HP:0040281",
                           "resources": ["ORPHANET", "OMIM"], "gene_count": 3 } ],
    "expanded":   [ { "hpo_id": "HP:0002354", "name": "Memory impairment",
                      "frequency": "HP:0040281", "gene_count": 3 } ],
    "top_genes":  [ { "symbol": "APP", "score": 1.37, "phenotype_count": 2,
                      "phenotypes": ["Memory impairment", "Dementia"] } ]
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
      ],

      // 表3（symptoms 有効時かつ症状注釈がある場合のみ）
      "symptom_weighted_percent": 79.5, "symptom_overlap_percent": 66.7,
      "symptom_matched_count": 4, "symptom_overlap_count": 3,
      "symptom_gene_count": 6,
      "symptom_target_self": "APP", "symptom_target_self_score": 1.37,
      "symptom_interpretation": "strong",
      "symptom_overlapping_genes": [ { "symbol": "MAPT", "score": 1.1 } ],
      "symptom_mean_percent": 73.0,            // 総合（ランキングに使用）
      "symptom_breadth_count": 5, "symptom_tested_count": 5,
      "symptom_breadth_percent": 100.0,
      "symptom_p_value": 5.53e-30,             // Fisher 統合
      "symptom_q_value": 2.76e-29,             // BH 補正
      "symptom_cells": [                        // 症状ごとのスコア
        { "hpo_id": "HP:0002354", "name": "Memory impairment", "percent": 100.0,
          "matched_count": 3, "gene_count": 3, "target_self": true }
      ]
    }
  ]
}
```

`results` は加重重複率の降順です。個別の遺伝子で PPI 取得に失敗した場合、
その要素は `{"gene": ..., "error": ...}` になり、他の遺伝子の結果は返ります。

ステータスコード: `400`（入力不備・遺伝子数超過）、`404`（疾患が見つからない・
関連遺伝子ゼロ）、`502`（Open Targets 側の障害）。

### `GET /api/hpo/phenotypes?q=<症状名 または HP ID>&limit=50`

`limit` は既定 50、最大 500。`total` は絞り込み前のヒット総数です。

HPO の表現型 (HP 用語) を検索します。遺伝子が紐づく用語のみが対象です。

```jsonc
{ "query": "HP:0002354",
  "results": [
    { "hpo_id": "HP:0002354", "ontology_id": "HP_0002354",
      "name": "Memory impairment", "gene_count": 12, "exact": true }
  ] }
```

### `GET /api/hpo/diseases?q=<疾患名>&limit=50`

HPO 自身の疾患レジストリを検索します。`limit` は既定 50、最大 500。

```jsonc
{ "query": "alzheimer", "total": 33,
  "results": [
    { "id": "ORPHA:1020", "name": "Alzheimer disease", "source": "ORPHANET",
      "phenotype_count": 3, "exact": false }
  ] }
```

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
| `NW_HUB_THRESHOLD` | `500` | ハブ判定の既定値 |
| `NW_MAX_NODES` | `1000` | PPIパートナー上限の既定値 |
| `NW_ENRICH_GENE_N` | `100` | エンリッチメントに使う疾患遺伝子の件数 |
| `NW_ENRICH_MAX_TERMS` | `50` | 取得するパスウェイ数の上限 |
| `NW_SYMPTOM_LIST_N` | `50` | 取得する症状(HPO)の件数 |
| `NW_SYMPTOM_MAX` | `20` | 遺伝子に展開する症状の件数（1件につきAPI呼び出し1回） |
| `NW_SYMPTOM_GENES_PER` | `50` | 症状1件あたりに取得する遺伝子数 |
| `NW_SYMPTOM_SOURCE` | `auto` | 症状の取得元（`auto` / `opentargets` / `hpo`） |
| `NW_SYMPTOM_BACKGROUND` | `20000` | 有意性検定の背景遺伝子数 |
| `PPI_CACHE_DIR` | `./ppi_cache` | PPI キャッシュの保存先 |
| `PPI_CACHE_DISABLED` | 無効 | `1` でディスクキャッシュを無効化 |
| `GUNICORN_WORKERS` / `GUNICORN_THREADS` | `2` / `8` | 常駐時のサーバ規模 |
| `GUNICORN_TIMEOUT` | `300` | リクエストタイムアウト（秒） |
| `LOG_LEVEL` | `info` | ログレベル |

外部APIのURL（`OT_API_URL` / `HGNC_API_URL` / `SIGNOR_TSV_URL` / `STRING_URL` /
`BIOGRID_URL` / `INTACT_URL` / `GPROFILER_URL` / `HPO_HPOA_URL` /
`HPO_PHENOTYPE_TO_GENES_URL`）も上書きできます。ミラーの利用やテストに使います。

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
| [HPO](https://hpo.jax.org/) | 症状（疾患→表現型→遺伝子）の直接取得 | 不要 | HPO ライセンス |

各ソースは best-effort です。1つが落ちても残りで解析を続行し、
HGNC に到達できない場合も入力シンボルのまま解析を進めます。

## 構成

```
app.py                    # Flask ルーティング・入力検証・並列実行
nw_overlap.py             # 表1: 遺伝子レベルの重複スコア
enrichment_overlap.py     # 表2: パスウェイレベルの重複スコア
symptom_genes.py          # 表3: 症状 → 遺伝子セットの構築
symptom_stats.py          # 表3: 超幾何検定・Fisher統合・BH補正
ppi_network.py            # PPIグラフ構築・パートナーランキング・ハブ判定
collectors/
├── _http.py              # リトライ付き共有 HTTP セッション
├── _cache.py             # ディスクキャッシュの場所
├── opentargets.py        # 疾患検索・疾患上位遺伝子
├── hgnc.py               # HGNC シンボル照合
├── hpo.py                # HPO 直接取得（症状・症状別遺伝子）
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
- 表3は症状1件につき Open Targets を1回呼びます（既定で最大20件、リクエストごと
  に1回だけ・全遺伝子で共有）。症状注釈のない疾患では表3は出ません。
- PPIパートナー上限を大きくすると、ハブ判定のための IntAct 問い合わせが
  パートナー数だけ発生します（30日キャッシュされ遺伝子間で共有されますが、
  初回は時間がかかります）。速度を優先する場合はハブ除外をオフにしてください。
- 有意性検定は各症状の遺伝子集合を互いに独立とみなしています。実際には症状間で
  遺伝子が重複するため、p 値はやや小さめに出ます。順位付けの指標として扱い、
  絶対値を過信しないでください。
- Open Targets が疾患として収載していない HPO 表現型は、`auto` では HPO の
  curated リンクで補完します。`opentargets` 固定にした場合はスキップされ、
  `unindexed_count` に計上されます。
- 自動の HPO 照合は疾患名の**完全一致**のみです。`dbXRefs` に OMIM / Orphanet の
  ID がなく名前も一致しない場合は、STEP 1 から疾患 (OMIM / Orphanet ID) を選ぶか、
  症状 (HP ID) を直接指定してください（どちらも部分一致で検索できます）。
- IntAct の相互作用数が取得できない遺伝子はハブ判定の対象外です
  （取得失敗を「ハブでない」とも「ハブである」とも扱いません）。
