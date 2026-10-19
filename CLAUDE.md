# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 概要

CSS 符号で横断的対角ゲート $U(\theta)=\bigotimes_j \mathrm{diag}(1,\omega^{\theta_j})$（$\omega=e^{2\pi i/2^m}$、角度は非一様可）が真にレベル $m$ の論理ゲートを与えるかを判定する、arXiv:2602.14499（Haruna）の判定原理の最小リファレンス実装。兄弟リポジトリ homology_and_transversality の判定原理を切り出した位置づけ。純 Python・外部依存なし・パッケージ化なし（フラットな 4 ファイル）。判定原理と結果表は `README.md`（英語）と `README_ja.md`（日本語、内容は同じ）に詳しい。README を変えるときは両方を同期する。

## コマンド

Windows では `python3` ではなく `python` を使う。

```
python codes.py        # 符号の自己検査のみ（H_X H_Z^T = 0、[[n,k,d]] を全列挙で確認）
python test_known.py   # 判定 + 既知事実テスト。結果表と PASS/FAIL を出力、FAIL があれば exit 1
```

テストは pytest ではなく `test_known.py` の `main()` 内の `report(...)` 呼び出しの列（全体で 3〜4 秒程度で、個別実行の仕組みはない）。単発で調べるなら:

```
python -c "from bockstein_check import analyze; from codes import steane; print(analyze(*steane(), 3))"
```

## アーキテクチャ

- `bockstein_check.py` — 判定本体。流れは `analyze(HX, HZ, m)`:
  1. `constraint_rows`: $\mathrm{rowspan}(H_X)$ の**全 $2^{r_X}$ 元** × 全 $2^k$ 論理類代表について、制約行 $s_j(1-2l_j)$（成分 $\in\{-1,0,1\}$）を作る。整数持ち上げ $a\oplus b=a+b-2ab$ が GF(2) 上非線形なので、生成元だけで張ってはいけない（この設計の要）。
  2. `solve_homogeneous_mod`: `diagonalize_int`（Smith 標準形型の整数対角化、$V$ だけ追跡）で $A\theta\equiv0 \pmod{2^m}$ の解モジュールの生成元を厳密に求める。$\theta$ の総当たりはしない。
  3. 各非自明論理類で、生成元上の論理位相 $\sum_j\theta_j l_j$ の gcd $d$ を取る。達成可能位相は $d\cdot\mathbb{Z}_{2^m}$、`feasible` ⟺ $d=1$（奇数位相が出る）。可なら奇数位相を 1 に正規化した例 `example` を返す。
  - `verify_theta` は線形代数を使わない総当たり位相列挙による**独立検証**。`analyze` が可と判定した例は必ずこれで検証する（`test_known.py` がそうしている）ので、ソルバを変えてもこの二重チェックは残す。
- `generator_form.py` — 同じ判定の生成元形式。制約を $h,\ 2h{*}k,\ 4h{*}k{*}k',\dots$（$h$ は $C_2$ 基底、$k$ は $\ker H_Z$ 基底、$*$ は成分積）で張り、最上位レベルは生成中に GF(2) 簡約、残りは $\mathbb{Z}_{2^m}$ 上の逐次核計算（`solve_levels`）で解く。ベクトルは int ビットマスク、$\mathbb{Z}_{2^m}$ ベクトルは $m$ 枚のビットプレーン。`analyze` と同じ dict を返す。大きい符号の例は `verify_theta_sampled`（抜き取り）で検証する。CLI: `python generator_form.py HX.txt HZ.txt [m]`。`test_known.py` が両形式の一致（3 符号と seed 固定のランダム小符号、$m=1..4$）を確認しているので、どちらかを変えたらここが壊れないことを確かめる。
- `codes.py` — 符号コンストラクタ（`steane`, `reed_muller_15`, `rotated_surface_9`、いずれも `(HX, HZ)` を GF(2) の 0/1 リストで返す）と `check_code`。新しい符号は `ALL_CODES` に `(名前, コンストラクタ, (n,k,d))` で足す。大きい族 `rotated_surface(d)`, `reed_muller(r)` もある（距離の全列挙が重いので `ALL_CODES` には入れない）。`test_known.py` は `codes[0..2]` を位置で Steane / RM15 / 表面符号として参照している。
- `test_known.py` — 各符号を $m=1,2,3$ で判定し、既知事実（Steane: S 可・T 不可、RM15: T 可、表面符号: S・T 不可）と照合。

## 注意

- `analyze` は全列挙なので計算量は $2^{r_X+k}$ に比例し、小さい符号専用（参照実装として残す）。大きい符号は `generator_form.analyze_generators` を使う（積の数は $O(r_X(\dim C_1)^{m-1})$、$n=185,\ r_X=90,\ m=3$ で数秒）。$m\ge4$ は積の数が急増する。
- 対象は単一量子ビット対角ゲートのテンソル積のみ（CZ 等の多体対角は対象外）。判定基準「ある非自明類で奇数位相」は README の "Implementation notes"（README_ja.md では「実装上の注意」）に挙げた解釈上の仮定に基づく。判定結果や仮定を変えたら README の結果表・注意も合わせて更新する。
- `test_known.py` の出力は ASCII のみに保つ（docstring に明記された方針）。
