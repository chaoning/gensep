# gensep — 方案 (plan)

自包含 C++ 工具,把 LDAK SumHer 的 `--sum-hers` / `--sum-cors` 估计 + block-jackknife
**全移植**过来,并在同一批 jackknife 块上派生 **VS / case-case h²(仅 formula 路径)/ upper-limit AUC**,
全部给出 **block-jackknife SE**。

落地决策(2026-05-24 与用户确认):
- 路线 A 的"全移植"版:不依赖 LDAK 二进制,tagging/summaries 解析 + 求解器 + jackknife 全在 gensep。
- 线性代数用 **Eigen(header-only)**;矩阵是 total×total(total≤3),无需 MKL/BLAS。
- 派生量的**点估计保持现状**:h₁/h₂ 用 sum-hers、rg 用 sum-cors;**只有 SE** 走 sum-cors 联合 jackknife。
- `--num-blocks 200` 固定。
- **case-case h² 仅算 formula 路径**(VS/(VS+4));gdis / direct 不在 gensep 范围 →
  不需要 cc 的 `.hers.jackknife`、不需要 50/50 缩放、不需要 `--Pcc`。
  (若仍要 fig2 的 gdis-vs-direct 验证,那部分留在 step5 用现有数,与 gensep 解耦。)

---

## 0. 已探明的关键事实

- **HumDef.tagging = 单类别**。表头 `Predictor A1 A2 Neighbours Tagging Weight MAF Categories Exp_Heritability Base`;
  559763 SNP 行 + 2 行 footer。`count = countrows - 2 - num_parts`,此处 num_parts=1。
  - `stags[j]` = 第5列 Tagging;`svars[0][j]` = 第10列 Base;`ssums[0]` = footer 行的 559763。
- 求解器 `total = num_parts + gcon + cept ≤ 3` → 正规方程极小,大循环只在 ~5.6e5 SNP 上做点积。
- LDAK jackknife = delete-one-block;SE 公式两处一致:
  `var = (B-1)*(mean(θ²) - mean(θ)²)`,即 `SE = √[(B-1)/B · Σ(θ_b - θ̄)²]`
  (对账 `sumfuns.c:947` 与 `sumfuns.c:1570`)。gensep 逐字复用。
- 当前流水线用预编译 `/home/chaoning/software/LDAK/ldak`(需 MKL),源码在 `code/LDAK/`。
  gensep 走 Eigen → **绕开 MKL 编译难题**,这是 C++ 重写相对打补丁的主要收益。

## 1. 要移植的 LDAK 函数(源在 code/LDAK/)

| LDAK 源 | 作用 | gensep 文件 |
|---|---|---|
| `read_tagfile` (parsefiles.c) | 解析 .tagging → stags / svars / ssums / preds / al1 al2 | `tagging.{hpp,cpp}` |
| `read_sumsfile` (parsefiles.c) | 解析 .summaries → snss / schis / srhos(按 pred+allele 对齐 tagging) | `summary.{hpp,cpp}` |
| `sumsa.c` QC 段 (L58–160) | 去 N=0、截断大 χ²、stags<1→1、schis==0→1e-6、按 SNP 对齐压缩 | `qc.{hpp,cpp}` |
| `solve_sums` (sumfuns.c:245) | 单 trait WLS+GC 迭代 + jackknife | `sumhers.{hpp,cpp}` |
| `solve_cors` (sumfuns.c:1047) | pair 联合 Her1/Her2/Coher/Cor + jackknife | `sumcors.{hpp,cpp}` |
| `get_factor` (regfuns.c:41) + `code/common/case_case_auc.py` | Lee / λ / VS / h²cc / AUC | `gensep.{hpp,cpp}` |
| `ldak.c` mode 146/147 调度 | CLI 装配 | `main.cpp` |

## 2. 派生量公式(复用 code/common/case_case_auc.py)

λᵢ = φ(Φ⁻¹(1−Kᵢ))/Kᵢ,tᵢ = Φ⁻¹(1−Kᵢ),δᵢ = λᵢ(λᵢ−tᵢ)。h₁/h₂ 为 liability-scale。

- **VS** = λ₁²h₁ + λ₂²h₂ − 2λ₁λ₂·rg·√(h₁h₂)
- **h²cc · formula** = VS/(VS+4)   ← gensep 唯一的 case-case h² 输出
- **AUC · exact** = Φ( VS / √(2VS − a₁²δ₁ − a₂²δ₂) ),aᵢ = λᵢhᵢ − λ_{3−i}·rg·√(h₁h₂)
- **AUC · approx** = Φ(√(VS/2))

（已移除:h²cc · gdis = H₁₅₀+H₂₅₀−2rg√(H₁₅₀H₂₅₀)、h²cc · direct = h_cc_obs·0.25/(P_cc(1−P_cc))。）

## 3. K / P 输入

- **observed h²/rg/coher + SE**:不需要 K、不需要 P。
- **VS / h²cc(formula) / AUC**:需要 **K₁,K₂**(λ、Lee)。
- **Lee obs→liab** `(K(1−K)/φ(t))²/(P(1−P))`:需要 **P**(样本病例比例)。formula 路径不含 50/50 → 不需要 Pcc。
- ⚠️ P 不在 summaries 里(只有每 SNP 有效 N),必须 CLI 传入:`--P1 --P2`。
- 接口:K/P 全可选。不给 → 只出 observed + SE;给齐 → 额外出 VS/h²cc_formula/AUC + SE。

## 4. CLI

```
gensep hers --summary g.summaries --tagfile HumDef.tagging --out PREFIX
   → PREFIX.hers (Component Heritability SE …)  ← SE 来自 likelihood Hessian(LDAK 默认 chisol=1),
     不做 jackknife、不产 .hers.jackknife:case-case h² 仅 formula 路径,h1/h2 只作点估计,
     其 SE 与 LDAK .hers 完全一致即可;所有派生量 SE 一律走 sum-cors 的 jackknife。

gensep cors --summary g1.summaries --summary2 g2.summaries --tagfile T \
             --num-blocks 200 --out PREFIX
   → PREFIX.cors, PREFIX.cors.jackknife (B 行: Her1 Her2 Coher Cor)

gensep sep  --cors-jack PREFIX.cors.jackknife \
             --h1obs .. --h2obs .. --rg ..       # 点估计(现状来源)
             --K1 .. --K2 .. --P1 .. --P2 .. --out PREFIX
   → PREFIX.gensep: VS h2cc_formula auc_exact auc_approx + 各 SE
```
（`sep` 也可直接并入 `cors`，由是否提供 K/P 决定是否输出 .gensep。）

## 5. block-jackknife + SE 传播(sep 子命令)

读 `PREFIX.cors.jackknife`(B 行 `[h1_obs_b, h2_obs_b, coher_b, rg_b]`),逐块:
1. `h1_liab_b = h1_obs_b·lee(K1,P1)`，`h2_liab_b = h2_obs_b·lee(K2,P2)`
2. `rg_b' = clip(rg_b, ±0.999)`
3. 代入 §2 得 `VS_b, h2cc_formula_b, auc_exact_b, auc_approx_b`
4. 无效块处理:`VS_b<0` 或 `denom_sq<0` → 该量该块 NaN,成对删除,记有效块数。

SE = √[(B−1)/B · Σ(θ_b − θ̄_b)²]，θ̄_b 用块均值(与 LDAK 一致)。
⚠️ 报告的**点估计**仍用 §对应现状来源(sum-hers h₁/h₂ + sum-cors rg)，
故 `VS_point ≠ mean(VS_b)`;SE 是离散度,与点估计并列报告,不取代之。
T1D 等 `h_liab>1` 不 clip(沿用既有决定),SE 照算。

**NaN/越界块处理(两套,刻意不同)**:
- **sum-cors 输出(.cors)= LDAK 忠实移植**:用 `jackknife_se_ldak`,对全部 200 块算,
  任一块 NaN(如某块 Her1·Her2<0 → Cor_b=NaN)→ SE=NaN,与 LDAK 的 `-nan` 完全一致。
- **sep 派生量(.gensep)= 新功能,刻意成对删除**:某块 `VS_b≤0`(VS 是方差;h²cc 近极点;AUC 无定义)
  → 该块 VS/h²cc/AUC **全部**记 NaN 丢弃;`auc_exact` 另在 `denom≤0` 时丢弃。
  点估计不受影响(VS<0 也照报)。报告 `B_used(VS,h2cc)` 与 `B_used(auc_exact)` 两个有效块数。

## 6. 数值保真(两大风险点,先验证)

1. **read_tagfile**:tagging 文件列/footer 解析、SNP+allele 对齐(含 amb/翻转)必须和 LDAK 一致。
2. **solve_sums/solve_cors 的 GC + scaling 迭代**:逐行照搬迭代(权重随拟合更新、gc、scale、intercept)。

**回归测试**:磁盘已有 LDAK 输出可对账,例如
`classification/simdata/end_to_end_subControl_1to1/.../rep_8/ldak_rg.cors`、`ldak_h2_cc.hers`。
gensep 在相同 summaries+tagging 上跑,核对 h²/rg/coher 及 SE,目标 |Δ|<1e-3 再切流水线。

## 7. 集成

- `step3_ldak.py`:`LDAK --sum-hers/--sum-cors` → `gensep hers/cors`(命令等价,加 `--num-blocks 200`)。
- `step5_metrics.py`:VS/h²cc_formula/AUC 的**点估计照旧自己算**,新增直接读 `.gensep` 的 **SE 字段**写进 result.json
  (`VS_se / h2cc_formula_se / auc_exact_se / auc_approx_se`)。gdis/direct 若保留则仍是无 SE 的点估计。
- `collect_results.py`:fig1/2/4 加误差棒。

## 8. 构建

- C++17,header-only Eigen(放 `code/gensep/third_party/eigen/` 或系统路径)。
- `Makefile`:`g++ -O3 -std=c++17 -I third_party/eigen *.cpp -o gensep`(无外部链接)。
- 可选 OpenMP 给 jackknife/SNP 循环加速(`-fopenmp`)。

## 9. 实施顺序(分阶段,每阶段可单测)

1. 目录 + Makefile + Eigen + 本 plan(本次)。
2. `tagging.cpp` + `summary.cpp` + `qc.cpp`:解析对齐,对一个 summaries 打印 stags/snss/schis 前几行核对 LDAK 日志。
3. `sumhers.cpp`(solve_sums + jackknife dump):对账 `.hers` 点估计与 SE。
4. `sumcors.cpp`(solve_cors + jackknife dump):对账 `.cors`。
5. `gensep.cpp`(sep):VS/h²cc/AUC + SE。
6. 切 step3/step5/collect,7→35 对回归跑。
