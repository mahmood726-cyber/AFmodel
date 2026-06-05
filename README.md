# AFibSynth — Anticoagulant Network Meta-Analysis (offline)

A single-file, **fully offline** dashboard for atrial-fibrillation stroke-prevention
evidence: it pools log odds ratios across DOAC / VKA / antiplatelet trials with a
multi-chain Gibbs **Bayesian network meta-analysis**, plus a DerSimonian–Laird
pairwise fallback, SUCRA rankings, a comparison-adjusted funnel plot, node-split
inconsistency, dose-response meta-regression, and an individualised net-clinical-benefit
projection driven by CHA₂DS₂-VASc / HAS-BLED inputs.

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network, no external CDN.

## Layout

```
index.html   single-file UI (loads engine.js; builds an MCMC web-worker)
engine.js    pure DETERMINISTIC statistical core — runs in Node and the browser
tests.js     Node test harness, 44 hand-derived assertions
LICENSE      Apache-2.0
```

## Architecture: pure core vs. stochastic sampler

This is a Bayesian app. Its central estimator is a **multi-chain Gibbs MCMC
sampler** (`runBayesianMCMC`, `randomNormal`, `calculateMCMCDiagnostics`,
`calculateSUCRA`) that draws from `Math.random()` — it is *not* pure and stays
inline inside the web-worker. Only the **deterministic** helpers were extracted
into `engine.js` so they have a single source of truth and can be unit-tested
under Node. The worker thread is built by concatenating `engine.js`'s source in
front of the worker code (the same file feeds both the page and the worker, so
there is no duplicated copy).

| Function (engine.js) | What it does |
|---|---|
| `runDLMetaAnalysis(data, kEff, kSe)` | DerSimonian–Laird random-effects pooling of logOR: τ² = max(0,(Q−(k−1))/C), I² = max(0,(Q−(k−1))/Q), pooled on the log scale |
| `runFrequentistNMA(trials, drugs)` | inverse-variance per-drug effect vs. reference, sign-oriented |
| `buildNetworkStructure(trials, drugs)` | network nodes + weighted edges |
| `runSubgroupAnalysis(…, age, renal)` | deterministic class-effect (anti-IIa vs anti-Xa) interaction model |
| `calculateHeterogeneityRisk(trials)` | rule-based I² heuristic (explicitly **not** ML) |
| `checkInconsistency(…)` | direct-vs-indirect node split |
| `assessPublicationBias(trials,…)` | Egger-style precision regression (k≥3 guard) |
| `runDoseResponseAnalysis(trials,…)` | OLS dose meta-regression |
| `generateTimeToEventCurves(…)` | deterministic survival curves |
| `assessQualityFlags(…)` | GRADE-style downgrade flags |
| `calculateNCB(…)` | individualised net clinical benefit |
| `generateLeagueTable(…)` | pairwise OR contrast matrix |

## Fixes applied during revival (2026-06-05)

- **Offline**: removed the Google Fonts `<link>`; the app now loads **no external
  resource** (system fonts fall back). Verified `grep -E 'https?://'` on the
  shipped HTML returns nothing.
- **Single source of truth**: extracted the deterministic statistics into a pure
  `engine.js`. The inline duplicates in the web-worker were deleted; the worker
  is now assembled from the same `engine.js` source, so the page and the worker
  share one copy.
- **Correctness bug fixed** in `assessQualityFlags`: the GRADE score starts at 4
  (HIGH) but the `gradeLabels` array has only four entries (indices 0–3), so a
  network with **zero downgrades** indexed `gradeLabels[4]` and displayed an
  **`undefined`** overall quality grade. The index is now clamped to `[0,3]`, so
  clean evidence correctly reports **HIGH**. Verified by tests.
- **Scaffold**: added `.nojekyll`, `README.md`, `.gitignore`, `E156-PROTOCOL.md`;
  renamed `AFmodel.html` → `index.html`. Dropped the unsupported "Advanced"
  marketing label from the project description.

## Tests

```
node tests.js
# 44 passed, 0 failed
```

Coverage includes a fully hand-worked DerSimonian–Laird case
(τ²≈0.0065794, I²≈20.835%, pooled logOR≈−0.290704, seRE≈0.110580 — derivation in
`tests.js`), the **k=1** edge case (asserts no NaN: τ²=0, I²=0), a
**two-identical-trial** case (τ²=0, I²=0, seRE=√(1/200)), the GRADE off-by-one
regression (clean → HIGH, not `undefined`), the publication-bias k<3 guard, the
sign-flip in the frequentist per-drug effect, and exact dose-regression slope /
intercept.

## Caveats

The reported "MCMC" sampler is a lightweight illustrative Gibbs routine, not a
production NUTS/Stan implementation, and several inputs (subgroup multipliers,
dose map, baseline-risk tables, Egger thresholds) are **hardcoded demonstration
constants**, not fitted from the bundled seven-trial database. DerSimonian–Laird
under-estimates τ² for small *k* (REML / Paule–Mandel preferred for k<10). Treat
all outputs as a transparent teaching / exploration aid, **not** a clinical
decision rule. Apache-2.0 licensed.
