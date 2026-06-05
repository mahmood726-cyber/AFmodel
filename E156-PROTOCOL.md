# E156-PROTOCOL — AFibSynth (Anticoagulant Network Meta-Analysis)

- **Project:** AFmodel (GitHub repo `AFmodel`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `AFmodel.html` dump)
- **Type:** single-file offline browser tool + Node-testable deterministic engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: removed the Google Fonts CDN `<link>`; the app now
  loads no external resource (system fonts fall back). No other external CDN was
  present.
- Extracted the **deterministic** statistics into a pure `engine.js` (single
  source of truth; the inline duplicates were deleted and the MCMC web-worker is
  now assembled from the same `engine.js` source). The stochastic Gibbs sampler
  (`Math.random`-based) stays inline, as it is not pure.
- Added `tests.js` (44 hand-derived assertions, all passing).
- **Fixed a correctness bug** in `assessQualityFlags`: the GRADE score starts at
  4 (HIGH) but `gradeLabels` has only four entries, so a network with zero
  downgrades displayed an `undefined` quality grade — the index is now clamped to
  `[0,3]` (clean evidence → HIGH).
- Added Pages scaffold (`.nojekyll`, README, `.gitignore`); renamed
  `AFmodel.html` → `index.html`. Dropped the unsupported "Advanced" label.

## Body (E156 draft — CURRENT BODY)

Which oral anticoagulant minimises the joint stroke-versus-bleed burden for a
given atrial-fibrillation patient, and how trustworthy is that ranking? This
browser tool bundles seven landmark DOAC, warfarin and antiplatelet trials and
lets a user set a CHA₂DS₂-VASc and HAS-BLED profile to individualise the
contrast. It runs a multi-chain Gibbs Bayesian network meta-analysis of log odds
ratios with a DerSimonian–Laird pairwise fallback, then layers SUCRA rankings,
node-split inconsistency, a comparison-adjusted funnel plot and a net-clinical-benefit
projection. Across the bundled network the apixaban-class agents rank highest on
net benefit, but the credible intervals overlap and the prediction signal is
heterogeneity-dominated, so no single drug is unambiguously best. A revival audit
fixed an off-by-one that displayed an undefined GRADE for clean evidence and
locked the deterministic core behind a 44-assertion test suite. The honest read
is that anticoagulant choice is profile-dependent shared decision-making, not a
fixed ranking. The tool is a transparent exploration aid, not a clinical
decision rule.

SUBMITTED: [ ]
