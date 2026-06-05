/*
 * Node tests for the AFibSynth deterministic engine. Run: node tests.js
 *
 * Every expected value is INDEPENDENTLY hand-computed (see comments / README),
 * NOT a re-run of the engine. The MCMC sampler is stochastic and lives inline in
 * the page; it is not imported here. Where a tested function touches randomness
 * indirectly we assert PROPERTIES, but all functions below are deterministic.
 */
const {
  runDLMetaAnalysis,
  runFrequentistNMA,
  buildNetworkStructure,
  runSubgroupAnalysis,
  calculateHeterogeneityRisk,
  checkInconsistency,
  assessPublicationBias,
  runDoseResponseAnalysis,
  generateTimeToEventCurves,
  assessQualityFlags,
  calculateNCB,
  generateLeagueTable
} = require('./engine.js');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
  if (cond) { pass++; console.log('  ok  ' + name); }
  else { fail++; console.log(' FAIL ' + name + (detail ? '  -> ' + detail : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol || 1e-3); }

// ---------------------------------------------------------------------------
// runDLMetaAnalysis — DerSimonian-Laird pooling on the log scale.
// ---------------------------------------------------------------------------

// Hand-worked heterogeneous case (full derivation in README):
//   t1: logOR=ln(0.7)=-0.3566749, SE=0.1 -> v=0.01, w=100
//   t2: logOR=ln(0.9)=-0.1053605, SE=0.2 -> v=0.04, w=25
//   swFE=125; muFE=-38.301525/125 = -0.30641220
//   Q = 100*(-0.0502628)^2 + 25*(0.2010512)^2 = 0.2526349 + 1.0105395 = 1.2631744
//   df=1; C = 125 - (10000+625)/125 = 125 - 85 = 40
//   tau2 = (1.2631744-1)/40 = 0.0065793590
//   I2   = 0.2631744/1.2631744 * 100 = 20.8347 %
//   wRE = [60.31021, 21.46887]; swRE=81.77908;
//   muRE = (60.31021*ln0.7 + 21.46887*ln0.9)/81.77908 = -23.77480/81.77908 = -0.2907040
//   seRE = sqrt(1/81.77908) = 0.110580
const het = runDLMetaAnalysis(
  [
    { logOR: Math.log(0.7), SE: 0.1 },
    { logOR: Math.log(0.9), SE: 0.2 }
  ],
  'logOR', 'SE'
);
ok('DL het: tau2 ~ 0.00657936', close(het.tau2, 0.00657936, 1e-6), 'got ' + het.tau2);
ok('DL het: I2 ~ 20.8347%', close(het.I2, 20.8347, 1e-3), 'got ' + het.I2);
ok('DL het: muRE ~ -0.2907040', close(het.mu, -0.2907040, 1e-5), 'got ' + het.mu);
ok('DL het: seRE ~ 0.110580', close(het.se, 0.110580, 1e-5), 'got ' + het.se);

// Two IDENTICAL trials: no heterogeneity, tau2=0, I2=0, pooled == the study.
//   Each logOR=ln(0.8)=-0.2231435, SE=0.1 -> v=0.01, w=100.
//   swFE=200; muFE=-0.2231435; Q=0; tau2=max(0,(0-1)/C)=0; I2=max(0,(0-1)/0*100)=0.
//   seRE = sqrt(1/200) = 0.0707107.
const ident = runDLMetaAnalysis(
  [
    { logOR: Math.log(0.8), SE: 0.1 },
    { logOR: Math.log(0.8), SE: 0.1 }
  ],
  'logOR', 'SE'
);
ok('DL identical: tau2 == 0', ident.tau2 === 0, 'got ' + ident.tau2);
ok('DL identical: I2 == 0', ident.I2 === 0, 'got ' + ident.I2);
ok('DL identical: mu == ln(0.8)', close(ident.mu, Math.log(0.8), 1e-12), 'got ' + ident.mu);
ok('DL identical: se == sqrt(1/200)', close(ident.se, Math.sqrt(1/200), 1e-12), 'got ' + ident.se);

// k=1 edge case: MUST NOT be NaN. Single study -> mu=its logOR, se=its SE,
//   tau2=0, I2=0 (the df=max(1,k-1) + max(0,...) guards make C=0 -> -Inf -> 0).
const single = runDLMetaAnalysis([{ logOR: Math.log(0.75), SE: 0.12 }], 'logOR', 'SE');
ok('DL k=1: mu not NaN', !isNaN(single.mu), 'got ' + single.mu);
ok('DL k=1: mu == ln(0.75)', close(single.mu, Math.log(0.75), 1e-12), 'got ' + single.mu);
ok('DL k=1: se == 0.12', close(single.se, 0.12, 1e-12), 'got ' + single.se);
ok('DL k=1: tau2 == 0 (not NaN)', single.tau2 === 0, 'got ' + single.tau2);
ok('DL k=1: I2 == 0 (not NaN)', single.I2 === 0, 'got ' + single.I2);

// ---------------------------------------------------------------------------
// runFrequentistNMA — inverse-variance per-drug effect, sign-oriented.
// ---------------------------------------------------------------------------
// Single trial A(treatment) vs R(comparator), logOR_stroke=-0.4, SE=0.1.
//   For drug A: y=+(-0.4)=-0.4, w=100, mean=-0.4, se=0.1 -> ci ~ [-0.596, -0.204].
//   For drug R: y=-(-0.4)=+0.4 (R is comparator), mean=+0.4.
const fnma = runFrequentistNMA(
  [{ treatment: 'A', comparator: 'R', logOR_stroke: -0.4, SE_stroke: 0.1 }],
  ['R', 'A']
);
const drugA = fnma.posteriors.find(p => p.drug === 'A');
const drugR = fnma.posteriors.find(p => p.drug === 'R');
ok('FNMA: A mean == -0.4', close(drugA.mean, -0.4, 1e-12), 'got ' + drugA.mean);
ok('FNMA: A ci_low == -0.4-1.96*0.1', close(drugA.ci_low, -0.4 - 1.96*0.1, 1e-12));
ok('FNMA: R mean == +0.4 (sign flip)', close(drugR.mean, 0.4, 1e-12), 'got ' + drugR.mean);
// Drug with no trials -> neutral default.
const fnma2 = runFrequentistNMA([], ['Z']);
ok('FNMA: empty trials -> mean 0 default', fnma2.posteriors[0].mean === 0);

// ---------------------------------------------------------------------------
// buildNetworkStructure — edge counts and node sizes.
// ---------------------------------------------------------------------------
// Two trials both A vs B, n_total 100 and 300.
//   One merged edge (count=2, n=400); nodeSize[A]=nodeSize[B]=(100+300)/2=200.
const net = buildNetworkStructure(
  [
    { treatment: 'A', comparator: 'B', n_total: 100 },
    { treatment: 'A', comparator: 'B', n_total: 300 }
  ],
  ['A', 'B']
);
ok('NET: single merged edge', net.edges.length === 1);
ok('NET: edge count == 2', net.edges[0].count === 2);
ok('NET: edge n == 400', net.edges[0].n === 400);
ok('NET: nodeSize[A] == 200', net.nodeSize['A'] === 200, 'got ' + net.nodeSize['A']);

// ---------------------------------------------------------------------------
// runSubgroupAnalysis — deterministic class-effect interaction model.
// ---------------------------------------------------------------------------
// FXa drug (Apixaban), normal age/renal -> modifier 1.05, ageEffect=renalEffect=1.
//   adjEffect = -0.2 * 1.05 * 1 * 1 = -0.21.
const sg = runSubgroupAnalysis(null, ['Apixaban'], '<65', '>60');
ok('SUBGRP: FXa normal -> effect -0.21', close(sg[0].effect, -0.21, 1e-12), 'got ' + sg[0].effect);
// DTI drug (Dabigatran 150mg), severe renal <30, age >=75:
//   modifier = 1.4 (DTI & <30); ageEffect=1.15, renalEffect=1.25.
//   adjEffect = -0.2 * 1.4 * (1/1.15) * (1/1.25) = -0.28/1.4375 = -0.1947826.
const sg2 = runSubgroupAnalysis(null, ['Dabigatran 150mg'], '≥75', '<30');
ok('SUBGRP: DTI severe -> effect ~ -0.194783',
   close(sg2[0].effect, -0.2*1.4/1.15/1.25, 1e-9), 'got ' + sg2[0].effect);
ok('SUBGRP: ci width == 0.4', close(sg[0].ci_high - sg[0].ci_low, 0.4, 1e-12));

// ---------------------------------------------------------------------------
// calculateHeterogeneityRisk — rule-based I^2 heuristic.
// ---------------------------------------------------------------------------
// 3 studies, avgN=2000 (<5000), avgFU=1 (<1.5), yearSpread=12 (>10):
//   nStudies<5 +20; avgN<5000 +15; avgFU<1.5 +10; yearSpread>10 +25 = 70.
const hr = calculateHeterogeneityRisk([
  { n_total: 1000, followup_years: 1, year: 2005 },
  { n_total: 2000, followup_years: 1, year: 2011 },
  { n_total: 3000, followup_years: 1, year: 2017 }
]);
ok('HET-RISK: predictedI2 == 70', hr.predictedI2 === 70, 'got ' + hr.predictedI2);
// All four rules fire (k<5 +20, avgN<5000 +15, avgFU<1.5 +10, spread>10 +25 = 70),
// which is the maximum the rule set can produce; the Math.min(85,..) cap never
// binds. Document the true ceiling rather than an unreachable 85.
const hr2 = calculateHeterogeneityRisk([
  { n_total: 100, followup_years: 0.5, year: 1990 },
  { n_total: 100, followup_years: 0.5, year: 2020 }
]);
ok('HET-RISK: all-rules max == 70 (cap of 85 never binds)',
   hr2.predictedI2 === 70, 'got ' + hr2.predictedI2);
ok('HET-RISK: importances sum to 1.00',
   close(hr.featureImportance.reduce((a, f) => a + f.importance, 0), 1.0, 1e-9));

// ---------------------------------------------------------------------------
// assessPublicationBias — Egger-style precision regression + <3 guard.
// ---------------------------------------------------------------------------
ok('PUBBIAS: <3 trials -> null eggerP',
   assessPublicationBias([{ logOR_stroke: 0, SE_stroke: 0.1 }], [], {}).eggerP === null);
// 3 trials. Perfectly flat (all effects equal) -> slope 0, intercept == meanE.
//   effects all -0.2 -> intercept = -0.2, |intercept|>0.15 -> eggerP 0.03.
const pb = assessPublicationBias([
  { logOR_stroke: -0.2, SE_stroke: 0.1 },
  { logOR_stroke: -0.2, SE_stroke: 0.2 },
  { logOR_stroke: -0.2, SE_stroke: 0.05 }
], [], {});
ok('PUBBIAS: flat intercept == -0.2', close(pb.intercept, -0.2, 1e-9), 'got ' + pb.intercept);
ok('PUBBIAS: |intercept|>0.15 -> eggerP 0.03', pb.eggerP === 0.03);
ok('PUBBIAS: funnelData length == 3', pb.funnelData.length === 3);

// ---------------------------------------------------------------------------
// runDoseResponseAnalysis — OLS of effect on dose.
// ---------------------------------------------------------------------------
// Two drugs at known doses with a clean linear effect.
//   Apixaban dose 5 -> filtered out (dose>5 required). Use Edoxaban 30 & 60.
//   effects: 30->-0.1, 60->-0.4. slope=(-0.4 - -0.1)/(60-30) = -0.3/30 = -0.01.
//   intercept = meanEffect - slope*meanDose = -0.25 - (-0.01)*45 = -0.25 + 0.45 = 0.2.
const dr = runDoseResponseAnalysis([
  { treatment: 'Edoxaban 30mg', logOR_stroke: -0.1, SE_stroke: 0.1 },
  { treatment: 'Edoxaban 60mg', logOR_stroke: -0.4, SE_stroke: 0.1 }
], []);
ok('DOSE: slope == -0.01', close(dr.slope, -0.01, 1e-9), 'got ' + dr.slope);
ok('DOSE: intercept == 0.2', close(dr.intercept, 0.2, 1e-9), 'got ' + dr.intercept);

// ---------------------------------------------------------------------------
// generateTimeToEventCurves — deterministic survival curve.
// ---------------------------------------------------------------------------
// Drug idx0: baseHazard=0.03, hr=exp(-0.15). At t=0 survival=100.
//   At t=12: survProb=exp(-0.03*exp(-0.15)*12/12)=exp(-0.03*0.860708)=exp(-0.0258212)=0.974509 ->97.4509.
const ttc = generateTimeToEventCurves([], ['A', 'B']);
ok('TTE: t=0 survival == 100', close(ttc[0].points[0].survival, 100, 1e-9));
ok('TTE: monotone non-increasing',
   ttc[0].points.every((p, i, a) => i === 0 || p.survival <= a[i-1].survival + 1e-9));
ok('TTE: t=12 survival ~ 97.4509',
   close(ttc[0].points[4].survival, Math.exp(-0.03*Math.exp(-0.15))*100, 1e-6),
   'got ' + ttc[0].points[4].survival);

// ---------------------------------------------------------------------------
// assessQualityFlags — GRADE-style downgrades.
// ---------------------------------------------------------------------------
// Clean evidence: low RoB, tight CIs, consistent, no pub bias -> HIGH, all 0.
const mcmcClean = { posteriors: [{ ci_low: -0.3, ci_high: 0.1 }, { ci_low: -0.2, ci_high: 0.2 }] };
const gradeHi = assessQualityFlags(mcmcClean, [], { eggerP: 0.5 },
  [{ rob: 'low' }, { rob: 'low' }]);
ok('GRADE: clean -> HIGH', gradeHi.overall === 'HIGH', 'got ' + gradeHi.overall);
ok('GRADE: clean -> rob 0', gradeHi.rob === 0);
// All four downgrades fire -> VERY LOW. >50% high RoB; avg inconsistency>0.3;
//   wide CIs (>0.8); eggerP<0.05.
const wide = { posteriors: [{ ci_low: -1, ci_high: 1 }] };
const gradeLo = assessQualityFlags(wide,
  [{ difference: 0.5 }], { eggerP: 0.01 },
  [{ rob: 'high' }, { rob: 'high' }]);
ok('GRADE: all downgrades -> VERY LOW', gradeLo.overall === 'VERY LOW', 'got ' + gradeLo.overall);
ok('GRADE: all downgrades -> rob -1, incon -1, pubbias -1',
   gradeLo.rob === -1 && gradeLo.incon === -1 && gradeLo.pubbias === -1);

// ---------------------------------------------------------------------------
// calculateNCB — stroke-benefit minus weighted bleed harm.
// ---------------------------------------------------------------------------
// One drug, posterior mean=ln(0.5)=-0.693147 -> rrStroke=0.5.
//   strokeRisk=0.04, bleedRisk=0.02.
//   strokeReduction = 0.04*(1-0.5) = 0.02.
//   rrBleed = 1 + (0.5-1)*0.5 = 0.75; bleedIncrease = 0.02*(0.75-1) = -0.005.
//   ncb = 0.02 - (-0.005*2.5) = 0.02 + 0.0125 = 0.0325.
const ncb = calculateNCB({ posteriors: [{ drug: 'X', mean: Math.log(0.5) }] }, 0.04, 0.02, 0);
ok('NCB: == 0.0325', close(ncb[0].ncb, 0.0325, 1e-9), 'got ' + ncb[0].ncb);
ok('NCB: single drug probBest == 100', ncb[0].probBest === 100);

// ---------------------------------------------------------------------------
// generateLeagueTable — pairwise OR contrasts.
// ---------------------------------------------------------------------------
// Two drugs, means -0.2 and 0.1. Cell[0][1] = drugs[1] vs drugs[0]:
//   diff = mean[1]-mean[0] = 0.1-(-0.2)=0.3 -> or=exp(0.3)=1.349859.
const lt = generateLeagueTable({ posteriors: [{ mean: -0.2 }, { mean: 0.1 }] }, ['A', 'B']);
ok('LEAGUE: diagonal OR == 1', lt[0][0].or === 1.0);
ok('LEAGUE: B vs A OR == exp(0.3)', close(lt[0][1].or, Math.exp(0.3), 1e-9), 'got ' + lt[0][1].or);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
