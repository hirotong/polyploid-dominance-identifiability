# Panel-conditioned diagnostic workflow - summary

Generated from `workflow_scorecard.csv` over 12 panels. All reference controls are conditioned on each panel's own genotype matrix (no universal n/M ruler). H2=0.6, n_perm/n_ref per run settings.

## Decision rule
Biological digenic (D2) dominance is **not** claimed unless a panel passes BOTH:
(i) panel-conditioned D2 recoverability (a planted s=0.5 D2 spike is the largest non-additive component), AND
(ii) real-trait D2 localization above its own permutation-null / additive-only baseline,
with structure control not removing the signal and a non-wide D2 credible interval. Prediction gain is reported separately and never used to claim biological dominance.

## 1. Substantial fitted D-labelled share (s_hat_AD >= 0.30)
Panels: Alfalfa, Potato-SolCAP, Potato-SLU, BB-UF-GP, BB-UF-fruit, BB-UF-call*, Panicum, Rose, Miscanthus, Sweetpotato, Sugarcane, Positive control (sim)

## 2. Exceed their own permutation null (delta_null > margin AND real share above perm q95)
Panels: Alfalfa, Potato-SLU, Positive control (sim)

## 3. Structure-sensitive (S_struct_real >= 0.10 or above panel structure-ruler q95)
Panels: _none_

## 4. Can recover a panel-conditioned D2 spike-in (recovery rate >= 0.80 at s=0.5)
Recoverable: Alfalfa, Potato-SLU, BB-UF-GP, BB-UF-fruit, BB-UF-call*, Panicum, Sweetpotato, Positive control (sim)
NOT recoverable (geometry-limited): Potato-SolCAP, Rose, Miscanthus, Sugarcane

## 5. Real traits resembling the D2 spike-in controls (P_D2 closer to spike than to additive-null AND L_D2 above additive-null q95)
Panels: Positive control (sim)

## 6. Material prediction gain from standard A+D (mean delta_AD >= 0.01 and CI lower bound > 0)
Panels: Alfalfa, Potato-SolCAP, Potato-SLU, Sugarcane, Positive control (sim)

## 7. Recommended reporting language per panel
- **Alfalfa** -> _prediction-useful_  (s_AD=0.61, delta_null=+0.06, z=1.9, D2 recov=0.83, P_D2=0.26, S_struct=0.01, dA+D=+0.052)
- **Potato-SolCAP** -> _prediction-useful_; unresolved by panel geometry (D2 spike not recoverable)  (s_AD=0.49, delta_null=-0.05, z=-1.5, D2 recov=0.73, P_D2=0.21, S_struct=0.00, dA+D=+0.020)
- **Potato-SLU** -> _prediction-useful_  (s_AD=0.53, delta_null=+0.07, z=3.1, D2 recov=1.00, P_D2=0.23, S_struct=0.05, dA+D=+0.012)
- **BB-UF-GP** -> _model-basis_  (s_AD=0.46, delta_null=-0.06, z=-3.2, D2 recov=1.00, P_D2=0.17, S_struct=0.03, dA+D=+0.003)
- **BB-UF-fruit** -> _model-basis_  (s_AD=0.42, delta_null=-0.22, z=-8.7, D2 recov=1.00, P_D2=0.17, S_struct=0.04, dA+D=-0.001)
- **BB-UF-call*** -> _model-basis_  (s_AD=0.41, delta_null=-0.08, z=-4.7, D2 recov=1.00, P_D2=0.15, S_struct=0.07, dA+D=-0.001)
- **Panicum** -> _model-basis_  (s_AD=0.50, delta_null=-0.01, z=-0.5, D2 recov=0.90, P_D2=0.22, S_struct=0.02, dA+D=+0.002)
- **Rose** -> _model-basis_; unresolved by panel geometry (D2 spike not recoverable)  (s_AD=0.50, delta_null=-0.02, z=-0.6, D2 recov=0.67, P_D2=0.24, S_struct=0.01, dA+D=-0.007)
- **Miscanthus** -> _model-basis_; unresolved by panel geometry (D2 spike not recoverable)  (s_AD=0.47, delta_null=-0.02, z=-0.7, D2 recov=0.73, P_D2=0.19, S_struct=0.02, dA+D=-0.000)
- **Sweetpotato** -> _model-basis_  (s_AD=0.48, delta_null=-0.03, z=-1.2, D2 recov=0.87, P_D2=0.11, S_struct=0.00, dA+D=+0.002)
- **Sugarcane** -> _prediction-useful_; unresolved by panel geometry (D2 spike not recoverable)  (s_AD=0.50, delta_null=-0.02, z=-0.7, D2 recov=0.77, P_D2=0.18, S_struct=0.00, dA+D=+0.011)
- **Positive control (sim)** -> _D2-supported_  (s_AD=0.76, delta_null=+0.27, z=6.8, D2 recov=1.00, P_D2=0.74, S_struct=0.00, dA+D=+0.192)

## Notes
- A panel that fails BOTH real-trait D2 localization and D2 spike recovery is **unresolved by panel geometry**, not evidence of absent biological dominance (sanity check 4).
- A panel with material prediction gain but failed localization is **prediction-useful**, not biologically D2-supported (sanity check 5).
- Where real s_hat_AD sits at the permutation null, the fitted D-labelled share is a model-basis artefact, not dominance (sanity check 6).

