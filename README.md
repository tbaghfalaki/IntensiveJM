# Dynamic Prediction in Intensive Longitudinal Studies
**Background:**
Intensive longitudinal data, characterized by frequent repeated measurements within individuals, provide detailed information on dynamic processes but present substantial statistical and computational challenges for dynamic risk prediction. We systematically compared alternative modelling strategies to assess their predictive performance and computational efficiency under different longitudinal data structures.

**Methods:**
We compared joint modelling (JM), landmarking, multiple two-stage modelling, Multivariate Functional Principal Component Cox (MFPCCox), Penalized Regression Calibration (Pencal), and TransformerJM. Within the JM framework, we additionally evaluated three measurement-reduction strategies: random selection, systematic sampling, and aggregation. Predictive performance was assessed using time-dependent area under the curve (AUC) and Brier score. Two simulation scenarios considered multivariate longitudinal markers with either linear or nonlinear trajectories and subject-specific variability. Methods were further evaluated using data from 201 patients with subarachnoid hemorrhage, including hourly systolic blood pressure and heart rate measurements, to predict cerebral vasospasm.

**Results:**
In the linear simulation scenario, JM approaches achieved predictive performance close to the oracle benchmark, while landmarking with linear mixed models was also competitive. TransformerJM provided similar discrimination with lower Brier scores and substantially lower computational cost. Under nonlinear trajectories with correlated markers and subject-specific variability, no method consistently dominated across landmark times and performance measures. The JM variants produced predictive performance very similar to the full JM, while measurement reduction substantially reduced computational time. In the clinical application, predictive performance varied markedly with landmark time. TransformerJM showed the strongest early performance, with AUCs of 0.891 and 0.773 and Brier scores of 0.014 and 0.018 at 48 and 72 hours, respectively. JM-Aggregation achieved an AUC of 0.864 at 72 hours, while MFPCCox also showed strong performance at early landmarks. Predictive performance generally deteriorated at later landmark times.

**Conclusions:**
No single method was uniformly optimal across simulation and clinical settings. JM performed particularly well when its longitudinal structure was correctly specified, whereas more flexible approaches were advantageous for complex nonlinear trajectories. TransformerJM offered an attractive balance of predictive accuracy and computational efficiency, particularly for early prediction, while measurement reduction provided an effective strategy for reducing JM computational burden without substantial loss of predictive performance. Selected methods should therefore consider longitudinal complexity, prediction time, and computational constraints.


# Describtion of this page
This page contains simulation study code corresponding to   scenario 1 presented in the paper.


### Reference 
Baghfalaki, T., Ganjali, M., & Jacqmin‑Gadda, H. others (2026). Methodological strategies for dynamic prediction in intensive longitudinal studies. *Revised*.
