# Confounder_Selection_w_Causal_Discovery
 This repository stores the codes for the work:
 
 ```  
@article{zhong_challenges_2021,
    title = {Challenges in automating confounder selection with causal discovery methods.},
    author = {{Yongqi Zhong} and {Sofia Triantafillou} and {Edward H. Kennedy} and {Maria M. Brooks} and {Lisa M. Bodnar} and {Ashley I. Naimi}},
    journal = {Epidemiology},
    year = {2021},
    note = {Submitted}
}
```

- Simulations: Plasmode simulations using data from the Effects of Aspirin in Gestation and Reproduction (EAGeR) trial
- Estimation: Estimating the average treatment effects (ATE) on the simulated data using different confounder selection scenarios
    - With knowledge of the true data generating mechanisms: est_wo_cd.R
    - Using the Min-Max Hill-Climbing (MMHC) causal discovery method (default setting): est_cd.R
    - Using the Min-Max Hill-Climbing (MMHC) causal discovery method (tuned setting): est_cd_mmhc_tuned.R
- Analysis: Summarizing results from the estimated ATEs
    - Calculate absolute bias and MSE of ATEs estimated by different confounder selection scenarios: summarise_wo_CD.R, summarise_CD.R, summarise_CD_mmhc_tuned.R
    - Calculate accuracy of MMHC algorithm: summarise_all.R