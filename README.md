# solar_flare


![workflow](images/workflow2.png)

![SHAP clusters](images/subcluster_breakdown_dir_exact.png)


- step1_handle_missing.R --imputation for missing data 
- step2_make_solar_data.R --create window-based features  
- step3_big_xgb_brute_force.R --fit variable selection model
- step4_xgb_slim_tune.R --tune final model with 112 predictors
- step4_xgb_slim_tune.sh 
- step5_fit_and_submit.R --submit model
- make_sumdt_window_func.R --utility functions including feature calculation by window
- shap_plot_funcs.R --plotting functions
