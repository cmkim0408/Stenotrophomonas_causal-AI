# Runtime profiling summary

- environment: 3.12.11 on Windows-10-10.0.19045-SP0
- CPUs: 16

| Stage | n | wall (s) | peak RSS (MB) |
|---|---|---|---|
| lhs_sample | 500 | 0.0003 | 145.8 |
| fba_batch | 50 | 2.0686 | 281.2 |
| targeted_fva | 10 | 0.4246 | 289.4 |
| classifier+shap | 242 | 0.2854 | 400.9 |
| regressor+shap | 242 | 0.0598 | 402.3 |
| pc_bootstrap | 25 | 0.3154 | 417.0 |
