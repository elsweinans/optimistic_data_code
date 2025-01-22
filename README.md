# Data and code used for the paper "Signals of complexity and fragmentation in accelerometer data" (submitted to PLOS One).

For this paper, we re-analyzed data that was obtained for the OPTIMISTIC trial (ClinicalTrials.gov number: NCT02118779). The raw data can be find in the folder 'raw_data'. The accelerometer data was obtained via a triaxial accelerometer (GENEActiv, ActivInsights Ltd, UK) which was ankle worn and gave a signal every 5 seconds.

This data was used to fill in the table 'dataBL.mat'. This table has rows for each indvidual, with a column for age, sex, correlation dimension (corrdims), kAR, and kRA, which can be calculated from the raw data. This matlab table can be used in 'analysis_fulldataset.m' to reproduce the figures in the manuscript.

For more information about the steps and the interpretation of results, we refer to our manuscript.
