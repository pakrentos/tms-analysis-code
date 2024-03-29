This git repository contains scripts written in FieldTrip (MATLAB) and Python that are used for various purposes related to sensor-level and statistical analysis of EEG data, such as extracting motor imagery brain response time (MIBRT, a measure of brain activity), estimation and analysis of source-level power distributions, reconstruction of functional connectivities in source space, and analyzing the properties of reconstructed functional connectivities.

The scripts in the FieldTrip folder:

#### EEG_eLoreta.m and EEG_eLoreta_BgrsCalc.m
Used for the estimation of source-level power distributions using the eLORETA method.

#### Stat_source_power_distributions.m
Used for the statistical analysis of source-level power distributions in the between-subjects design.

#### Source power in cluster.m
Used for the calculation of source power in the selected cluster revealed in the statistical test.

#### Source power in zone.m
Used for the calculation of source power in the selected zone.

#### PLV_connectivity.m
Used for the calculation of PLV connectivity measure between the selected ROIs in the source space.

The scripts in the Python folder:

#### additional_group_preproc.ipynb
Used for data preprocessing

#### plv_coherence_stat.ipynb
Statistical analysis of the PLV network characteristic

#### sham_preproc.ipynb
Sham condition data preprocessing

#### stfr_analysis.ipynb
Spatio-time-frequency sensor-level analysis for preprocessed EEG data. Written in MATLAB language

#### TMS_mbrt_by_stat_final_feb-2022.ipynb
MBRT calculation using statistical methods (one-sample t-test for each subject)

#### TMS_mbrt_by_stat_plots_final_feb-2022.ipynb
Plotting the results of the TMS_mbrt_by_stat_final_feb-2022.ipynb for manual verification

#### TMS_src_stat.ipynb and src_power_stat.ipynb
Statistical analysis of reconstructed source power data 
