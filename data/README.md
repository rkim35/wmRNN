## Getting the data
The experimental data used in the study can be downloaded from http://crcns.org/data-sets/pfc/pfc-3/

## Neuronal timescales
To compute and extract neuronal timescales from the experimental dataset, use `autocorr_analysis_pfc3.m`. The script uses the excel file ("SummaryDatabase.xlsx" included in the dataset) to first identify a group of neurons. Next, it loads the spike train data and computes the neuronal timescale for each neuron.

The script is written to group neurons by three variables: learning condition, task type, and prefrontal cortex area. Refer to the script for more info on these variables.


