## Getting data
First, download the trained RNNs from [here](https://osf.io/md4wg/). For the DMS RNN model, download and unzip "DMS\_OSF.zip". The folder contains `.mat` files where each file corresponds to one trained RNN for the DMS task. Each `.mat` file contains several important parameters including neuronal timescales (extracted using the method similar to the one developed by [Murray et al.](https://www.nature.com/articles/nn.3862)).

NOTE: The DMS data folder (DMS\_OSF) contains both poor and good performance DMS RNNs.

## Neuronal timescales
The main function for extracting parameters from a group of trained RNNs is called `fnc_get_models.m`. Run the following MATLAB code to extract parameters from the good performance DMS RNNs:

```
% Load the good performance DMS RNNs
task_type = 'xor'; % Task name (either 'xor' or 'afc')
num_mods = 40; % number of models to load
fr_lim = [0 100]; % firing rate threshold ([min max]) in Hz
perf_thr = 0.95; % include only RNNs whose task performance > perf_thr
model_dir = <DIRECTORY TO THE DMS_OSF FOLDER>;
dms_out = fnc_get_models(model_dir, task_type, fr_lim, true, perf_thr, num_mods);
```

The above code will return all the DMS RNNs whose average task performance is greater than 95% (i.e., 0.95). `dms_out` is a struct variable containing all the parameters pulled from the trained RNNs. To plot the distribution of the neuronal timescales from this group of RNNs (i.e., Fig. 2c center), run the following code:

```
figure('Units', 'Inch', 'Outerposition', [0 0 5 4]);
[~, edges] = histcounts(log10(dms_out.taus), [0.7:0.1:2.7]);
histogram(dms_out.taus, 10.^edges, 'Normalization', 'probability', 'FaceColor', 'b');
set(gca, 'XScale', 'log'); axis tight; hold on;
```

To plot the autocorrelation curves (i.e., Fig. 2d center), run the following:

```
figure('Units', 'Inch', 'Outerposition', [0 0 5 4]);
axis tight; hold on;
plot(dms_out.auto_c(1:1:end, :)', 'b');
plot(nanmean(dms_out.auto_c(:, :)), 'wo-', 'linewidth', 2, 'markers', 8, 'MarkerFace', 'b')
xlim([0.5, 12])
```

## Cross-temporal discriminability analysis
The main script for the cross-temporal analysis (i.e., Fig. 3) is `cross_temporal_disc.m`. See the Methods in the manuscript for more details. Running the script will generate a MATLAB file (`cross_temporal_results_quantile.mat`). The following code can be used to plot the cross-temporal discriminability matrices:

```
load("cross_temporal_results_quantile.mat")

stim_on = 0.15;
stim_dur = 0.25;
delay = 0.75;
sim_t = t(1:100:end) - stim_on;

% SHORT TAU GROUP
figure('Units', 'Inch', 'Outerposition', [0 0 4 4]);
imagesc(sim_t, sim_t, corrs_short, [0, 1.00]);
set(gca, 'YDir', 'normal');
xlim([-stim_on, stim_dur+delay]); ylim([-stim_on, stim_dur+delay]); hold on;
plot([0 0], ylim, 'w-', 'linewidth', 3);
plot(xlim, [0 0], 'w-', 'linewidth', 3);
plot([stim_dur stim_dur], ylim, 'w-', 'linewidth', 3);
plot(xlim, [stim_dur stim_dur], 'w-', 'linewidth', 3);
axis square; hold on;

% LONG TAU GROUP
figure('Units', 'Inch', 'Outerposition', [0 0 4 4]);
imagesc(sim_t, sim_t, corrs_long, [0, 1.00]);
set(gca, 'YDir', 'normal');
xlim([-stim_on, stim_dur+delay]); ylim([-stim_on, stim_dur+delay]); hold on;
plot([0 0], ylim, 'w-', 'linewidth', 3);
plot(xlim, [0 0], 'w-', 'linewidth', 3);
plot([stim_dur stim_dur], ylim, 'w-', 'linewidth', 3);
plot(xlim, [stim_dur stim_dur], 'w-', 'linewidth', 3);
axis square; hold on;
```


