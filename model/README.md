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
[~, edges] = histcounts(log10(dms_out.taus), [0.7:0.1:2.7]);                                             histogram(dms_out.taus, 10.^edges, 'Normalization', 'probability', 'FaceColor', 'b');                    set(gca, 'XScale', 'log'); axis tight; hold on;     
```

To plot the autocorrelation curves (i.e., Fig. 2d center), run the following:

```
figure('Units', 'Inch', 'Outerposition', [0 0 5 4]);
axis tight; hold on;
plot(dms_out.auto_c(1:1:end, :)', 'b');
plot(nanmean(dms_out.auto_c(:, :)), 'wo-', 'linewidth', 2, 'markers', 8, 'MarkerFace', 'b')
xlim([0.5, 12])
```



