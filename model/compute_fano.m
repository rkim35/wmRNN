% Name: Robert Kim
% Date: 01-25-2020
% Email: rkim@salk.edu
% compute_fano.m
% Description: Script to compute the Fano factor for the model data

clear; clc;
addpath('/cnl/chaos/ROBERT/wm_intrinsic_timescales/code/matlab/sim');


% DMS RNNs
task_type = 'xor';
num_mods = 40;
fr_lim = [0 100];
perf_thr = 0.95;
model_dir = fullfile(pwd, 'DMS_OSF')
dms_out = fnc_get_models(model_dir, task_type, fr_lim, true, perf_thr, num_mods);

% AFC RNNs
task_type = 'afc';
num_mods = 40;
fr_lim = [0 100];
perf_thr = 0.95;
model_dir = fullfile(pwd, 'AFC_OSF')
afc_out = fnc_get_models(model_dir, task_type, fr_lim, true, perf_thr, num_mods);

% Compute the fano factor
spks = afc_out.spks;
afc_fanos = nan(1, length(spks));
afc_means = nan(1, length(spks));
afc_vars = nan(1, length(spks));
afc_consecs = nan(1, length(spks));
for i = 1:length(spks)
  curr_spk = spks{i}; % [50 x 20]
  consec = zeros(size(curr_spk, 1), 1);
  for jj = 1:size(curr_spk, 1)
    non_firing_bins = find(curr_spk(jj, :) == 0);
    if isempty(non_firing_bins)
      consec(jj) = size(curr_spk, 2);
    else
      nf = [];
      for jjj = 1:length(non_firing_bins)
        if jjj == 1
          nf = [nf non_firing_bins(jjj)-1];
        elseif jjj > 1 & jjj < length(non_firing_bins)
          nf = [nf non_firing_bins(jjj+1) - non_firing_bins(jjj) - 1];
        else
          nf = [nf size(curr_spk, 2) - non_firing_bins(jjj)];
        end
      end
      consec(jj) = sum(nf);
    end
  end

  sum_spk = sum(curr_spk, 2);
  afc_means(i) = nanmean(sum_spk);
  afc_vars(i) = nanvar(sum_spk);
  afc_fanos(i) = nanvar(sum_spk)/nanmean(sum_spk);
  afc_consecs(i) = nanvar(consec);
end

spks = dms_out.spks;
dms_fanos = nan(1, length(spks));
dms_means = nan(1, length(spks));
dms_vars = nan(1, length(spks));
dms_consecs = nan(1, length(spks));
for i = 1:length(spks)
  curr_spk = spks{i}; % [50 x 20]
  consec = zeros(size(curr_spk, 1), 1);
  for jj = 1:size(curr_spk, 1)
    non_firing_bins = find(curr_spk(jj, :) == 0);
    if isempty(non_firing_bins)
      consec(jj) = size(curr_spk, 2);
    else
      nf = [];
      for jjj = 1:length(non_firing_bins)
        if jjj == 1
          nf = [nf non_firing_bins(jjj)-1];
        elseif jjj > 1 & jjj < length(non_firing_bins)
          nf = [nf non_firing_bins(jjj+1) - non_firing_bins(jjj) - 1];
        else
          nf = [nf size(curr_spk, 2) - non_firing_bins(jjj)];
        end
      end
      consec(jj) = sum(nf);
    end
  end

  sum_spk = sum(curr_spk, 2);
  dms_means(i) = nanmean(sum_spk);
  dms_vars(i) = nanvar(sum_spk);
  dms_fanos(i) = nanvar(sum_spk)/nanmean(sum_spk);
  dms_consecs(i) = nanvar(consec);
end









