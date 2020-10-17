% Name: Robert Kim
% Date: 07-25-2019
% Email: rkim@salk.edu
% cross_temporal_disc.m
% Description: Discriminality analysis for the simulated data for long and short
% tau subgroups.

clear; clc;

% Trained model directory
task_path = fullfile(pwd, 'DMS_OSF');

perf_threshold = 0.95;
task_type = 'xor';
wcard = '*Taus_*';
stable_models = return_stable(task_path, wcard, perf_threshold, task_type);
disp(['Performance threshold set to ' num2str(perf_threshold)]);

% Used to get the intrinsic timescale quantiles
fr_lim = [0 100];
num_mods = 40;
ns_out = fnc_get_models(task_path, task_type, fr_lim, true, perf_threshold, num_mods);
low_qt = quantile(ns_out.taus, 0.25)
high_qt = quantile(ns_out.taus, 0.75)

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
shuffle_weights = false;
if shuffle_weights == true
  disp('SHUFFLING WEIGHTS')
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% Number of first stim classes (i.e 2)
num_classes = 2;

% Number of trials per neuron
num_trials = 12;

% Network size
num_neus = 200;

% Number of inputs (1 for 2afc and 2 for xor)
num_inputs = 2;

data = struct();
for i = 1:num_classes
  data(i).rs = [];
  data(i).neu = [];
end

% Two groups of data
dataA = data;
dataA_long = data;
dataA_short = data;

dataB = data;
dataB_long = data;
dataB_short = data;

% Input stim
input_stim = [-1, 1];

for md = 1:length(ns_out.decays)
  md
  model_name = stable_models{md};

  model_path = fullfile(task_path, model_name);

  model_path

  load(model_path);

  % Scaling factor and sampling rate
  scaling_factor = opt_scaling_factor;
  down_sample = 1;

  if ~exist('use_this_W')
    use_this_W = [];
  else
    disp('using the frozen weights!')
  end

  shuffle_w_in = [];

  if ~isempty(shuffle_w_in)
    disp('W_IN SHUFFLED!!!!!!!!!!!!!');
  end

  % Run the trained model once to get some params
  stims = struct();
  stims.mode = 'none';
  u = zeros(num_inputs, 411);
  [W, REC, spk, rs, all_fr, out, params] = LIF_network_fnc(model_path, scaling_factor,...
  u, stims, down_sample, shuffle_weights, use_this_W, shuffle_w_in);
  dt = params.dt;
  t = dt:dt:params.T;
  syn_decay = params.td*1000;
  T = length(out)/10; % total number of time-points

  % Now generate multi-trial data for each class (condition)
  for cn = 1:num_classes
    tmpA = [];
    tmpB = [];

    % Get long taus
    long_N = find(taus_decay_ms > high_qt); % quantile

    % Get short taus
    short_N = find(taus_decay_ms < low_qt); % quantile

    % For each trial
    parfor nt = 1:num_trials
      u = zeros(num_inputs, 411);
      if num_inputs == 2
        u(1, 31:80) = input_stim(cn);
        u(2, 231:280) = input_stim(cn);
      elseif num_inputs == 1
        u(1, 31:80) = input_stim(cn);
        u(1, 231:280) = input_stim(cn);
      end

      stims = struct();
      stims.mode = 'none';
      [W, REC, spk, rs, all_fr, out, params] = LIF_network_fnc(model_path, scaling_factor,...
      u, stims, down_sample, shuffle_weights, use_this_W, shuffle_w_in);

      % Gaussian filter width for smoothing the spike data
      gau_width = 50; % in ms
      gau_win = gau_width*(1/dt)/1000; 
      gau = gausswin(gau_win);
      spk_gaus = spk(auto_N, :);
      for ne = 1:size(spk_gaus, 1)
        spk_gaus(ne, :) = filter(gau, 1, spk_gaus(ne, :))*(1000/gau_width);
      end

      if mod(nt, 2) == 1 % odd trials
        tmpA = [tmpA; spk_gaus(:, 1:100:end)];
      else
        tmpB = [tmpB; spk_gaus(:, 1:100:end)];
      end
    end

    mean_dataA = zeros(length(auto_N), size(tmpA, 2));
    mean_dataB = zeros(length(auto_N), size(tmpB, 2));
    for mm = 1:length(auto_N)
      mean_dataA(mm, :) = nanmean(tmpA(mm:length(auto_N):end, :));
      mean_dataB(mm, :) = nanmean(tmpB(mm:length(auto_N):end, :));
    end

    dataA(cn).rs = [dataA(cn).rs; mean_dataA];
    dataB(cn).rs = [dataB(cn).rs; mean_dataB];

    dataA_long(cn).rs = [dataA_long(cn).rs; mean_dataA(long_N, :)];
    dataB_long(cn).rs = [dataB_long(cn).rs; mean_dataB(long_N, :)];

    dataA_short(cn).rs = [dataA_short(cn).rs; mean_dataA(short_N, :)];
    dataB_short(cn).rs = [dataB_short(cn).rs; mean_dataB(short_N, :)];
  end
  clearvars -except dataA* dataB* input_stim num_trials num_classes num_inputs ...
  data task_* wcard stable_models shuffle_weights num_neus t ns_out ...
  low_qt high_qt num_mods bad_model_names shuffle_w_in;
end

% -------------------------------------------------------
% Cross-temporal discriminability analysis
% -------------------------------------------------------
disp('Performing the cross-temporal analysis...')
t_pts = size(dataA(1).rs, 2);
corrs = zeros(t_pts, t_pts);
parfor t1 = 1:t_pts
  t1
  for t2 = 1:t_pts
    currA_diff = dataA(1).rs(:, t1) - dataA(2).rs(:, t1);
    currB_diff = dataB(1).rs(:, t2) - dataB(2).rs(:, t2);
    curr_corr = corr(currA_diff, currB_diff);
    corrs(t1, t2) = curr_corr;
  end
end

t_pts = size(dataA(1).rs, 2);
corrs_long = zeros(t_pts, t_pts);
parfor t1 = 1:t_pts
  t1
  for t2 = 1:t_pts
    currA_diff = dataA_long(1).rs(:, t1) - dataA_long(2).rs(:, t1);
    currB_diff = dataB_long(1).rs(:, t2) - dataB_long(2).rs(:, t2);
    curr_corr = corr(currA_diff, currB_diff);
    corrs_long(t1, t2) = curr_corr;
  end
end

t_pts = size(dataA(1).rs, 2);
corrs_short = zeros(t_pts, t_pts);
parfor t1 = 1:t_pts
  t1
  for t2 = 1:t_pts
    currA_diff = dataA_short(1).rs(:, t1) - dataA_short(2).rs(:, t1);
    currB_diff = dataB_short(1).rs(:, t2) - dataB_short(2).rs(:, t2);
    curr_corr = corr(currA_diff, currB_diff);
    corrs_short(t1, t2) = curr_corr;
  end
end

if ~isempty(shuffle_w_in) 
  out_fname = ['cross_temporal_results_quantile_' 'w_in_retrained.mat'];
else
  out_fname = 'cross_temporal_results_quantile.mat';
end

save(fullfile(task_path, out_fname), '-v7.3');



