% Name: Robert Kim
% Date: 08-09-2019
% Email: rkim@salk.edu
% fnc_get_models.m
% Description: Function to get important parameters including the estimated
% neuronal timescales from all the trained RNNs

function out = fnc_get_models(model_dir, task_type, fr_lim, pop_fit, perf_thr, num_mods)
% INPUT
%   model_dir: directory containing all the trained models
%   task_type: either 'xor' or 'mante'
%   fr_lim: firing rate limits ([min max]; inclusive)
%   pop_fit: true for fitting the population autocorrelation
%   perf_thr: performance threshold (0.5 -> 50% acc.)
%   num_mods: number of models to load
%
% OUTPUT
%   out: struct variable containing 

if nargin < 4
  pop_fit = true;
  perf_thr = 0.95;
  num_mods = NaN;
elseif nargin < 5
  perf_thr = 0.95;
  num_mods = NaN;
elseif nargin < 6
  num_mods = NaN;
end

dirs = dir(fullfile(model_dir, '*.mat*'));

all_taus = [];        % neuronal timescales (concatenated across all the networks)
all_tds = [];         % trained synaptic decay time-constants (concatenated across all the networks)
all_tds0 = [];        % initial synaptic decay time-constants (concatenated across all the networks)
all_auto_c = [];      % auto-correlation curves (concatenated across all the networks)
mean_auto_c = [];     % mean auto-correlation curves (averaged within each network first)
all_decays = [];      % neuronal timescales (concatenated across all the networks)
all_decays_mean = []; % mean neuronal timescales (averaged within each network first)
all_inh = [];         % binary vector indicating which units with neuronal timescales are inhibitory
all_frs = [];         % average firing rates (concatenated across all the networks)
all_trs = [];         % number of training trials required to train each network
perfs = [];           % average task performance for each spiking network
all_scaling_factors = [];   % scaling factor used for each network (for conversion to spiking network)

% Some variables for debugging
all_spks = {};        
neu_ind = {};
neu_tau_ind = {};

% Total number of connections (pooled from all the networks) 
% separated by the four connection types
% [ii, ie, ei, ee] which represents [I->I, E->I, I->E, E->E]
num_connections = zeros(1, 4); 

% Synaptic weights separated by synaptic connection types
% (pooled from all the networks)
w_inhs = [];
w_excs = [];
w_iis = [];
w_ees = [];
w_ies = [];
w_eis = [];

% For storing .mat filenames
fnames = {};


counter = 0;

% Loop over each of the .mat files in the specified directory
for d = 1:length(dirs);
  if ~isempty(strfind(dirs(d).name, 'Tau')) 
    load(fullfile(model_dir, dirs(d).name));

    if strcmpi(task_type, 'xor')
      % For the DMS (aka XOR) model, use the stable performance
      % measure (stable_perfs) for the task performance measure
      % Briefly, the stable performance measure refers to the performance
      % on the long delay (750-ms) DMS task. Remember that each DMS network
      % was first trained to perform the easier, short delay (50-ms) DMS task.
      perf_measure = nanmean(stable_perfs);
    else
      % For the AFC task, there is no "stable_perfs" and just use
      % the best task performance computed from the lambda grid search
      % (please refer to our previous method for the lambda grid search method
      % at https://github.com/rkim35/spikeRNN)
      perf_measure = max(all_perfs);
    end

    if length(perf_thr) == 1
      if exist('taus_decay_ms') & perf_measure > perf_thr
        perfs = [perfs, perf_measure];
        all_trs = [all_trs, tr];
        all_taus = [all_taus, taus_decay_ms];
        all_tds = [all_tds, syn_decay(auto_N)'];
        syn_decay0 = (1./(1+exp(-taus_gaus0))*(taus(2) - taus(1))+taus(1))*5;
        all_tds0 = [all_tds0, syn_decay0(auto_N)'];
        all_auto_c = [all_auto_c; new_auto_c];
        counter = counter + 1;
        all_decays= [all_decays mean_decay];
        all_decays_mean = [all_decays_mean; nanmean(new_auto_c)];
        mean_auto_c = [mean_auto_c; nanmean(auto_c)];
        all_inh = [all_inh, inh(auto_N)'];
        all_frs = [all_frs, auto_N_fr'];
        all_scaling_factors = [all_scaling_factors, opt_scaling_factor];

        W = w*m.*som_m/opt_scaling_factor;
        w_inh = W(:, inh); w_inh = w_inh(:);
        w_inh = w_inh(w_inh~=0);
        w_inhs = [w_inhs, nanmean(w_inh)];

        w_exc = W(:, exc); w_exc = w_exc(:);
        w_exc = w_exc(w_exc~=0);
        w_excs = [w_excs, nanmean(w_exc)];

        w_ee = W(exc, exc); w_ee= w_ee(:);
        w_ee = w_ee(w_ee~=0);
        w_ees = [w_ees, nanmean(w_ee)];

        w_ii = W(inh, inh); w_ii= w_ii(:);
        w_ii = w_ii(w_ii~=0);
        w_iis = [w_iis, nanmean(w_ii)];

        w_ie = W(inh, exc); w_ie= w_ie(:);
        w_ie = w_ie(w_ie~=0);
        w_ies = [w_ies, nanmean(w_ie)];

        w_ei = W(exc, inh); w_ei= w_ei(:);
        w_ei = w_ei(w_ei~=0);
        w_eis = [w_eis, nanmean(w_ei)];

        % Number of connections by each synaptic type
        N_wii = W(inh, inh); N_wii = N_wii(:); N_wii = length(find(N_wii));
        N_wee = W(exc, exc); N_wee = N_wee(:); N_wee = length(find(N_wee));
        N_wei = W(exc, inh); N_wei = N_wei(:); N_wei = length(find(N_wei));
        N_wie = W(inh, exc); N_wie = N_wie(:); N_wie = length(find(N_wie));
        num_connections(1) = num_connections(1) + N_wii;
        num_connections(2) = num_connections(2) + N_wie;
        num_connections(3) = num_connections(3) + N_wei;
        num_connections(4) = num_connections(4) + N_wee;

        neu_ind = [neu_ind, auto_N];
        neu_tau_ind = [neu_tau_ind, taus_decay_ms];

        auto_spks = cell(1, length(taus_decay_ms));
        for ijk = 1:length(taus_decay_ms)
          auto_spks{ijk} = trial_spks(auto_N(ijk):N:end, :);
        end
        all_spks = [all_spks, auto_spks];

        fnames = [fnames, dirs(d).name];

        clearvars -except dirs model_dir task_type perfs perf_thr num_mods ...
        all_* mean_auto_c counter ns_* sf_* whi_* mante_* fr_lim pop_fit neu_* ...
        num_connections w_inhs w_excs w_iis w_ees w_eis w_ies fnames;
      end
    elseif length(perf_thr) == 2
      if exist('taus_decay_ms') & perf_measure > perf_thr(1) & perf_measure < perf_thr(2)
        perfs = [perfs, perf_measure];
        all_taus = [all_taus, taus_decay_ms];
        all_trs = [all_trs, tr];
        all_tds = [all_tds, syn_decay(auto_N)'];
        syn_decay0 = (1./(1+exp(-taus_gaus0))*(taus(2) - taus(1))+taus(1))*5;
        all_tds0 = [all_tds0, syn_decay0(auto_N)'];
        all_auto_c = [all_auto_c; new_auto_c];
        counter = counter + 1;
        all_decays= [all_decays mean_decay];
        all_decays_mean = [all_decays_mean; nanmean(new_auto_c)];
        mean_auto_c = [mean_auto_c; nanmean(auto_c)];
        all_inh = [all_inh, inh(auto_N)'];
        all_frs = [all_frs, auto_N_fr'];
        all_scaling_factors = [all_scaling_factors, opt_scaling_factor];

        W = w*m.*som_m/opt_scaling_factor;
        w_inh = W(:, inh); w_inh = w_inh(:);
        w_inh = w_inh(w_inh~=0);
        w_inhs = [w_inhs, nanmean(w_inh)];

        w_exc = W(:, exc); w_exc = w_exc(:);
        w_exc = w_exc(w_exc~=0);
        w_excs = [w_excs, nanmean(w_exc)];

        w_ee = W(exc, exc); w_ee= w_ee(:);
        w_ee = w_ee(w_ee~=0);
        w_ees = [w_ees, nanmean(w_ee)];

        w_ii = W(inh, inh); w_ii= w_ii(:);
        w_ii = w_ii(w_ii~=0);
        w_iis = [w_iis, nanmean(w_ii)];

        w_ie = W(inh, exc); w_ie= w_ie(:);
        w_ie = w_ie(w_ie~=0);
        w_ies = [w_ies, nanmean(w_ie)];

        w_ei = W(exc, inh); w_ei= w_ei(:);
        w_ei = w_ei(w_ei~=0);
        w_eis = [w_eis, nanmean(w_ei)];

        % Number of connections by each synaptic type
        N_wii = W(inh, inh); N_wii = N_wii(:); N_wii = length(find(N_wii));
        N_wee = W(exc, exc); N_wee = N_wee(:); N_wee = length(find(N_wee));
        N_wei = W(exc, inh); N_wei = N_wei(:); N_wei = length(find(N_wei));
        N_wie = W(inh, exc); N_wie = N_wie(:); N_wie = length(find(N_wie));
        num_connections(1) = num_connections(1) + N_wii;
        num_connections(2) = num_connections(2) + N_wie;
        num_connections(3) = num_connections(3) + N_wei;
        num_connections(4) = num_connections(4) + N_wee;

        neu_ind = [neu_ind, auto_N];
        neu_tau_ind = [neu_tau_ind, taus_decay_ms];

        auto_spks = cell(1, length(taus_decay_ms));
        for ijk = 1:length(taus_decay_ms)
          auto_spks{ijk} = trial_spks(auto_N(ijk):N:end, :);
        end
        all_spks = [all_spks, auto_spks];

        fnames = [fnames, dirs(d).name];

        clearvars -except dirs model_dir task_type perfs perf_thr num_mods ...
        all_* mean_auto_c counter ns_* sf_* whi_* mante_* fr_lim pop_fit neu_* ...
        num_connections w_inhs w_excs w_iis w_ees w_eis w_ies fnames;
      end
    end
    if ~isnan(num_mods) & counter == num_mods
      break
    end
  end
end

if ~isempty(fr_lim)
  fr_units = find(all_frs >= fr_lim(1) & all_frs <= fr_lim(2));
  all_taus = all_taus(fr_units);
  all_tds = all_tds(fr_units);
  all_auto_c = all_auto_c(fr_units, :);
  all_inh = all_inh(fr_units);
  all_frs = all_frs(fr_units);
  all_spks = all_spks(fr_units);
end

% Exponential decay fitting
if pop_fit == true
  if length(find(isnan(all_auto_c(:, end)))) == size(all_auto_c, 1)
    mean_auto = nanmean(all_auto_c(:, 1:end-1));
  else
    mean_auto = nanmean(all_auto_c);
  end
  xdata = 1:size(mean_auto, 2);
  ydata = mean_auto(1:end);
  fun = @(x, xdata)x(1)*(exp(x(2)*xdata)+x(3));
  x0 = [0, 0, 0];
  options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
  x = lsqcurvefit(fun, x0, xdata, ydata, [], [], options);
  xdata2 = linspace(xdata(1), xdata(end), length(xdata)*2);
  taus_decay = 1/-x(2);
  ydata2 = fun(x, xdata2);
end

out.taus = all_taus;
out.auto_c = all_auto_c;
out.tds = all_tds;
out.tds0 = all_tds0;
out.decays = all_decays;
out.decays_mean = all_decays_mean;
out.inh = all_inh;
if pop_fit == true
  out.fit_x = xdata2;
  out.fit_y = ydata2;
  out.pop_tau = taus_decay;
end
out.frs = all_frs;
out.trs = all_trs;
out.spks = all_spks;
out.perfs = perfs;
out.neu_ind = neu_ind;
out.neu_tau_ind = neu_tau_ind;
out.all_scaling_factors = all_scaling_factors;
out.num_connections = num_connections;
out.w_inhs = w_inhs;
out.w_excs = w_excs;
out.w_iis = w_iis;
out.w_ees = w_ees;
out.w_eis = w_eis;
out.w_ies = w_ies;
out.fnames = fnames;



