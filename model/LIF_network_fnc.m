% Name: Robert Kim
% Date: 10-22-2018
% Email: robert.f.kim@gmail.com
% LIF_network_RK.m
% Description: Script to construct a network of LIF units and
% transfer the trained weights from the continuous rate model
% NOTE: LIF network implementation taken from LIFFORCESINE.m from NicolaClopath2017
% NOTE: edited to add an option to remove auto-connections

function [W, REC, spk, rs, all_fr, out, params] = LIF_network_fnc(model_path,...
scaling_factor, u, stims, downsample, shuffle_weights, use_this_W, shuffle_w_in, syn_type, remove_autoconn)
% FUNCTION LIF_network_fnc
% INPUT
%   - model_path: trained model full path (directory + filename)
%   - scaling_factor: scaling factor for transferring weights from rate to spk
%   - u: input stimulus to be used
%   - stims: struct for artificial stimulations (to model optogenetic stim)
%       - mode: "none", "exc" (depolarizing), or "inh" (hyperpolarizing)
%       - dur: [stim_onset stim_offset]
%       - units: vector containing unit indices to be stimulated
%
% OUTPUT
%   - W: recurrent connectivity matrix scaled by the scaling factor (N x N)
%   - REC: membrane voltage from all the units (N x t)
%   - spk: binary matrix indicating spikes (N x t)
%   - rs: firing rates from all the units (N x t)
%   - all_fr: average firing rates from all the units (N x 1)
%   - out: network output (1 x t)
%   - params: struct containing sampling rate info

%------------------------------------------------------
% Extract the number of units and the connectivity
% matrix from the trained continuous rate model
%------------------------------------------------------

if nargin == 9
  remove_autoconn = false;
elseif nargin == 8
  syn_type = [];
  remove_autoconn = false;
elseif nargin  == 7
  shuffle_w_in = [];
  syn_type = [];
  remove_autoconn = false;
elseif nargin  == 6
  use_this_W = [];
  shuffle_w_in = [];
  syn_type = [];
  remove_autoconn = false;
end

%load(model_path, 'w_in', 'w', 'w0', 'N', 'm', 'som_m', 'w_out', ...
%'inh', 'exc');

load(model_path, 'w_in', 'w', 'w0', 'N', 'm', 'som_m', 'w_out', ...
'inh', 'exc', 'taus_gaus0', 'taus_gaus', 'taus', 'tau');

% Number of neurons and the trained connectivity weight
% matrix (extracted from the trained continuous rate model)
N = double(N); 
%w = w*m.*som_m;
%w = w0*m.*som_m;

if ~isempty(syn_type)
  shuffle_weights = false;
end

if strcmpi(syn_type, 'ee_inc')
  w(exc, exc) = w(exc, exc)*1.3;
elseif strcmpi(syn_type, 'ee_dec')
  w(exc, exc) = w(exc, exc)*0.7;
elseif strcmpi(syn_type, 'ii_inc')
  w(inh, inh) = w(inh, inh)*1.3;
elseif strcmpi(syn_type, 'ii_dec')
  w(inh, inh) = w(inh, inh)*0.7;
elseif strcmpi(syn_type, 'ie_inc')
  w(inh, exc) = w(inh, exc)*1.3;
elseif strcmpi(syn_type, 'ie_dec')
  w(inh, exc) = w(inh, exc)*0.7;
elseif strcmpi(syn_type, 'ei_inc')
  w(exc, inh) = w(exc, inh)*1.3;
elseif strcmpi(syn_type, 'ei_dec')
  w(exc, inh) = w(exc, inh)*0.7;
end


% Shuffle nonzero weights
if shuffle_weights == true
  % Initial random weights
  %w0(inh, inh) = w0(inh, inh)*2;
  %w = w0*m.*som_m;

  % **Random weights**
%{
  spr_w = sprandn(N, N, 0.05);
  spr_w = full(spr_w);
  spr_w = spr_w/sqrt(N*0.05)*1.5;
  spr_w = max(0, spr_w);
  %spr_w(inh, inh) = spr_w(inh, inh)*4;
  w = spr_w*m.*som_m;
%}

%{
  % Random with dense I -> E
  spr_w_ei = sprandn(length(find(exc)), length(find(inh)), 0.35);
  %spr_w_ei = sprandn(length(find(exc)), length(find(inh)), 0.05);
  spr_w_ei = full(spr_w_ei);

  %spr_w_ii = sprandn(length(find(inh)), length(find(inh)), 0.05);
  spr_w_ii = sprandn(length(find(inh)), length(find(inh)), 0.25);
  spr_w_ii = full(spr_w_ii);

  %spr_w_e = sprandn(N, length(find(exc)), 0.10);
  spr_w_e = sprandn(N, length(find(exc)), 0.003);
  spr_w_e = full(spr_w_e);

  w = zeros(N, N);
  w(exc, inh) = spr_w_ei;
  w(inh, inh) = spr_w_ii;
  w(:, exc) = spr_w_e;
  w = max(0, w);
  %w(inh, inh) = w(inh, inh)*4;
  w = w*m.*som_m;
%}

  % Diminish inh unit weights
  w(inh, inh) = w(inh, inh)*0.5;
  w = w*m.*som_m;

  % Diminish exc unit weights
  %w(:, exc) = w(:, exc)/4;
  %w(:, exc) = 0;
  %w = w*m.*som_m;

  % Diminish only the strongly inhibited units
  %w_tmp = w*m.*som_m;
  %indx = fnc_find_most_inh(w_tmp, inh);
  %w(indx, :) = w(indx, :)/2;
  %w = w*m.*som_m;


%{
  % Create random E->E connections
  exc_ind = find(exc); inh_ind = find(inh);
  w_ee = w(exc, exc); w_ee = w_ee(:);
  w_ee = w_ee(w_ee ~= 0);
  w_ie = w(inh, exc); w_ie = w_ie(:);
  w_ie = w_ie(w_ie ~= 0);
  for nn = 1:length(exc_ind)
    for nnn = 1:length(inh_ind)
      if rand < 0.04
        %w(exc_ind(nnn), exc_ind(nn)) = nanmean(w_ee)*2;
        w(inh_ind(nnn), exc_ind(nn)) = nanmean(w_ie)*2;
      end
    end
  end
  w = w*m.*som_m;
%}

%{
  % Create random E->all connections
  exc_ind = find(exc); inh_ind = find(inh);
  w_ee = w(:, exc); w_ee = w_ee(:);
  w_ee = w_ee(w_ee ~= 0);
  for nn = 1:length(exc_ind)
    for nnn = 1:N
      if rand < 0.03
        w(nnn, exc_ind(nn)) = nanmean(w_ee)*2;
      end
    end
  end
  w = w*m.*som_m;
%}

  %nw = w(:);
  %nw = nw(nw ~= 0);
  %w = zeros(N*N, 1);
  %rand_perm = randperm(N*N);
  %w(rand_perm(1:length(nw))) = nw;
  %w = reshape(w, N, N);
  %w = w*m.*som_m;

%{
  % Use this for shuffling
  w_original = w*m.*som_m;
  inhN = find(inh);
  excN = find(exc);
  rand_exc = randperm(length(excN));
  for nn = 1:N
  %for nn = 1:length(excN)
  %for nn = 1:round(length(excN)/2)
    curr_w_col = w(:, nn);
    %curr_w_col = w(:, excN(nn));
    %curr_w_col = w(:, excN(rand_exc(nn))); % for shuffling half of exc
    next_w_col = zeros(N, 1);
    curr_w_col_nonzero = curr_w_col(find(curr_w_col ~= 0));
    rand_p = randperm(N);
    next_w_col(rand_p(1:length(curr_w_col_nonzero))) = curr_w_col_nonzero;
    %w(:, excN(rand_exc(nn))) = next_w_col;
    %w(:, excN(nn)) = next_w_col;
    w(:, nn) = next_w_col;
  end
  w = w*m.*som_m;
  w(inh, inh) = w_original(inh, inh);
  w(inh, inh) = w(inh, inh)*2;
  %w(exc, inh) = w_original(exc, inh);
  %w(inh, exc) = w_original(inh, exc);
%}

%{
  % INH -> INH/EXC shuffling
  inhN = find(inh);
  excN = find(exc);
  for nn = 1:length(inhN)
  %for nn = 1:round(length(inhN)/5)
    curr_w_col = w(:, inhN(nn));
    next_w_col = curr_w_col;

    % INH -> EXC
    next_w_col(excN) = 0;
    curr_w_col_nonzero = curr_w_col(find(curr_w_col ~= 0 & exc));
    rand_p = randperm(length(excN));
    next_w_col(excN(rand_p(1:length(curr_w_col_nonzero)))) = curr_w_col_nonzero;

    % INH -> INH
    %next_w_col(inhN) = 0;
    %curr_w_col_nonzero = curr_w_col(find(curr_w_col ~= 0 & inh));
    %rand_p = randperm(length(inhN));
    %next_w_col(inhN(rand_p(1:length(curr_w_col_nonzero)))) = curr_w_col_nonzero;

    w(:, inhN(nn)) = next_w_col;
  end
  w = w*m.*som_m;
%}

%{
  % EXC -> INH/EXC shuffling
  inhN = find(inh);
  excN = find(exc);
  for nn = 1:length(excN)
    curr_w_col = w(:, excN(nn));
    next_w_col = curr_w_col;

    % EXC -> EXC
    next_w_col(excN) = 0;
    curr_w_col_nonzero = curr_w_col(find(curr_w_col ~= 0 & exc));
    rand_p = randperm(length(excN));
    next_w_col(excN(rand_p(1:length(curr_w_col_nonzero)))) = curr_w_col_nonzero;

    % EXC -> INH
    %next_w_col(inhN) = 0;
    %curr_w_col_nonzero = curr_w_col(find(curr_w_col ~= 0 & inh));
    %rand_p = randperm(length(inhN));
    %next_w_col(inhN(rand_p(1:length(curr_w_col_nonzero)))) = curr_w_col_nonzero;

    w(:, excN(nn)) = next_w_col;
  end
  w = w*m.*som_m;
%}

%{
  % Use this for removing reciprocal connections
  w = w*m.*som_m;
  tmp_w = w.*w';
  for qq = 1:N
    tmp_ind = find(tmp_w(qq, :) ~= 0);
    if ~isempty(tmp_ind)
      for qqq = 1:length(tmp_ind)
        if tmp_ind(qqq) > 0
          w(qq, qqq) = 0;
        end
      end
    end
    %w(qq, tmp_ind) = 0;
  end
%}

%{
    % Similar tunes
    disinh_strength = 0.7;
    f_st = strfind(model_path, '/');
    f_path = model_path(1:f_st(end));
    f_name = model_path(f_st(end)+1:end);
    load(fullfile(f_path, '..', 'selectivity', f_name), 'inh_pos_ind', 'inh_neg_ind');

    w(inh_pos_ind, inh_pos_ind) = w(inh_pos_ind, inh_pos_ind)*disinh_strength;
    w(inh_neg_ind, inh_neg_ind) = w(inh_neg_ind, inh_neg_ind)*disinh_strength;
    w = w*m.*som_m;
%}

%{
    % Dissimilar tunes
    disinh_strength = 0.7;
    f_st = strfind(model_path, '/');
    f_path = model_path(1:f_st(end));
    f_name = model_path(f_st(end)+1:end);
    load(fullfile(f_path, '..', 'selectivity', f_name), 'inh_pos_ind', 'inh_neg_ind');

    w(inh_pos_ind, inh_neg_ind) = w(inh_pos_ind, inh_neg_ind)*disinh_strength;
    w(inh_neg_ind, inh_pos_ind) = w(inh_neg_ind, inh_pos_ind)*disinh_strength;
    w = w*m.*som_m;
%}

else
  w = w*m.*som_m;
end

%w00 = w0/sqrt(N*0.20)*1.5;
%w = w00*m.*som_m;
%w = w0;

% Remove auto-connections (i.e. self connections)
if remove_autoconn == true
  for ii = 1:size(w, 1)
    w(ii, ii) = 0;
  end
end


if isempty(use_this_W)
  OMEGA = w/scaling_factor;
else
  OMEGA = use_this_W/scaling_factor;
end

%for id = 1:N
  %OMEGA(id, id) = 0;
%end


%OMEGA(dspns, soms) = 0;
%OMEGA(chats, soms) = 0;


% Inhibitory and excitatory neurons
inh_ind = find(inh);
exc_ind = find(exc);

% Altering connectivity matrix
%OMEGA(inh_ind, exc_ind) = OMEGA(inh_ind, exc_ind)*0.10;
%OMEGA(exc_ind, inh_ind) = OMEGA(exc_ind, inh_ind)*0.15;
%OMEGA(exc_ind, exc_ind) = OMEGA(exc_ind, exc_ind)*0.50;
%OMEGA(inh_ind, inh_ind) = OMEGA(inh_ind, inh_ind)*0.10;
%OMEGA(:, inh_ind) = OMEGA(:, inh_ind)*0.135;
%OMEGA(:, inh_ind) = OMEGA(:, inh_ind)*0.15;
%rand_inh = randperm(length(inh_ind));
%OMEGA(:, inh_ind(rand_inh(1:25))) = OMEGA(:, inh_ind(rand_inh(1:25)))*0.15;
%OMEGA(:, inh_ind) = OMEGA(:, inh_ind)*5;
%OMEGA(:, exc_ind) = OMEGA(:, exc_ind)*0.15;
%OMEGA(go_pref_units, :) = OMEGA(go_pref_units, :)*0;
%OMEGA(nogo_pref_units, :) = OMEGA(nogo_pref_units, :)*0;
%OMEGA(go_pref_units, :) = 0.05;
%OMEGA(nogo_pref_units, :) = 0.05;
%OMEGA(:, go_pref_units) = OMEGA(:, go_pref_units)*1.50;
%OMEGA(:, nogo_pref_units) = OMEGA(:, nogo_pref_units)*1.50;

%som_inh_ind = inh_ind(1:som_N);
%OMEGA(exc_ind, som_inh_ind) = OMEGA(exc_ind, som_inh_ind)*0.50;
%OMEGA(:, som_inh_ind) = OMEGA(:, som_inh_ind)*0.75;
%OMEGA(som_inh_ind, exc_ind) = OMEGA(som_inh_ind, exc_ind)*0;

W = OMEGA;

%for nn = 1:N
  %for mm = 1:N
    %if OMEGA(nn, mm) > -0.01
      %OMEGA(nn, mm) = 0;
    %end
  %end
%end


% Shuffling the input weights
if ~isempty(shuffle_w_in)
  w_in = shuffle_w_in;
end


%u = u(:, 1:4:end);
u = u(:, 1:downsample:end);
%w_in(inh, 1) = w_in(inh, 1)*0.20;
%w_in(inh, 2) = w_in(inh, 2)*0.20;
%w_in(exc, :) = w_in(exc, :)*5.50;
%w_in = [w_in, w_in];
ext_stim = w_in*u;
%ext_stim = w_in*u*3;
%ext_stim = w_in*zeros(size(u));
%ext_stim = w_in*rand(size(u));

%------------------------------------------------------
% LIF network parameters
%------------------------------------------------------
dt = 0.00005*downsample;   % sampling rate
%dt = 0.00005*4;   % sampling rate
T = (size(u, 2)-1)*dt*100;
nt = round(T/dt);
tref = 0.002;   % refractory time constant (in sec)
tm = 0.010;      % membrane time constant (in sec)
vreset = -65;   % voltage reset (in mV)
vpeak = -40;    % voltage peak (in mV) for linear LIF
%vpeak = 30;    % voltage peak (in mV) for quadratic LIF
%rng(1);
%td = 0.02; tr = 0.005; % synaptic filters ORIGINAL
%td = 0.18; tr = 0.005; % synaptic filters for SINE
%td = 0.125; tr = 0.002; % synaptic filters for XOR 

%td = ones(N, 1)*0.125; tr = 0.002;

%td = (1./(1+exp(-taus))*double(tau)+4)*5/1000; tr = 0.002;

if length(taus) > 1
    td = (1./(1+exp(-taus_gaus))*(taus(2) - taus(1))+taus(1))*5/1000; 
    td0 = (1./(1+exp(-taus_gaus0))*(taus(2) - taus(1))+taus(1))*5/1000;
    tr = 0.002;
else
    td = taus*5/1000; 
    td0 = td;
    tr = 0.002;
end

%tr = 0;

%td(inh) = 0.020;
%td(exc) = 0.020;
%td(:) = 0.020;
%td(:) = 0.120;

%td = td0;
%td = scramble_list(td); td = td';
%td = rand(N, 1)/10 + 20;
%td = ones(N, 1)*30/1000;

%td = 0.050; tr = 0.001; % synaptic filters for XOR 
%td = 0.080; tr = 0.001; % synaptic filters for XOR 
%td_inh = 0.14; tr_inh = 0.005; % synaptic filters

% Synaptic parameters
IPSC = zeros(N,1);      % post synaptic current storage variable
h = zeros(N,1);         % storage variable for filtered firing rates
r = zeros(N,1);         % second storage variable for filtered rates
hr = zeros(N,1);        % third variable for filtered rates
JD = 0*IPSC;            % storage variable required for each spike time
tspike = zeros(4*nt,2); % storage variable for spike times
ns = 0;                 % number of spikes, counts during simulation

v = vreset + rand(N,1)*(30-vreset); % initialize voltage with random distribtuions
%v = vreset;
v_ = v;   % v_ is the voltage at previous time steps
v0 = v;   % store the initial voltage values

% Record REC (membrane voltage), Is (input currents), 
% spk (spike raster), rs (firing rates) from all the units
REC = zeros(nt,N);  % membrane voltage (in mV) values
Is = zeros(N, nt);  % input currents from the ext_stim
IPSCs = zeros(N, nt); % IPSC over time
%IPSCs = zeros(N, 1); % IPSC over time
spk = zeros(N, nt); % spikes
rs = zeros(N, nt);  % firing rates
hs = zeros(N, nt); % filtered firing rates

% used to set the refractory times
tlast = zeros(N,1); 

% set the BIAS current, can help decrease/increase firing rates. 0 is fine.
BIAS = vpeak; % for linear LIF
%BIAS = 0; % for quadratic LIF
%BIAS = vreset;

%OMEGA(depol, :) = OMEGA(depol, :)*0;
%OMEGA(:, depol) = OMEGA(:, depol)*1.5;
%OMEGA(depol, 190) = 0;

%pv = scramble_list([1:N]);
%pv = pv(1:9);
%OMEGA(pv, depol) = 0;
%OMEGA(:, input_pv_chat) = 0;
%OMEGA(pv, input_pv_chat) = 0;
%OMEGA(pv, depol) = OMEGA(pv, depol)*3.5;
%OMEGA(inhibit_pvs, pv) = 0;
%OMEGA(depol, inhibit_pvs) = 0;
%OMEGA(inhibit_pvs, depol) = 0;
%OMEGA(gaba, depol) = 0;
%OMEGA(inhibit_pvs, gaba) = 0;
%OMEGA(pv, inhibit_pvs) = 0;

%------------------------------------------------------
% Start the simulation
%------------------------------------------------------
for i = 1:nt
    IPSCs(:, i) = IPSC; % record the IPSC over time (comment out if not used to save time)

    I = IPSC + BIAS; % synaptic current

    % Apply external input stim if there is any
    I = I + ext_stim(:, round(i/100)+1);
    Is(:, i) = ext_stim(:, round(i/100)+1);

    % LIF voltage equation with refractory period
    dv = (dt*i>tlast + tref).*(-v+I)/tm; % linear LIF
    %dv = (dt*i>tlast + tref).*(v.^2+I)/tm; % quadratic LIF
    %v = v + dt*(dv) + randn(N, 1)/4;
    v = v + dt*(dv) + randn(N, 1)/10;
    %v = v + dt*(dv) + randn(N, 1)/100;
    %v = v + dt*(dv);

    %if i > 400 & i < 440
      %v(314) = -42;
    %end
    %if i > 1121 & i < 1159
      %v(427) = -42;
    %end

    %v(find(td*1000 > 90)) = v(find(td*1000 > 90)) - 0.5;

    % Artificial stimulation/inhibition
    if strcmpi(stims.mode, 'exc')
      if i >= stims.dur(1) & i < stims.dur(2)
        if rand < 0.50
          v(stims.units) = v(stims.units) + 0.5;
        end
      end
    elseif strcmpi(stims.mode, 'inh')
      if i >= stims.dur(1) & i < stims.dur(2)
        if rand < 0.50
          v(stims.units) = v(stims.units) - 0.5;
        end
      end
    end

    %if i > 2000 & i < 3000
    %%if i > 2000 & i < 10000
      %if rand < 0.70
        %v([depol]) = ...
        %v([depol]) + 0.5;
      %end
      %%v([go_pref_units, nogo_pref_units]) = -44;
    %end

    %aaa = 190;
    %aaa = find(OMEGA(25, :) < 0);
    %%aaa = depol;
    %aaa = pv;
    %aaa = inhibit_pvs;
    %if i > 2000 & i < 4400
    %if i > 2000 & i < 18000
      %if rand < 0.10
        %v([aaa]) = ...
        %v([aaa]) + 0.5;
      %end
      %%v([go_pref_units, nogo_pref_units]) = -44;
    %end

    %aaa = gaba;
    %aaa = inhibit_pvs;
    %if i > 2000 & i < 4400
    %%if i > 2000 & i < 18000
      %if rand < 0.30
        %v([aaa]) = ...
        %v([aaa]) + 0.5;
      %end
      %%v([go_pref_units, nogo_pref_units]) = -44;
    %end

    %aaa = pv;
    %if i > 5400 & i < 5600
      %if rand < 0.70
        %v([aaa]) = ...
        %v([aaa]) + 0.5;
      %end
    %end

    % find the neurons that have spiked
    index = find(v>=vpeak);  

    %!!!!!!!! Raise the threshold for INH neurons !!!!!!!!!!
    %v_exc = v(exc);
    %v_inh = v(inh);
    %index_exc = find(v_exc>=vpeak);  
    %if rand > 0.1
      %index_inh = find(v_inh>=vpeak+0.15);
    %else
      %index_inh = find(v_inh>=vpeak);
    %end
    %index = [exc_ind(index_exc); inh_ind(index_inh)];

    %!!!!!!!! Raise the threshold for EXC neurons !!!!!!!!!!
    %v_exc = v(exc);
    %v_inh = v(inh);
    %index_exc = find(v_exc>=vpeak+0.15);  
    %index_inh = find(v_inh>=vpeak);
    %index = [exc_ind(index_exc); inh_ind(index_inh)];

    %% subset of inhibitory
    %sub_inh = inh_ind(1:20);
    %sub_inh = scramble_list(inh_ind);
    %sub_inh = sub_inh(1:20);
    %rm_indx = [];
    %if length(index) > 0
      %for ij = 1:length(index)
        %if find(index(ij) == sub_inh)
          %rm_indx = [rm_indx, ij];
        %end
      %end
    %end
    %index(rm_indx) = [];

    % Activating random inhibitory neurons
    %if i > 5000 & i < 6000
      %if rand > 0.85
        %rand_inh = scramble_list(inh_ind);
        %%index = [index; rand_inh(1:2)'];
        %index = [index; inh_ind(2)'];
        %index = unique(index);
      %end
    %end

    % Activating random excitatory neurons
    %if i > 5000 & i < 6000
      %if rand > 0.93
        %rand_exc = scramble_list(exc_ind);
        %%index = [index; rand_exc(1:2)'];
        %index = [index; exc_ind(2)'];
        %index = unique(index);
      %end
    %end

    %if i > 10000 & i < 30000
      %rm_indx = [];
      %if length(index) > 0
        %for ii = 1:length(index)
          %if length(find(index(ii) == exc_ind)) > 0
            %rm_indx = [rm_indx, ii];
          %end
        %end
      %end
      %index(rm_indx) = [];
    %end

    % store spike times, and get the weight matrix column sum of spikers
    if length(index)>0
      %if length(find(exc(index) == 1)) == length(index)
      %if length(find(index == 21)) == 0  
        JD = sum(OMEGA(:,index),2); %compute the increase in current due to spiking
        tspike(ns+1:ns+length(index),:) = [index, [0*index+dt*i]];
        ns = ns + length(index);  % total number of psikes so far
      %end
    end

    % used to set the refractory period of LIF neurons
    tlast = tlast + (dt*i -tlast).*(v>=vpeak);      

    % if the rise time is 0, then use the single synaptic filter,
    % otherwise (i.e. rise time is positive) use the double filter
    if tr == 0
        IPSC = IPSC.*exp(-dt./td)+JD*(length(index)>0)./(td);
        r = r.*exp(-dt./td) + (v>=vpeak)./td;
        rs(:, i) = r;
    else
        %IPSC = IPSC*exp(-dt/td) + h*dt;
        %h = h*exp(-dt/tr) + JD*(length(index)>0)/(tr*td);  %Integrate the current
        %hs(:, i) = h;
        %h_exc = h.*exc;

        IPSC = IPSC.*exp(-dt./td) + h*dt;
        h = h*exp(-dt/tr) + JD*(length(index)>0)./(tr*td);  %Integrate the current
        hs(:, i) = h;

        %h = h*exp(-dt/td_inh) + JD*(length(index)>0)/(tr_inh*td_inh);  %Integrate the current
        %h_inh = h.*inh;

        %h = h_exc + h_inh;
        
        %r = r*exp(-dt/td) + hr*dt;
        r = r.*exp(-dt./td) + hr*dt;
        %r_exc = r.*exc;

        %r = r*exp(-dt/tr_inh) + hr*dt;
        %r_inh = r.*inh;
        
        %r = r_exc + r_inh;

        %hr = hr*exp(-dt/tr) + (v>=vpeak)/(tr*td);
        hr = hr*exp(-dt/tr) + (v>=vpeak)./(tr.*td);
        %hr_exc  = hr.*exc;

        %hr = hr*exp(-dt/td_inh) + (v>=vpeak)/(tr_inh*td_inh);
        %hr_inh  = hr.*inh;

        %hr = hr_exc + hr_inh;

        rs(:, i) = r;
    end

    % record the spikes
    spk(:, i) = v >= vpeak;
    
    v = v + (30 - v).*(v>=vpeak);

    % record the membrane voltage tracings from all the units
    REC(i,:) = v(1:N); 

    %reset with spike time interpolant implemented.
    v = v + (vreset - v).*(v>=vpeak); 
end

time = 1:1:nt;
%TotNumSpikes = ns
%%tspike = tspike(1:1:ns,:);
%M = tspike(tspike(:,2)>0, :);
%AverageRate = length(M)/(N*T)

% Plot the population response
out = w_out/scaling_factor*rs;

% "Essential" neurons
%essential_neus = unique([depol, pv, inhibit_pvs, gaba]);
%out2 = w_out(essential_neus)*rs(essential_neus, :)/40;

%figure; axis tight; hold on;
%plot(pop_response); 
%%plot(1:100:nt, target);
%plot(0:100:nt, u);


% Compute average firing rate for each population (excitatory/inhibitory)
inh_fr = zeros(size(inh_ind));
for i = 1:length(inh_ind)
  inh_fr(i) = length(find(spk(inh_ind(i), :)>0))/T;
end

exc_fr = zeros(size(exc_ind));
for i = 1:length(exc_ind)
  exc_fr(i) = length(find(spk(exc_ind(i), :)>0))/T;
end

all_fr = zeros(size(exc));
for i = 1:N
  all_fr(i) = length(find(spk(i, :)>0))/T;
end

REC = REC';

params = {};
params.dt =  dt;
params.T = T;
params.nt = nt;
params.w_out = w_out;
params.td = td;
params.td0 = td0;
params.IPSCs = IPSCs;



