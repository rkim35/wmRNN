% Name: Robert Kim
% Date: 07-10-2019
% Email: rkim@salk.edu
% return_stable.m
% Description: Function to return a list of stable RNNs

function [fnames] = return_stable(task_dir, wcard, perf_threshold, task_type, max_tr)
% INPUT
%   task_dir: directory for the trained models
%   wcard: wildcard for getting files
%   perf_threshold: performance threshold
%   task_type: either 'xor' or 'non-xor'
%   max_tr: max number of trials

if nargin < 5
  max_tr = 5999;
end

mat_files = dir(fullfile(task_dir, wcard));

fnames = {};
for i = 1:length(mat_files)
  if ~isdir(fullfile(task_dir, mat_files(i).name)) & ~isempty(strfind(mat_files(i).name, '.mat'))
    curr_mat = fullfile(task_dir, mat_files(i).name);

    if strcmpi(task_type, 'xor')
      load(curr_mat, 'stable_perfs', 'tr');
      perf_measure = mean(stable_perfs);
    else
      %load(curr_mat, 'all_perfs', 'tr');
      load(curr_mat);
      if isnan(opt_scaling_factor)
        continue;
      end
      perf_measure = max(all_perfs);
    end

    if length(perf_threshold) > 1
      if perf_measure > perf_threshold(1) & perf_measure < perf_threshold(2) & tr < max_tr
        fnames = [fnames, mat_files(i).name];
        clearvars -except mat_files fnames task_dir wcard perf_threshold task_type max_tr;
      else
        clearvars -except mat_files fnames task_dir wcard perf_threshold task_type max_tr;
      end

    elseif length(perf_threshold) == 1
      if perf_measure > perf_threshold & tr < max_tr
        fnames = [fnames, mat_files(i).name];
        clearvars -except mat_files fnames task_dir wcard perf_threshold task_type max_tr;
      else
        clearvars -except mat_files fnames task_dir wcard perf_threshold task_type max_tr;
      end
    end

  end
end


