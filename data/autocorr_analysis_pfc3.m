% Name: Robert Kim
% Date: 01-19-2020
% Email: rkim@salk.edu
% autocorr_analysis_pfc3.m
% Description: Script to compute neuronal timescales from the PFC3 dataset
% NOTE: make sure the 'SummaryDatabase.xlsx' (included in the dataset) is
% included in the output directory (out_dir).

clear; clc;

data_dir = DATA DIRECTORY;
out_dir = OUTPUT DIRECTORY

% learning condition ('pre', 'post')
conds = {'pre', 'post'};

% task type ('feature', 'spatial')
types = {'feature', 'spatial'};

% prefrontal cortex area ('dorsal', 'ventral')
areas = {'dorsal'};


for aaa = 1:length(conds)
  for bbb = 1:length(types)
    for ccc = 1:length(areas)
      out_fname = ['autocorr_' conds{aaa} '_train_' types{bbb} '_' areas{ccc} '.mat'];

      xls_file = fullfile(out_dir, 'SummaryDatabase.xlsx');

      % Read the xls file
      if strcmpi(conds{aaa}, 'pre') & strcmpi(types{bbb}, 'feature')
        [xls_num, xls_txt] = xlsread(xls_file, 'Pre_FeatureSets');
      elseif strcmpi(conds{aaa}, 'pre') & strcmpi(types{bbb}, 'spatial')
        [xls_num, xls_txt] = xlsread(xls_file, 'Pre_SpatialSet');
      elseif strcmpi(conds{aaa}, 'post') & strcmpi(types{bbb}, 'feature')
        [xls_num, xls_txt] = xlsread(xls_file, 'Post_FeatureSets');
      elseif strcmpi(conds{aaa}, 'post') & strcmpi(types{bbb}, 'spatial')
        [xls_num, xls_txt] = xlsread(xls_file, 'Post_SpatialSet');
      end

      post_txt = xls_txt(:, 1);
      area_txt = xls_txt(:, 5);
      cell_txt = xls_num;

      if ~exist(fullfile(out_dir, out_fname))
        mat_files = dir(fullfile(data_dir, '*.mat'));

        counter = 1;
        fr = [];
        neu_spk = {};
        neu_nums = [];
        spk_t = 0:50:1000;
        for i = 1:length(mat_files) % for each neuron
          i
          curr_fname = mat_files(i).name;
          cell_num = findstr(curr_fname, '.mat');
          cell_num = curr_fname(cell_num-4:cell_num-1);
          if isempty(str2num(cell_num))
            clearvars -except data_dir out_dir mat_files counter fr...
            spk_t neu_spk xls_file post_txt out_fname cell_txt ...
            area_txt neu_nums conds areas types aaa bbb ccc;
            continue;
          end

          if ~isempty(find(contains(post_txt, curr_fname(1:8))))
            curr_area = area_txt(find(cell_txt == str2num(cell_num)));
            if strcmpi(areas{ccc}, 'dorsal')
              target_area = 'ventral';
            else
              target_area = 'dorsal';
            end

            if strcmpi(curr_area, target_area)
              clearvars -except data_dir out_dir mat_files counter fr...
              spk_t neu_spk xls_file post_txt out_fname cell_txt ...
              area_txt neu_nums conds areas types aaa bbb ccc;
              continue;
            else
              load(fullfile(mat_files(i).folder, mat_files(i).name));
              if ~exist('MatData')
                clearvars -except data_dir out_dir mat_files counter fr...
                spk_t neu_spk xls_file post_txt out_fname cell_txt ...
                area_txt neu_nums conds areas types aaa bbb ccc;
                continue
              end
              if ~isempty(MatData)
                spk = [];
                % for each trial type (8 for features and 9 for spatial)
                for tr = 1:length(MatData.class) 
                  curr_data = MatData.class(tr).ntr;
                  % for each trial
                  for j = 1:length(curr_data) 
                    curr_cue = curr_data(j).Cue_onT;
                    curr_TS = curr_data(j).TS;
                    curr_fix = curr_data(j).fix; % fixation firing rate
                    if curr_cue >= 1 & curr_fix > 1
                      curr_preTS = curr_TS(curr_TS <= curr_cue)*1000;
                      [v, b] = hist(curr_preTS, spk_t);
                      spk = [spk; v];
                    end
                  end
                end
                neu_spk{counter} = spk;
                neu_nums = [neu_nums, str2num(cell_num)];
                counter = counter + 1;
                clearvars -except data_dir out_dir mat_files counter fr...
                spk_t neu_spk xls_file post_txt out_fname cell_txt ...
                area_txt neu_nums conds areas types aaa bbb ccc;
              end
              clearvars -except data_dir out_dir mat_files counter fr...
              spk_t neu_spk xls_file post_txt out_fname cell_txt ...
              area_txt neu_nums conds areas types aaa bbb ccc;
            end
          end
        end

        % Compute autocorrelation
        auto_c = nan(length(neu_spk), 15);
        for ij = 1:length(neu_spk)
            ij
            ex_neu = neu_spk{ij};
            if size(ex_neu, 1) > 10 % minimum number of trials for each neuron
              for jjj = 1:size(auto_c, 2)
                  curr_c = [];
                  for iii = 1:size(ex_neu, 2) - jjj
                      curr_corr = corr(ex_neu(:, iii), ex_neu(:, iii+jjj));
                      curr_c = [curr_c, curr_corr];
                  end
                  auto_c(ij, jjj) = nanmean(curr_c);
              end
            else
              continue;
            end
        end

        save(fullfile(out_dir, out_fname), '-v7.3');
      else
        load(fullfile(out_dir, out_fname));
      end

      new_auto_c = nan(size(auto_c));

      % Compute the taus
      if ~exist('taus')
        taus = nan(size(neu_spk));
        taus_amp = nan(size(neu_spk));
        for n = 1:length(neu_spk)
          n
          xdata = 1:size(auto_c, 2);
          ydata = auto_c(n, :);
          %xdata = 1:size(auto_c, 2)-1;
          %ydata = auto_c(n, 2:end);

          mm = 2;
          delta = 0;
          while delta >= 0
            delta = ydata(mm) - ydata(mm-1);
            mm = mm + 1;
          end
          mm = mm - 2;

          if length(find(isnan(ydata))) == 0 & mm < 4
            xdata = 1:size(auto_c, 2) - (mm-1);
            ydata = auto_c(n, mm:end);

            new_auto_c(n, xdata) = ydata;

            fun = @(x, xdata)x(1)*(exp(x(2)*xdata)+x(3));
            x0 = [0, 0, 0];
            options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
            x = lsqcurvefit(fun, x0, xdata, ydata, [], [], options);
            taus(n) = -1/x(2);
            taus_amp(n) = x(1);
          end
        end
        taus = taus*50;
        save(fullfile(out_dir, out_fname), 'new_auto_c', 'taus', 'taus_amp', '-append');
      end
    end
  end
end


%figure; hold on; axis tight;
%plot(nanmean(auto_c(taus > 0 & taus <=500, 2:end)), 'ro-', 'linewidth', 2,...
%'markers', 8, 'MarkerFace', 'w');

%figure; hold on; axis tight;
%plot(nanmean(auto_c(taus >= 0 & taus <=200, :)), 'ro-', 'linewidth', 2,...
%'markers', 8, 'MarkerFace', 'w');


