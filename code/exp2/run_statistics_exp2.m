[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
timing_analysis_file = fullfile(file_path, sprintf('../../results/%s/timing_analysis.mat', exp_name)); % get path to results
load(timing_analysis_file)

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

percent_correct = nan(length(results(1).percent_correct),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)
    percent_correct(:,i) = results(i).percent_correct;
    keep_subs(i) = results(i).pass;
end

percent_correct = percent_correct(:, logical(keep_subs));
n_subs = size(percent_correct,2);

datatable = array2table(percent_correct');
group_info = results(1).group_info;
within_table = table(categorical(group_info(:,1)), categorical(group_info(:,2)), 'VariableNames', {'SNR', 'Time'});
rm = fitrm(datatable,'Var1-Var20~1','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'SNR*Time'); % C is used to handle epsilon correction for violations of sphericity

report_main_effect(rm, C{2}, ranovatbl, 'SNR')
report_main_effect(rm, C{3}, ranovatbl, 'Time')
report_interaction(rm, C{4}, ranovatbl, 'SNR', 'Time')

percent_correct_avg = (percent_correct(1:10,:) + percent_correct(11:20,:)) / 2;
x = linspace(250,2500,10);
y = mean(percent_correct_avg,2)';
est_params_avg = fit_elbow_func(x,y);
actual_elbow_point = est_params_avg(3);
elbow_points = est_params(:,3);

fprintf('Elbow point: %.0f ms, 95%% CI: [%.0f, %.0f] ms\n', actual_elbow_point, prctile(elbow_points, 2.5), prctile(elbow_points, 97.5))