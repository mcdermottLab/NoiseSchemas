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

azimuth_errors = nan(length(results(1).mean_azimuth_error),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)
    azimuth_errors(:,i) = results(i).mean_azimuth_error;
    keep_subs(i) = results(i).pass;
end

azimuth_errors = azimuth_errors(:, logical(keep_subs));
n_subs = size(azimuth_errors,2);
n_conditions = 5;

datatable = array2table(azimuth_errors');
within_table = table(categorical(1:n_conditions)', 'VariableNames', {'Time'});
rm = fitrm(datatable,'Var1-Var5~1','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'Time');
report_main_effect(rm, C{2}, ranovatbl, 'Time')

x = linspace(250,2250,5);
y = mean(azimuth_errors,2)';
est_params_avg = fit_elbow_func_localization(x,y);
actual_elbow_point = est_params_avg(3);
elbow_points = est_params(:,3);

fprintf('Elbow point: %.0f ms, 95%% CI: [%.0f, %.0f] ms\n', actual_elbow_point, prctile(elbow_points, 2.5), prctile(elbow_points, 97.5))