[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
timing_analysis_file = fullfile(file_path, sprintf('../../results/%s/timing_analysis.mat', exp_name)); % get path to save timing analysis

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

distcomp.feature( 'LocalUseMpiexec', false ); % added to help the pool work...
num_workers = 20;
parallel_pool = parpool(num_workers);

tic
x = linspace(250,2250,5);
niters = 10000;
est_params = nan(niters,3); % iteration by params (slope, int, elbow point)
rng(1)
parfor n = 1:niters      
    p = randi(n_subs,n_subs,1);
    y = mean(azimuth_errors(:,p),2)';
    est_params_avg = fit_elbow_func_localization(x,y);  
    est_params(n,:,:) = est_params_avg;
end
toc
save(timing_analysis_file, 'est_params')
delete(parallel_pool)