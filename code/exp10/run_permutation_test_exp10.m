[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
permutation_test_file = fullfile(file_path, sprintf('../../results/%s/permutation_test.mat', exp_name)); % get path to save permutation test results

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
n_subs = size(d_primes,2);

% load comparison data from exp 1
comparison_exp = 'exp1';
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', comparison_exp)); % get path to results
load(results_file)
d_primes_comparison = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes_comparison(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end
d_primes_comparison1 = d_primes_comparison(:, logical(keep_subs));
d_primes_comparison1 = (d_primes_comparison1(1:10, :) + d_primes_comparison1(11:20, :)) / 2;

% load comparison data from exp 7
comparison_exp = 'exp7';
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', comparison_exp)); % get path to results
load(results_file)
d_primes_comparison = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes_comparison(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end
d_primes_comparison7 = d_primes_comparison(:, logical(keep_subs));

d_primes_comparison = [d_primes_comparison1, d_primes_comparison7];
n_subs_comparison = size(d_primes_comparison, 2);

% fit elbow function to actual data
x = linspace(250,2500,10);
x_test = 250:2500;

y = mean(d_primes,2)';
est_params = fit_elbow_func(x,y);
y_pred = elbow_func(x_test, est_params);

y_comparison = mean(d_primes_comparison,2)';
est_params_comparison = fit_elbow_func(x,y_comparison);
y_pred_comparison = elbow_func(x_test, est_params_comparison);

timescale_statistic = est_params(3) - est_params_comparison(3);
magnitude_statistic = (y_pred(end)-y_pred(1)) - (y_pred_comparison(end)-y_pred_comparison(1));

all_d_primes = [d_primes, d_primes_comparison];
rng(1)
niters = 10000;
timescale_nullstats = nan(niters, 1);
magnitude_nullstats = nan(niters, 1);

distcomp.feature( 'LocalUseMpiexec', false ); % added to help the pool work...
num_workers = 20;
parallel_pool = parpool(num_workers);

tic
parfor n = 1:niters
    rp = randperm(n_subs + n_subs_comparison);
    perm_dp = all_d_primes(:, rp(1:n_subs));
    perm_y = mean(perm_dp,2)';
    perm_est_params = fit_elbow_func(x,perm_y);
    perm_y_pred = elbow_func(x_test, perm_est_params);

    perm_dp_comparison = all_d_primes(:, rp((1+n_subs):end));
    perm_y_comparison = mean(perm_dp_comparison,2)';
    perm_est_params_comparison = fit_elbow_func(x,perm_y_comparison);
    perm_y_pred_comparison = elbow_func(x_test, perm_est_params_comparison);

    timescale_nullstats(n) = perm_est_params(3) - perm_est_params_comparison(3);
    magnitude_nullstats(n) = (perm_y_pred(end)-perm_y_pred(1)) - (perm_y_pred_comparison(end)-perm_y_pred_comparison(1));
end
toc

save(permutation_test_file, 'timescale_nullstats', 'timescale_statistic', 'magnitude_nullstats', 'magnitude_statistic')
delete(parallel_pool)