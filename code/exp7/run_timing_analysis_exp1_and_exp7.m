[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name

timing_analysis_file = fullfile(file_path, sprintf('../../results/%s/timing_analysis_exp1_and_exp7.mat', exp_name)); % get path to save timing analysis

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

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

d_primes = [d_primes_comparison1, d_primes_comparison7];
n_subs = size(d_primes, 2);

distcomp.feature( 'LocalUseMpiexec', false ); % added to help the pool work...
num_workers = 20;
parallel_pool = parpool(num_workers);

tic
x = linspace(250,2500,10);
niters = 10000;
est_params = nan(niters,3); % iteration by params (slope, int, elbow point)
rng(1)
parfor n = 1:niters   
    p = randi(n_subs,n_subs,1);
    y = mean(d_primes(:,p),2)';
    est_params_avg = fit_elbow_func(x,y);  
    est_params(n,:) = est_params_avg;
end
toc
save(timing_analysis_file, 'est_params')
delete(parallel_pool)