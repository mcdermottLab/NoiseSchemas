[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
comparison_exp = 'exp6';
permutation_test_file = fullfile(file_path, sprintf('../../results/%s/matched_permutation_test_%s_vs_%s.mat', exp_name, exp_name, comparison_exp)); % get path to save permutation test results

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

d_primes_nonrep = d_primes;
n_nonrep = n_subs;

prev_exp_name = sprintf('exp%d', str2num(exp_name(4:end))-1);
prev_results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', prev_exp_name)); % get path to results
load(prev_results_file)

d_primes_rep = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes_rep(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes_rep = d_primes_rep(:, logical(keep_subs));
n_rep = size(d_primes_rep,2);
 
% MATCHING PROCEDURE
% match participants across experiments
indices = [6, 8, 10]; % corresponds to timepoints 1500, 2000, 2500 ms
test_indices = setdiff(1:10, indices);
max_d = 0.025; % maximum acceptable difference for subjects
offsets = -0.2:0.01:0.2; % range of offsets to use
deltas = nan(length(offsets), length(indices));
for ind = 1:length(indices) % length(indices)-fold cross-val to select offset
    inds_to_use = setdiff(1:length(indices), ind);
    for k = 1:length(offsets)   
        mean_scores = mean(d_primes_rep(indices(inds_to_use),:),'omitnan');
        mean_scores_comparison = mean(d_primes_nonrep(indices(inds_to_use),:),'omitnan');
        offset = offsets(k);
        matched_subs = [];
        matched_subs_comparison = [];
        [~, rp] = sort(abs(mean_scores - mean([mean_scores, mean_scores_comparison]))); % sort participants by those closest to the average performance across experiments
        for i = 1:length(mean_scores)
            index = rp(i);
            [min_val, comparison_index] = min(abs(mean_scores_comparison - mean_scores(index) + offset)); % for ith participant, find a matching subject from comparison experiment
            if min_val < max_d
                matched_subs(end+1) = index;
                matched_subs_comparison(end+1) = comparison_index;
                mean_scores_comparison(comparison_index) = nan;
            end
        end    
        % compute performance difference (delta) across the matched groups in validation data
        delta = mean(d_primes_nonrep(indices(ind),matched_subs_comparison),'all','omitnan') - mean(d_primes_rep(indices(ind),matched_subs),'all','omitnan');    
        deltas(k, ind) = delta;
    end
end

[~,k] = min(mean(abs(deltas),2)); % find which offset results in the smallest across group difference (delta)
% re-do matching with the optimal offset
mean_scores = mean(d_primes_rep(indices,:),'omitnan');
mean_scores_comparison = mean(d_primes_nonrep(indices,:),'omitnan');
offset = offsets(k);
matched_subs = [];
matched_subs_comparison = [];
[~, rp] = sort(abs(mean_scores - mean([mean_scores, mean_scores_comparison])));
for i = 1:length(mean_scores)
    index = rp(i);
    [min_val, comparison_index] = min(abs(mean_scores_comparison - mean_scores(index) + offset));
    if min_val < max_d
        matched_subs(end+1) = index;
        matched_subs_comparison(end+1) = comparison_index;
        mean_scores_comparison(comparison_index) = nan;
    end
end

d_primes_matched = d_primes_rep(:, matched_subs);
n_subs = size(d_primes_matched,2);

d_primes_matched_comparison = d_primes_nonrep(:, matched_subs_comparison);
n_subs_comparison = size(d_primes_matched_comparison,2);

% fit elbow function to actual data
x = linspace(250,2500,10);
x = x(test_indices);
x_test = 250:2500;
y = mean(d_primes_matched(test_indices,:),2,'omitnan')';
est_params = fit_elbow_func(x,y);
y_pred = elbow_func(x_test, est_params);

y_comparison = mean(d_primes_matched_comparison(test_indices,:),2,'omitnan')';
est_params_comparison = fit_elbow_func(x,y_comparison);
y_pred_comparison = elbow_func(x_test, est_params_comparison);

magnitude_statistic = (y_pred(end)-y_pred(1)) - (y_pred_comparison(end)-y_pred_comparison(1));

all_d_primes_matched = [d_primes_matched, d_primes_matched_comparison];

niters = 10000;
rng(1)
magnitude_nullstats = nan(niters, 1);

distcomp.feature( 'LocalUseMpiexec', false ); % added to help the pool work...
num_workers = 20;
parallel_pool = parpool(num_workers);

tic
parfor n = 1:niters
    rp = randperm(n_subs + n_subs_comparison);
    perm_dp = all_d_primes_matched(test_indices, rp(1:n_subs));
    perm_y = mean(perm_dp,2,'omitnan')';
    perm_est_params = fit_elbow_func(x,perm_y);
    perm_y_pred = elbow_func(x_test, perm_est_params);

    perm_dp_comparison = all_d_primes_matched(test_indices, rp((1+n_subs):end));
    perm_y_comparison = mean(perm_dp_comparison,2,'omitnan')';
    perm_est_params_comparison = fit_elbow_func(x,perm_y_comparison);
    perm_y_pred_comparison = elbow_func(x_test, perm_est_params_comparison);

    magnitude_nullstats(n) = (perm_y_pred(end)-perm_y_pred(1)) - (perm_y_pred_comparison(end)-perm_y_pred_comparison(1));
end
toc

save(permutation_test_file, 'magnitude_nullstats', 'magnitude_statistic', 'matched_subs', 'matched_subs_comparison', 'offset')
delete(parallel_pool)