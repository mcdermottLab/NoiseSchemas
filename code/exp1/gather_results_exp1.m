[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
data_path = fullfile(file_path, sprintf('../../results/%s/raw/', exp_name)); % get path to raw data
save_path = fullfile(file_path, sprintf('../../results/%s/', exp_name)); % get path to save results

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

files = dir(fullfile(data_path, 'Batch*.csv')); % get list of files
n_trials = 330; % exp1 has 320 main trials plus 10 catch trials
threshold = 0.6; % performance threshold for average d'
results = struct; % create an empty struct to store data in
k = 0; % counter for number of subjects with usable data
for f = 1:length(files)
    batch_file_name = fullfile(data_path, files(f).name);   
    batch_info = readtable(batch_file_name,'Delimiter',',', 'HeaderLines',0); % reads in the batch csv from mturk
    for sub = 1:size(batch_info,1) % loop through all subjects
        if isnan(batch_info.Answer_BonusAmount(sub)) % if they don't have any bonus, then they failed the headphone check and didn't do the experiment OR didn't finish it
        elseif strcmp(batch_info.Answer_Q6hearingloss{sub}, 'yesHL') % do not include if they self-report hearing loss
        else
            k = k + 1; % increase counter since subject has usable data
            responses = zeros(n_trials,5);  % preallocate an array
            for trial = 1:n_trials % loop through all the trials
                    stimID_var = sprintf('Answer_stimID_for_expt_1trial_%d',trial); % store the variable name to index the table for the given trial
                    responses(trial,1:4) = str2num(strrep(batch_info{sub,stimID_var}{:}, '_', ' ')); % extract 4 stimulus parameters (see below)
                    resp_var = sprintf('Answer_Resp_for_expt_1trial_%d',trial); % store the variable name to index the table for the given trial                    
                    responses(trial,5) = batch_info{sub,resp_var}; % extract subject response
            end
            results(k).responses = responses; % store stimulus info and subject responses organized along columns as: | bg/fg paring | SNR | fg onset time | fg on/off | subject response |             
            results(k).demographics = {batch_info.Answer_Q1Gender{sub}, batch_info.Answer_Q2age(sub), batch_info.Answer_Q7music{sub}}; % store demographic info

            catch_trials = responses(:,1) > 900; % get list of catch trials            
            non_catch_trials = ~catch_trials; % get list of non-catch (i.e., regular experimental) trials
            
            [G, ID1, ID2] = findgroups(responses(non_catch_trials,2),responses(non_catch_trials,3)); % get the group info
            group_info = [ID1 ID2];
            hit_rates = splitapply(@compute_hit_rate, responses(non_catch_trials, 4), responses(non_catch_trials, 5), G);  % get hit rates and false alarm rates for all conditions (SNR by onset time)      
            mean_false_alarm_rate = compute_false_alarm_rate(responses(non_catch_trials, 4), responses(non_catch_trials, 5)); % compute false alarms over all "conditions" since they only apply with foregrounds present

            results(k).hit_rates = hit_rates; % store hit rates
            results(k).mean_false_alarm_rate = mean_false_alarm_rate; % store mean false alarm rate
            results(k).d_primes = norminv(hit_rates) - norminv(mean_false_alarm_rate); % store d'
            results(k).overall_performance = mean(results(k).d_primes); % store average d'
            if results(k).overall_performance >= threshold
                results(k).pass = true;
            else
                results(k).pass = false;
            end
            results(k).group_info = group_info; % store grouping info
        end
    end 
end

save_name = fullfile(save_path, 'results.mat');
save(save_name, 'results')