[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
data_path = fullfile(file_path, sprintf('../../results/%s/raw/', exp_name)); % get path to raw data
save_path = fullfile(file_path, sprintf('../../results/%s/', exp_name)); % get path to save results

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

files = dir(fullfile(data_path, 'sub-*.csv')); % get list of files
n_trials = 320; % exp has 320 trials
repcode = repmat(1:2, 1, n_trials/2)'; % 1==repeated bg and 2==non-repeated bg
threshold = 0.8; % performance threshold for average d' (increased threshold since experiment was easier)
results = struct; % create an empty struct to store data in
k = 0; % counter for number of subjects with usable data

for f = 1:length(files)
    file_name = fullfile(data_path, files(f).name);   
    data_table = readtable(file_name,'Delimiter',',', 'HeaderLines',0); % reads in the batch csv from prolific
    if size(data_table,1)~=1
        data_table = data_table(1,:);
    end
    if ~ismember(data_table.Properties.VariableNames, 'BonusAmount') % if they don't have any bonus, then they failed the headphone check and didn't do the experiment OR didn't finish it
    elseif strcmp(data_table.Q6hearingloss{:}, 'yesHL') % do not include if they self-report hearing loss
    else        
        responses = zeros(n_trials,6);  % preallocate an array
        for trial = 1:n_trials % loop through all the trials                
            stimID_var = sprintf('stimID_for_expt_1trial_%d',trial); % store the variable name to index the table for the given trial                
            resp_var = sprintf('Resp_for_expt_1trial_%d',trial); % store the variable name to index the table for the given trial
            has_data = sum(ismember(data_table.Properties.VariableNames, stimID_var)) + sum(ismember(data_table.Properties.VariableNames, resp_var)) == 2;
            if has_data
                responses(trial,1:5) = str2num(strrep(data_table.(stimID_var){1},'_',' ')); % access the stimID from the table and convert into an array
                responses(trial,6) = data_table.(resp_var); % get the subject response                    
            end
        end
        if sum(sum(responses, 2)==0) <= round(0.1 * n_trials) % subject is missing at most 10% of trials
            k = k + 1; % increase counter since subject has usable data
            results(k).responses = responses; % store stimulus info and subject responses organized along columns as: | bg | bg exemplar | fg |  SNR | fg onset time | fg on/off | subject response |             
            results(k).demographics = {data_table.Q1Gender{:}, data_table.Q2age, data_table.Q7music{:}}; % store demographic info

            responses_to_use = responses(:, [1 3:6]); % get relevant stimulus/condition info
            responses_to_use(:, 2) = repcode;
            responses_to_use(sum(responses,2)==0,:) = []; % remove trials with no data                                        

            [G, ID1, ID2] = findgroups(responses_to_use(:,2),responses_to_use(:,3)); % get the group info
            group_info = [ID1 ID2];
            hit_rates = splitapply(@compute_hit_rate, responses_to_use(:, 4), responses_to_use(:, 5), G);  % get hit rates and false alarm rates for all conditions (SNR by onset time)      
            mean_false_alarm_rates = splitapply(@compute_false_alarm_rate, responses_to_use(:, 4), responses_to_use(:, 5), responses_to_use(:, 2)); % compute false alarms over all "conditions" since they only apply with foregrounds present
            temp_false_alarm_rates = repmat(mean_false_alarm_rates, 1, 10)';

            results(k).hit_rates = hit_rates; % store hit rates
            results(k).mean_false_alarm_rates = mean_false_alarm_rates; % store mean false alarm rate
            results(k).d_primes = norminv(hit_rates) - norminv(temp_false_alarm_rates(:)); % store d'
            results(k).overall_performance = mean(results(k).d_primes); % store average d'
            if results(k).overall_performance >= threshold
                results(k).pass = true;
            else
                results(k).pass = false;
            end
            results(k).group_info = group_info; % store condition/grouping info (helpful for computing stats)
        end
    end
end

save_name = fullfile(save_path, 'results.mat');
save(save_name, 'results')