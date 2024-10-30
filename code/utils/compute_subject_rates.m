function [hit_rates, false_alarm_rates, group_info] = compute_subject_rates(responses)
    [G, ID1, ID2] = findgroups(responses(:,2),responses(:,3)); % get the group info
    [hit_rates, false_alarm_rates] = splitapply(@compute_rates, responses, G);
    group_info = [ID1, ID2];
end

function [hit_rate, false_alarm_rate] = compute_rates(data)
    on_trials = data(:,4) == 2;
    hit_rate = mean(data(on_trials,5)-1);

    off_trials = ~on_trials;
    false_alarm_rate = mean(data(off_trials,5)-1);
end