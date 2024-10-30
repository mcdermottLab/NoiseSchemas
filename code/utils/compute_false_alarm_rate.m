function false_alarm_rate = compute_false_alarm_rate(correct_answers, responses)
    % correct_answers: 1 => signal absent, 2 => signal present
    off_trials = correct_answers == 1;
    n_trials = sum(off_trials);
    false_alarm_rate = mean(responses(off_trials)-1);
    false_alarm_rate = correct_extreme_rates(false_alarm_rate, n_trials);
end