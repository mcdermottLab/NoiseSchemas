function hit_rate = compute_hit_rate(correct_answers, responses)
    % correct_answers: 1 => signal absent, 2 => signal present
    on_trials = correct_answers == 2;
    n_trials = sum(on_trials);
    hit_rate = mean(responses(on_trials)-1);
    hit_rate = correct_extreme_rates(hit_rate, n_trials);
end