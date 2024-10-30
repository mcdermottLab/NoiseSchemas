function corrected_rates = correct_extreme_rates(rates, n_trials)
    % apply correction (suggested by Stanislaw & Todorov 1999) to 
    % extreme hit / false alarm rates (1 or 0) in order to avoid infinities
    % in d' calculations    
    corrected_rates = rates;
    corrected_rates(corrected_rates==1) = (n_trials - 0.5) / n_trials;
    corrected_rates(corrected_rates==0) = 0.5 / n_trials;
end