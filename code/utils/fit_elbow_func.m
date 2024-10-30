function parameters = fit_elbow_func(x,y)
    % this function is used for fitting elbow functions to dprime data
    % if used for something else, the initial params estimate (params0) might be off a bit 
    % (i.e. the random part should be scaled for the problem at hand)
    
    error_function = @(params) sum(abs(y - elbow_func(x,params)));
    n_seeds = 100;
    sigmas = [0.01, 0.5, 500];
    options = optimset('Display','off');
    all_params = nan(n_seeds,3);
    all_fvals = nan(n_seeds,1);
    b = regress(y', [x', ones(size(x'))]);
    for n = 1:n_seeds
        params0 = [randn(1)*sigmas(1) + b(1), randn(1)*sigmas(2) + b(2), randn(1)*sigmas(3) + mean(x)];
        [est_params,fval] = fmincon(error_function,params0,[],[],[],[],[0, min(y) - range(y)/min(diff(x)) * x(1), x(1)],[2*range(y)/min(diff(x)), max(y), x(end)],[],options);
        all_params(n,:) = est_params;
        all_fvals(n) = fval;
    end
    [~,index] = min(all_fvals);
    parameters = all_params(index,:);
end