function y = elbow_func(x,params)
    y = NaN(size(x));    
    for i = 1:length(x)
        if x(i) < params(3)
            y(i) = params(1)*x(i) + params(2);
        else
            y(i) = params(1)*params(3) + params(2);
        end
    end
end