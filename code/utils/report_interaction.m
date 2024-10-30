function report_interaction(rm, C, ranovatbl, var1, var2)
    m_results = mauchly(rm, C);
    if (m_results.DF ~= 0) && (m_results.pValue < 0.05)
        fprintf('Sphericity assumption violated! Using GG correction.\n')
        pvalname = 'pValueGG';
        correction_name = 'GreenhouseGeisser';
    else
        fprintf('Sphericity assumption satisfied.\n')
        pvalname = 'pValue';
        correction_name = 'Uncorrected';
    end
      
    epsilon_table = epsilon(rm, C);
    eps = epsilon_table.(correction_name);
    row_index = find(strcmp(ranovatbl.Properties.RowNames, sprintf('(Intercept):%s:%s', var1, var2)));
    F = ranovatbl(row_index,:).F;
    df1 = ranovatbl(row_index,:).DF * eps;
    df2 = ranovatbl(row_index+1,:).DF * eps;
    pvalue = ranovatbl(row_index,:).(pvalname);
    if pvalue < 0.001
        pvalue_str = 'p<0.001';
    else
        pvalue_str = sprintf('p=%.2f', pvalue);
    end
    partial_eta_squared = ranovatbl(row_index,:).SumSq / (ranovatbl(row_index,:).SumSq + ranovatbl(row_index+1,:).SumSq);
    
    fprintf('Interaction Between %s & %s\n', var1, var2)
    fprintf('F(%.2f,%.2f)=%.2f, %s\n', df1, df2, F, pvalue_str)
    fprintf('Partial Eta Squared: %.2f\n\n', partial_eta_squared)

%     fprintf('DEBUGGING:\n')
%     fprintf('Table p-value: %.6f\n', pvalue)
%     pval_calc = fcdf(F, df1, df2, 'upper');
%     fprintf('Calc. p-value: %.6f\n\n', pval_calc)
end