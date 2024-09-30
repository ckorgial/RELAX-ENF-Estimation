% List of scheme names and MSE values
schemes = {'single', 'E-single', 'MLE', 'WMLE', 'E-MLE', 'E-WMLE', 'S-MLE', 'S-WMLE', 'P-MLE', 'P-WMLE', ...
           'LAD-single', 'LAD-E-single', 'LAD-MLE', 'LAD-WMLE','LAD-E-MLE', 'LAD-E-WMLE', ...
           'LAD-S-MLE', 'LAD-S-WMLE', 'LAD-P-MLE', 'LAD-P-WMLE', 'RELAX-single', 'RELAX-E-single', 'RELAX-MLE', ...
           'RELAX-WMLE', 'RELAX-E-MLE', 'RELAX-E-WMLE', 'RELAX-S-MLE', 'RELAX-S-WMLE', 'RELAX-P-MLE', 'RELAX-P-WMLE'};
MSE_values = {MSE_single, MSE_E_single, MSE_MLE, MSE_WMLE, MSE_E_MLE, MSE_E_WMLE, MSE_S_MLE, MSE_S_WMLE, ...
              MSE_P_MLE, MSE_P_WMLE, MSE_LAD_single, MSE_LAD_E_single, MSE_LAD_E_MLE, MSE_LAD_E_WMLE, ...
              MSE_LAD_MLE, MSE_LAD_WMLE, MSE_LAD_S_MLE, MSE_LAD_S_WMLE, MSE_LAD_P_MLE, MSE_LAD_P_WMLE, ...
              MSE_single_relax, MSE_E_single_relax, MSE_MLE_relax, MSE_WMLE_relax, MSE_E_MLE_relax, MSE_E_WMLE_relax, ...
              MSE_S_MLE_relax, MSE_S_WMLE_relax, MSE_P_MLE_relax, MSE_P_WMLE_relax};

% Extract the RELAX schemes and non-RELAX schemes
relax_indices = find(contains(schemes, 'RELAX'));
non_relax_indices = find(~contains(schemes, 'RELAX'));

% Initialize matrix to store p-values
p_values = zeros(length(relax_indices), length(non_relax_indices));

% Perform paired t-tests for RELAX schemes against all non-RELAX schemes
for i = 1:length(relax_indices)
    for j = 1:length(non_relax_indices)
        [~, p_value] = ttest(MSE_values{relax_indices(i)}, MSE_values{non_relax_indices(j)});
        p_values(i, j) = p_value;
    end
end

% Extract the scheme names for RELAX and non-RELAX schemes
relax_schemes = schemes(relax_indices);
non_relax_schemes = schemes(non_relax_indices);

% Create a table to display p-values
p_values_table = array2table(p_values, 'VariableNames', non_relax_schemes, 'RowNames', relax_schemes);

% Display the table of p-values
disp('Paired t-test p-values (RELAX vs. non-RELAX schemes):');
disp(p_values_table);

% Save the table of p-values to a CSV file
writetable(p_values_table, 'p_values_table_relax_vs_non_relax.csv', 'WriteRowNames', true);

% Determine which comparisons are statistically significant
significance_threshold = 0.05;
significant_differences = p_values < significance_threshold;

% Create a table to display significant differences
significant_table = array2table(significant_differences, 'VariableNames', non_relax_schemes, 'RowNames', relax_schemes);

% Display the table of significant differences
disp('Significant differences (p < 0.05) (RELAX vs. non-RELAX schemes):');
disp(significant_table);

% Save the table of significant differences to a CSV file
writetable(significant_table, 'significant_differences_table_relax_vs_non_relax.csv', 'WriteRowNames', true);
