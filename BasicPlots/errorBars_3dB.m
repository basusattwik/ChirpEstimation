% Data from the LaTeX table
means_CGLMC = [10.17, 40.85, -68.96, 109.84, 50.08, 62.58, -104.46, 112.26];
stds_CGLMC = [0.55, 1.88, 3.91, 2.40, 0.55, 1.27, 7.05, 3.91];
abs_error_CGLMC = [0.17, 0.85, 1.04, 0.16, 0.08, 2.58, 14.46, 7.26];

means_LMC = [10.05, 39.98, -80.43, 134.39, 50.55, 60.33, -102.16, 129.99];
stds_LMC = [0.97, 8.71, 24.51, 21.57, 0.92, 4.74, 7.19, 17.15];
abs_error_LMC = [0.05, 0.02, 10.43, 24.39, 0.55, 0.33, 12.16, 24.99];

means_NALMC = [9.93, 41.20, -77.04, 126.36, 50.27, 57.34, -85.59, 106.80];
stds_NALMC = [1.58, 10.93, 15.63, 24.37, 1.20, 9.02, 18.26, 8.20];
abs_error_NALMC = [0.07, 1.20, 7.04, 16.36, 0.27, 2.66, 4.41, 1.80];

% Organize data for boxplot
data = [
    means_CGLMC', stds_CGLMC', abs_error_CGLMC';  % CG-LMC
    means_LMC', stds_LMC', abs_error_LMC';        % LMC
    means_NALMC', stds_NALMC', abs_error_NALMC'   % NA-LMC
];

% Grouping variables
group = [repmat({'CG-LMC'}, 8, 1); repmat({'LMC'}, 8, 1); repmat({'NA-LMC'}, 8, 1)];
category_labels = [repmat({'Mean'}, 8, 1); repmat({'Std'}, 8, 1); repmat({'Abs Error'}, 8, 1)];

% Reshape data for plotting
reshaped_data = reshape(data, [], 1);
method_labels = [repmat({'CG-LMC'}, 8*3, 1); repmat({'LMC'}, 8*3, 1); repmat({'NA-LMC'}, 8*3, 1)];
data_labels = [repmat({'Mean'}, 8, 1); repmat({'Std'}, 8, 1); repmat({'Abs Error'}, 8, 1)];

% Creating the boxplot
figure;
boxplot(reshaped_data, {method_labels, data_labels}, 'factorgap', 5, 'colorgroup', method_labels, 'colors', 'brg', 'Widths', 0.5);

xlabel('Method and Data Type');
ylabel('Values');
title('Box Plot for CG-LMC, LMC, and NA-LMC (Mean, Std, Abs Error)');
grid on;

% Add a legend for the different groups
legend('CG-LMC', 'LMC', 'NA-LMC', 'Location', 'BestOutside');
