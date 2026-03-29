%% sygklish

data = load("convergence.txt");
data2 = load('convergence_b.txt');
data3 = load('convergence_c.txt');
data4 = load('convergence_d.txt');

iter = data(:,1);
iterb = data2(:,1);
iterc = data3(:,1);
iterd = data4(:,1);

max_diff = data(:,2);
max_diffb = data2(:,2);
max_diffc = data3(:,2);
max_diffd = data4(:,2);


figure;
hold on; grid on; box on;

% FVS
plot(iter,  max_diff,  '--r', 'LineWidth', 1.2);   % FVS 1ης τάξης
plot(iterb, max_diffb, '-r',  'LineWidth', 1.2);   % FVS 2ης τάξης

% Roe
plot(iterc, max_diffc, '--',  'Color', [0 0 0.6], 'LineWidth', 1.2); % Roe 1ης τάξης
plot(iterd, max_diffd, '-',   'Color', [0 0 0.6], 'LineWidth', 1.2); % Roe 2ης τάξης

set(gca, 'YScale', 'log');

xlabel('Αριθμός επαναλήψεων');
ylabel('Residual (max |ΔU|)');
title('Σύγκλιση αριθμητικών μεθόδων');

legend( ...
    'FVS 1ης τάξης', ...
    'FVS 2ης τάξης', ...
    'Roe 1ης τάξης', ...
    'Roe 2ης τάξης', ...
    'Location','northeast' ...
);

set(gca,'FontSize',12);
title('Σύγκλιση ΜΕΘΟΔΩΝ');
