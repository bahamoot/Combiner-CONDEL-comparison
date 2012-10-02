function compare_roc_result_with_condel(working_dir)

load(strcat(working_dir, '/data/combiner_params_0.2897.mat'));

pathogenic_ratio     = 0.5;

condel_result = importdata(strcat(working_dir, '/bin/tmp/condel_score'), '\t');
features      = condel_result.data(:, 1:7)';
targets       = condel_result.data(:, 8)';

neutral_SNP_index    = find(targets == 0);
neutral_SNP_features = features(:, neutral_SNP_index);
neutral_SNP_targets  = targets(:, neutral_SNP_index);
[number_of_features number_of_neutral_samples] = size(neutral_SNP_features);


[neutral_out1, neutral_out2, neutral_out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, neutral_SNP_features(1:6, :));
combiner_neutral_error = (sum(sum(abs(neutral_out3 - neutral_SNP_targets) ./ (number_of_neutral_samples .* ones(size(neutral_SNP_targets))))));
condel_neutral_error = (sum(sum(abs(neutral_SNP_features(7, :) - neutral_SNP_targets) ./ (number_of_neutral_samples .* ones(size(neutral_SNP_targets))))));

combiner_false_positive_group = find(neutral_out3 > pathogenic_ratio);
combiner_true_negative_group  = find(neutral_out3 <= pathogenic_ratio);
condel_false_positive_group   = find(neutral_SNP_features(7, :) > pathogenic_ratio);
condel_true_negative_group    = find(neutral_SNP_features(7, :) <= pathogenic_ratio);

combiner_number_of_false_positive = size(combiner_false_positive_group, 2);
combiner_number_of_true_negative  = size(combiner_true_negative_group, 2);
condel_number_of_false_positive   = size(condel_false_positive_group, 2);
condel_number_of_true_negative    = size(condel_true_negative_group, 2);


pathogenic_SNP_index    = find(targets == 1);
pathogenic_SNP_features = features(:, pathogenic_SNP_index);
pathogenic_SNP_targets  = targets(:, pathogenic_SNP_index);
number_of_pathogenic_samples = size(pathogenic_SNP_features, 2);

[pathogenic_out1, pathogenic_out2, pathogenic_out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, pathogenic_SNP_features(1:6, :));
combiner_pathogenic_error = (sum(sum(abs(pathogenic_out3 - pathogenic_SNP_targets) ./ (number_of_pathogenic_samples .* ones(size(pathogenic_SNP_targets))))));
condel_pathogenic_error = (sum(sum(abs(pathogenic_SNP_features(7, :) - pathogenic_SNP_targets) ./ (number_of_pathogenic_samples .* ones(size(pathogenic_SNP_targets))))));

combiner_false_negative_group = find(pathogenic_out3 <= pathogenic_ratio);
combiner_true_positive_group  = find(pathogenic_out3 > pathogenic_ratio);
condel_false_negative_group   = find(pathogenic_SNP_features(7, :) <= pathogenic_ratio);
condel_true_positive_group    = find(pathogenic_SNP_features(7, :) > pathogenic_ratio);

combiner_number_of_false_negative = size(combiner_false_negative_group, 2);
combiner_number_of_true_positive  = size(combiner_true_positive_group, 2);
condel_number_of_false_negative   = size(condel_false_negative_group, 2);
condel_number_of_true_positive    = size(condel_true_positive_group, 2);


roc_range = 0.0:0.001:1;
FP_rates = [];
TP_rates = [];

for roc_pathogenic_ratio = roc_range
    FP_rates = [FP_rates, sum([neutral_SNP_features(2:number_of_features-1,:); neutral_out3] > roc_pathogenic_ratio, 2) / number_of_neutral_samples];
    TP_rates = [TP_rates, sum([pathogenic_SNP_features(2:number_of_features-1,:); pathogenic_out3] > roc_pathogenic_ratio, 2) / number_of_pathogenic_samples];
end;

min_gerp = min(features(1,:));
max_gerp = max(features(1,:));

gerp_roc_range = min_gerp:(max_gerp-min_gerp)/1000:max_gerp;
FP_rate =[];
TP_rate = [];

for gerp_pathogenic_ratio = gerp_roc_range
    FP_rate = [FP_rate, sum([neutral_SNP_features(1,:)] > gerp_pathogenic_ratio, 2) / number_of_neutral_samples];
    TP_rate = [TP_rate, sum([pathogenic_SNP_features(1,:)] > gerp_pathogenic_ratio, 2) / number_of_pathogenic_samples];
end;

FP_rates = [FP_rate; FP_rates];
TP_rates = [TP_rate; TP_rates];

min_Condel = min(features(7,:));
max_Condel = max(features(7,:));    

Condel_roc_range = min_Condel:(max_Condel-min_Condel)/1000:max_Condel;
FP_rate =[];
TP_rate = [];

for Condel_pathogenic_ratio = Condel_roc_range
    FP_rate = [FP_rate, sum([neutral_SNP_features(7,:)] > Condel_pathogenic_ratio, 2) / number_of_neutral_samples];
    TP_rate = [TP_rate, sum([pathogenic_SNP_features(7,:)] > Condel_pathogenic_ratio, 2) / number_of_pathogenic_samples];
end;

FP_rates = [FP_rates; FP_rate];
TP_rates = [TP_rates; TP_rate];


f = figure;
hold on
plot(FP_rates(1,:), TP_rates(1,:), '--m');
plot(FP_rates(2,:), TP_rates(2,:), '--b');
plot(FP_rates(3,:), TP_rates(3,:), '--g');
plot(FP_rates(4,:), TP_rates(4,:), '--k');
plot(FP_rates(5,:), TP_rates(5,:), '--c');
plot(FP_rates(6,:), TP_rates(6,:), '--r');
plot(FP_rates(7,:), TP_rates(7,:), 'k');
plot(FP_rates(8,:), TP_rates(8,:), 'r');
hold off

ylabel('True positive rate', 'FontName', 'Arial', 'FontSize', 18);
xlabel('False positive rate', 'FontName', 'Arial', 'FontSize', 18);
title('ROC curve', 'FontName', 'Arial', 'FontSize', 18);
h = legend('GERP++', 'SIFT', 'polyphen2', 'phylop', 'LRT', 'Mutation Taster', 'Combiner', 'Condel');
set(h, 'FontName','Arial');
set(h, 'FontSize',18);
set(h, 'Position', [0.5 0.35 0.45 0.15]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6]);
roc_filename = strcat(working_dir, '/result/roc.bmp');
print(f, '-dbmp16m', roc_filename);

subplot(3,2,3);
neutral_distribution = hist(neutral_SNP_features(3,:), [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_SNP_features(3,:), [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
hold off
ylabel('samples', 'FontName', 'Arial', 'FontSize', 18);
xlabel('score', 'FontName', 'Arial', 'FontSize', 18);
title('Polyphen2 score distribution', 'FontName', 'Arial', 'FontSize', 18);

subplot(3,2,4);
neutral_distribution = hist(neutral_SNP_features(6,:), [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_SNP_features(6,:), [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
hold off
ylabel('samples', 'FontName', 'Arial', 'FontSize', 18);
xlabel('score', 'FontName', 'Arial', 'FontSize', 18);
title('Mutation Taster score distribution', 'FontName', 'Arial', 'FontSize', 18);

Condel_distribution_range = min_Condel:(max_Condel-min_Condel)/50:max_Condel;

subplot(3,2,5);
neutral_distribution = hist(neutral_SNP_features(7,:), Condel_distribution_range);
plot(Condel_distribution_range, neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_SNP_features(7,:), Condel_distribution_range);
plot(Condel_distribution_range, pathogenic_distribution, '-.r');
hold off
ylabel('samples', 'FontName', 'Arial', 'FontSize', 18);
xlabel('score', 'FontName', 'Arial', 'FontSize', 18);
title('Condel score distribution', 'FontName', 'Arial', 'FontSize', 18);

subplot(3,2,6);
neutral_distribution = hist(neutral_out3, [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_out3, [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
hold off
ylabel('samples');
xlabel('score');
ylabel('samples', 'FontName', 'Arial', 'FontSize', 18);
xlabel('score', 'FontName', 'Arial', 'FontSize', 18);
title('Combiner score distribution', 'FontName', 'Arial', 'FontSize', 18);

h = legend('neutral variants', 'pathogenic variants');
set(h, 'FontName','Arial');
set(h, 'FontSize',18);
set(h, 'Position', [0.265 0.65 0.5 0.1]);

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12]);
score_distribution_filename = strcat(working_dir, '/result/score_distribution_combiner.bmp');
print(f, '-dbmp16m', score_distribution_filename);

close(f);

size_neutral_training_data    = sum(training_targets == 0);
size_pathogenic_training_data = sum(training_targets == 1);
size_training_data = size_neutral_training_data + size_pathogenic_training_data;
size_neutral_validation_data    = sum(validating_targets == 0);
size_pathogenic_validation_data = sum(validating_targets == 1);
size_validation_data = size_neutral_validation_data + size_pathogenic_validation_data;
size_neutral_test_data    = sum(testing_targets == 0);
size_pathogenic_test_data = sum(testing_targets == 1);
size_test_data = size_neutral_test_data + size_pathogenic_test_data;
size_neutral_total_data    = size_neutral_training_data + size_neutral_validation_data + size_neutral_test_data;
size_pathogenic_total_data = size_pathogenic_training_data + size_pathogenic_validation_data + size_pathogenic_test_data;
size_total_data = size_neutral_total_data + size_pathogenic_total_data;


disp(' ');
disp(sprintf('%32s%-15s%-15s%-15s', '', 'Combiner', 'Condel', 'Total'));
disp(sprintf('%-30s: %-15s%-15s%-15s', 'Training dataset', sprintf('%d', size_neutral_training_data), sprintf('%d', size_pathogenic_training_data), sprintf('%d', size_training_data)));
disp(sprintf('%-30s: %-15s%-15s%-15s', 'Validation dataset', sprintf('%d', size_neutral_validation_data), sprintf('%d', size_pathogenic_validation_data), sprintf('%d', size_validation_data)));
disp(sprintf('%-30s: %-15s%-15s%-15s', 'Test dataset', sprintf('%d', size_neutral_test_data), sprintf('%d', size_pathogenic_test_data), sprintf('%d', size_test_data)));
disp(sprintf('%-30s: %-15s%-15s%-15s', 'Total', sprintf('%d', size_neutral_total_data), sprintf('%d', size_pathogenic_total_data), sprintf('%d', size_total_data)));

combiner_accuracy          = (combiner_number_of_true_negative + combiner_number_of_true_positive) / (combiner_number_of_true_negative + combiner_number_of_true_positive + combiner_number_of_false_negative + combiner_number_of_false_positive);
combiner_sensitivity       = combiner_number_of_true_positive / (combiner_number_of_true_positive + combiner_number_of_false_negative);
combiner_specificity       = combiner_number_of_true_negative / (combiner_number_of_true_negative + combiner_number_of_false_positive);
combiner_balanced_accuracy = (combiner_sensitivity + combiner_specificity) / 2;

combiner_balanced_error    = (combiner_neutral_error + combiner_pathogenic_error) / 2;

condel_accuracy            = (condel_number_of_true_negative + condel_number_of_true_positive) / (condel_number_of_true_negative + condel_number_of_true_positive + condel_number_of_false_negative + condel_number_of_false_positive);
condel_sensitivity         = condel_number_of_true_positive / (condel_number_of_true_positive + condel_number_of_false_negative);
condel_specificity         = condel_number_of_true_negative / (condel_number_of_true_negative + condel_number_of_false_positive);
condel_balanced_accuracy   = (condel_sensitivity + condel_specificity) / 2;

condel_balanced_error      = (condel_neutral_error + condel_pathogenic_error) / 2;

disp(' ');
disp('************ Training configuration *************');
disp(sprintf('%-30s: %-15s', 'step size', sprintf('%6.4f', step_size)));
disp(sprintf('%-30s: %-15s', 'iteration', sprintf('%d', epoch_time)));
disp(sprintf('%-30s: %2d to %2d', '#hidden nodes in 1st layer', min(hidden_nodes_A_layer), max(hidden_nodes_A_layer)));
disp(sprintf('%-30s: %2d to %2d', '#hidden nodes in 2nd layer', min(hidden_nodes_B_layer), max(hidden_nodes_B_layer)));
disp(' ');
disp('****** Precision measurement configuration ******');
disp(sprintf('%-30s: %d', 'neutral samples', number_of_neutral_samples));
disp(sprintf('%-30s: %d', 'pathogenic samples', number_of_pathogenic_samples));
disp(sprintf('%-30s: %6.4f', 'pathogenic ratio', pathogenic_ratio));
disp(' ');
disp('*********** Type I and Type II errors ***********');
disp(sprintf('%32s%-15s%-15s', '', 'Combiner', 'Condel'));
disp(sprintf('%-30s: %-15s%-15s', 'false positive samples', sprintf('%d', condel_number_of_false_positive), sprintf('%d', combiner_number_of_false_positive)));
disp(sprintf('%-30s: %-15s%-15s', 'true negative samples', sprintf('%d', condel_number_of_true_negative), sprintf('%d', combiner_number_of_true_negative)));
disp(sprintf('%-30s: %-15s%-15s', 'false negative samples', sprintf('%d', condel_number_of_false_negative), sprintf('%d', combiner_number_of_false_negative)));
disp(sprintf('%-30s: %-15s%-15s', 'true positive samples', sprintf('%d', condel_number_of_true_positive), sprintf('%d', combiner_number_of_true_positive)));
disp(' ');
disp('******************* Precision *******************');
disp(sprintf('%32s%-15s%-15s', '', 'Combiner', 'Condel'));
disp(sprintf('%-30s: %-15s%-15s', 'accuracy', sprintf('%6.4f', condel_accuracy), sprintf('%6.4f', combiner_accuracy)));
disp(sprintf('%-30s: %-15s%-15s', 'sensitivity', sprintf('%6.4f', condel_sensitivity), sprintf('%6.4f', combiner_sensitivity)));
disp(sprintf('%-30s: %-15s%-15s', 'specificity', sprintf('%6.4f', condel_specificity), sprintf('%6.4f', combiner_specificity)));
disp(sprintf('%-30s: %-15s%-15s', 'balanced accuracy', sprintf('%6.4f', condel_balanced_accuracy), sprintf('%6.4f', combiner_balanced_accuracy)));
disp(' ');
disp('************** Probabilistic error **************');
disp(sprintf('%32s%-15s%-15s', '', 'Combiner', 'Condel'));
disp(sprintf('%-30s: %-15s%-15s', 'Neutral', sprintf('%6.4f', condel_neutral_error), sprintf('%6.4f', combiner_neutral_error)));
disp(sprintf('%-30s: %-15s%-15s', 'Pathogenic', sprintf('%6.4f', condel_pathogenic_error), sprintf('%6.4f', combiner_pathogenic_error)));
disp(sprintf('%-30s: %-15s%-15s', 'Balanced', sprintf('%6.4f', condel_balanced_error), sprintf('%6.4f', combiner_balanced_error)));


