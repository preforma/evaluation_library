pool_new = pool;
run_delete = pool;
run_delete([runSet.Properties.RowNames],:) = [];
run_delete = pool;
run_delete([runSet.Properties.RowNames],:) = [];
pool_new([run_delete.Properties.RowNames],:) = [];
[accuracy, AUC, LAM, consistency] = PREFORMA_Measures(pool_new, runSet);
a = vertcat(accuracy{:,1});
mean(cell2mat(a))


b = vertcat(AUC{:,1});
mean(cell2mat(b))


stem(vertcat(accuracy{:,1}{:}),'filled', 'LineWidth', 1);
hold on;
set(gca,'Xticklabel',accuracy.Properties.RowNames,'XTick',1:numel(accuracy.Properties.RowNames));
xtickangle(gca, 45);
xlim([0 numel(AUC.Properties.RowNames)+1]);
ylim([0 1.05]);

xlabel('Classes', 'FontSize', 16); % x-axis label
ylabel('Accuracy', 'FontSize', 16); % y-axis label
title('Accuracy of text media type classes', 'FontSize', 20);
mu =mean(cell2mat(a));
hline = refline([0 mu]);
hline.Color = 'r';


figure();
stem(vertcat(AUC{:,1}{:}),'filled', 'LineWidth', 1);
hold on;
set(gca,'Xticklabel',AUC.Properties.RowNames,'XTick',1:numel(AUC.Properties.RowNames));
xtickangle(gca, 45);
xlim([0 numel(AUC.Properties.RowNames)+1]);
ylim([0 1.05]);
xlabel('Classes', 'FontSize', 16); % x-axis label
ylabel('AUC', 'FontSize', 16); % y-axis label
title('AUC of text media type classes', 'FontSize', 20);
mu =mean(cell2mat(a));
hline = refline([0 mu]);
hline.Color = 'r';
