%% Setup
clear; home;

rho=[25, 35]; %mm

% X=[mua1, mua2, musp1, musp2, ztop] (1/mm, 1/mm, 1/mm, 1/mm, mm)
X=[...
    0.010, 0.010, 1.0, 1.0, 10;...
    0.010, 0.010, 0.5, 1.5, 10;...
    0.010, 0.010, 1.5, 0.5, 10;...
    0.005, 0.015, 1.0, 1.0, 10;...
    0.015, 0.005, 1.0, 1.0, 10];

[Y]=twoLayEffHomoOptProp(X, rho);

%% Plot
figure(1); clf;
subplot(3, 1, 1);
bar([X(:, 1:2), Y(:, 1)]);
xlim([0.5, 5.5]);
set(gca, 'XTickLabel', {});
ylabel('\mu_a (mm^{-1})');
legend('Top', 'Bottom', 'Homo Rec');

subplot(3, 1, 2);
bar([X(:, 3:4), Y(:, 2)]);
xlim([0.5, 5.5]);
set(gca, 'XTickLabel', {});
ylabel('\mu''_s (mm^{-1})');

subplot(3, 1, 3);
bar(X(:, 5), 'k');
xlim([0.5, 5.5]);
xlabel('Case');
ylabel('z_{top}');
