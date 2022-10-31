%% Setup 
clear; home;

rng(123);

rSrc=[5.5, 5, 0;...
    39.5, 5, 0];
rDetA=[16.5, 5, 10];
rDetB=[28.5, 5, 10];
rDet=[rDetA; rDetB];

L=[45, 10, 10];

opts.omega=2*pi*100e6; %rad/sec
opts.BC='EBC';
opts.lMax=3;
opts.mMax=3;
opts.nMax=3;

optProp.mua=0.01;
optProp.musp=1;
optProp.nin=1.3;
optProp.nout=1;

%% Run T
Ts=Tcuv(rSrc(1, :), rDet(1, :), L, optProp, opts);
Tl=Tcuv(rSrc(1, :), rDet(2, :), L, optProp, opts);

C=(rand(4, 1)*0.95+0.05).*exp(1i*(rand(4, 1)-0.5));

T1A=C(1)*C(3)*Tcuv(rSrc(1, :), rDet(1, :), L, optProp, opts);
T1B=C(1)*C(4)*Tcuv(rSrc(1, :), rDet(2, :), L, optProp, opts);
T2A=C(2)*C(3)*Tcuv(rSrc(2, :), rDet(1, :), L, optProp, opts);
T2B=C(2)*C(4)*Tcuv(rSrc(2, :), rDet(2, :), L, optProp, opts);

% diffI(muaInd, muspInd)=log(abs(TB)/abs(TA));
% diffP(muaInd, muspInd)=angle(TB)-angle(TA);

%% Plot
cols=lines(2);

figure(120); clf;
subaxis(1, 2, 1, 's', 0.1,...
    'mt', 0.09, 'mb', 0.33, 'mr', 0.03);
bar([log(abs(T1B)/abs(T1A)),...
    mean([log(abs(T1B)/abs(T1A)), log(abs(T2A)/abs(T2B))]),...
    log(abs(T2A)/abs(T2B))], 'FaceColor', 'b'); hold on;
h1=yline(log(abs(Tl)/abs(Ts)), '--k'); 
text(-0.1, -0.7, '\textbf{(a)}', 'Interpreter', 'latex');
hold off;
ylabel('$\ln|\widetilde{T}_{long}|-\ln|\widetilde{T}_{short}|$',...
    'Interpreter', 'latex');
ax=gca;
ax.XTick=1:3;
ax.XTickLabel={...
    '$\ln|$SR$\{\widetilde{T}\}|_{\mathrm{\texttt{1AB}}}$',...
    '$\ln|$DR$\{\widetilde{T}\}|_{\mathrm{\texttt{1AB2}}}$',...
    '$\ln|$SR$\{\widetilde{T}\}|_{\mathrm{\texttt{2BA}}}$'};
ax=gca;
for i=1:length(ax.YTickLabel)
    ax.YTickLabel{i}=['$' ax.YTickLabel{i} '$'];
    ax.YTick(i)=ax.YTick(i);
end
legend(h1, sprintf('Theoretical\nValue'), 'location', 'northeast');

subaxis(1, 2, 2);
bar([angle(T1B)-angle(T1A),...
    mean([angle(T1B)-angle(T1A), angle(T2A)-angle(T2B)]),...
    angle(T2A)-angle(T2B)], 'FaceColor', 'r'); hold on;
yline(angle(Tl)-angle(Ts), '--k'); 
text(-0, 0.34, '\textbf{(b)}', 'Interpreter', 'latex');
hold off;
ylabel('$\angle\widetilde{T}_{long}-\angle\widetilde{T}_{short}$ (rad)',...
    'Interpreter', 'latex');

ax=gca;
ax.XTick=1:3;
ax.XTickLabel={...
    '$\angle$SR$\{\widetilde{T}\}_{\mathrm{\texttt{1AB}}}$',...
    '$\angle$DR$\{\widetilde{T}\}_{\mathrm{\texttt{1AB2}}}$',...
    '$\angle$SR$\{\widetilde{T}\}_{\mathrm{\texttt{2BA}}}$'};
% ax.XTickLabel={...
%     '$\angle\widetilde{T}_{1B}-\angle\widetilde{T}_{1A}$',...
%     'Average',...
%     '$\angle\widetilde{T}_{2A}-\angle\widetilde{T}_{2B}$'};

h=gcf;
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1.5);
set(findall(gcf,'-property','Markersize'),'Markersize',5);
set(gcf, 'Units', 'centimeters', 'Innerposition', [1, 1,...
    [13.85, 5.7]]);
figName=['exampleTcal'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');