%% Setup
clear; home;

nNoise=101;

rSrc=[5.5, 5, 0];
rDetA=[16.5, 5, 10];
rDetB=[28.5, 5, 10];
rDet=[rDetA; rDetB];

L=[45, 10, 10];

opts.omega=2*pi*100e6; %rad/sec
opts.BC='EBC';
opts.lMax=3;
opts.mMax=3;
opts.nMax=3;

mua_ex=[0.005, 0.01, 0.02];
musp_ex=[0.5, 1, 2];
ni=1.3;
no=1;

optProp.nin=ni;
optProp.nout=no;

IerrFrac=0.001;
Perr=0.1*pi/180;

kappa=1.2e-3;

%% Run Fits
optsCost=opts;
optsCost.nin=ni;
optsCost.nout=no;
optsCost.L=L;
optsCost.kappa=kappa;

optsCostNOkappa=optsCost;
optsCostNOkappa.kappa=1;

% mu0=[0.01, 1];

muaRec=NaN(nNoise, length(mua_ex), length(musp_ex));
muspRec=NaN(nNoise, length(mua_ex), length(musp_ex));
% muaRecNOkappa=NaN(nNoise, length(mua_ex), length(musp_ex));
% muspRecNOkappa=NaN(nNoise, length(mua_ex), length(musp_ex));
for muaInd=1:length(mua_ex)
    for muspInd=1:length(musp_ex)
        
        optProp.mua=mua_ex(muaInd);
        optProp.musp=musp_ex(muspInd);
        
        T=NaN(size(rDet, 1), 1);
        for i=1:size(rDet, 1)
            T(i)=Tcuv(rSrc, rDet(i, :), L, optProp, opts);
        end
        Terr=(abs(T)*IerrFrac)*exp(1i*Perr);
        
        for noiseInd=1:nNoise
%             Tmeas=(abs(T)+randn(size(T)).*abs(T)*IerrFrac).*...
%                 exp(1i*(angle(T)+randn(size(T))*Perr));
            Tmeas1=(abs(T)+randn(size(T)).*abs(Terr)).*...
                exp(1i*(angle(T)+randn(size(T)).*angle(Terr)));
            Tmeas2=(abs(T)+randn(size(T)).*abs(Terr)).*...
                exp(1i*(angle(T)+randn(size(T)).*angle(Terr)));
            Tmeas=[Tmeas1, Tmeas2];
            
            mu0=[(mua_ex(2)-mua_ex(1))*rand+mua_ex(1),...
                (musp_ex(2)-musp_ex(1))*rand+musp_ex(1)];
            
            options=optimoptions('fmincon',...
                'Algorithm', 'interior-point');
            
            mu=fmincon(@cuvTcost,...
                    mu0, [], [], [], [], [0, 0], [], [], options,...
                    Tmeas, Terr, rSrc, rDet, optsCost);
            muaRec(noiseInd, muaInd, muspInd)=mu(1);
            muspRec(noiseInd, muaInd, muspInd)=mu(2);
            
%             mu=fmincon(@cuvTcost,...
%                     mu0, [], [], [], [], [0, 0], [], [], options,...
%                     Tmeas, Terr, rSrc, rDet, optsCostNOkappa);
%             muaRecNOkappa(noiseInd, muaInd, muspInd)=mu(1);
%             muspRecNOkappa(noiseInd, muaInd, muspInd)=mu(2);
        end
    end
end

%% Plot
cols=hsv(length(mua_ex)*length(musp_ex));

grp={}; n=0;
for i=1:length(mua_ex)
    for j=1:length(musp_ex)
        n=n+1;
        strTmp=sprintf(...
            '$\\mu_{a,true}=%.3f$ mm$^{-1}$\n$\\mu''_{s,true}=%.2f$ mm$^{-1}$',...
            mua_ex(i), musp_ex(j));
        grp=[grp; repmat({strTmp}, size(muaRec, 1), 1)];
    end
end

figure(100); clf;
h1=scatterhist(muaRec(:), muspRec(:),...
    'Group', grp, 'Kernel', 'on',...
    'Color', cols, 'LineStyle', '-', 'LineWidth', 2,...
    'Marker', '.'); hold on;
text(-0.01, 2.2, '\textbf{(a)}');
text(0.002, 2.1, '\textbf{(b)}');
text(0.001, -0.15, '\textbf{(c)}');
xlim([0, 0.026]);
ylabel('$\mu''_s$ (mm$^{-1}$)');
xlabel('$\mu_a$ (mm$^{-1}$)');

title(sprintf(...
    ['$\\sigma_{|\\widetilde{T}|}/|\\widetilde{T}|=%.3f$, '...
    '$\\sigma_{\\angle\\widetilde{T}}=%.1f$ mrad'],...
    IerrFrac, 1000*Perr));

h1(1).Position(2)=h1(1).Position(2)+0.06;
h1(1).Position(3)=h1(1).Position(3)+0.09;

h=gcf;
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1.5);
set(findall(gcf,'-property','Markersize'),'Markersize',6);
set(gcf, 'Units', 'centimeters', 'Innerposition', [1, 1,...
    [13.85, 18.5]]);
figName=['exampleErrors'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');

%% Table
fprintf('\n');
for i=1:length(mua_ex)
    for j=1:length(musp_ex)
%         tmp=cov(muaRec(:, i, j), muspRec(:, i, j));
        tmp2=corrcoef(muaRec(:, i, j), muspRec(:, i, j));
        fprintf('%.4f & %.2f & %.5f & %.3f & %.3f & %.3f & %.5f\\\\\n',...
            mua_ex(i), musp_ex(j),...
            std(muaRec(:, i, j)), std(muaRec(:, i, j))/mean(muaRec(:, i, j)),...
            std(muspRec(:, i, j)), std(muspRec(:, i, j))/mean(muspRec(:, i, j)),...
            tmp2(1, 2));
    end
end