%% Setup
clear; home;

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
ni_all=linspace(1, 2, 101);
no_all=linspace(1, 2, length(ni_all));
ni_ass=1.3;
no_ass_all=[1, 2];
% no_ass=2;

IerrFrac=0.001;
Perr=0.1*pi/180;

kappa=1.2e-3;

%% Run Fits
optsCost=opts;
optsCost.nin=ni_ass;
optsCost.L=L;
optsCost.kappa=kappa;

muaRec=NaN(length(ni_all), 2, length(mua_ex), length(musp_ex), length(no_ass_all));
muspRec=NaN(length(ni_all), 2, length(mua_ex), length(musp_ex), length(no_ass_all));
for noAss_ind=1:length(no_ass_all)

    no_ass=no_ass_all(noAss_ind);
    optsCost.nout=no_ass;
    
    for muaInd=1:length(mua_ex)
        for muspInd=1:length(musp_ex)
            
            optProp.mua=mua_ex(muaInd);
            optProp.musp=musp_ex(muspInd);
            
            for InOutInd=1:2
                for nInd=1:length(ni_all)
                    switch InOutInd
                        case 1
                            optProp.nin=ni_all(nInd);
                            optProp.nout=no_ass;
                        case 2
                            optProp.nin=ni_ass;
                            optProp.nout=no_all(nInd);
                        otherwise
                    end
                    
                    T=NaN(size(rDet, 1), 1);
                    for i=1:size(rDet, 1)
                        T(i)=Tcuv(rSrc, rDet(i, :), L, optProp, opts);
                    end
                    Terr=(abs(T)*IerrFrac)*exp(1i*Perr);
                    Tmeas=T;
                    
                    mu0=[optProp.mua, optProp.musp];
                    
                    options=optimoptions('fmincon',...
                        'Algorithm', 'interior-point');
                    
                    mu=fmincon(@cuvTcost,...
                            mu0, [], [], [], [], [0, 0], [], [], options,...
                            Tmeas, Terr, rSrc, rDet, optsCost);
                    muaRec(nInd, InOutInd, muaInd, muspInd, noAss_ind)=mu(1);
                    muspRec(nInd, InOutInd, muaInd, muspInd, noAss_ind)=mu(2);
                end
            end
        end
    end
end

%% Plot
figure(110); clf;
for noAss_ind=1:length(no_ass_all)
    cols=hsv(length(mua_ex)*length(musp_ex));
    muaLim=[-0.001, 0.09];
    muspLim=[0, 3.5];
    nLim=[0.98, 2.02];
    
    legs={}; n=0;
    for i=1:length(mua_ex)
        for j=1:length(musp_ex)
            n=n+1;
            strTmp=sprintf(...
                '$\\mu_{a,true}=%.3f$ mm$^{-1}$\n$\\mu''_{s,true}=%.2f$ mm$^{-1}$',...
                mua_ex(i), musp_ex(j));
            legs=[legs; strTmp];
        end
    end
    legs=[legs; 'Assumed $n$';...
        sprintf('$n_{i,assumed}=%.1f$\n$n_{o,assumed}=%.1f$', ni_ass, no_ass_all(1));...
        sprintf('$n_{i,assumed}=%.1f$\n$n_{o,assumed}=%.1f$', ni_ass, no_ass_all(2));];
    
    for InOutInd=2:-1:1
        switch InOutInd
            case 1
                n_ass=ni_ass;
                n_all=ni_all;
            case 2
                n_ass=no_ass_all(noAss_ind);
                n_all=no_all;
            otherwise
        end

        switch noAss_ind
            case 1
                sym='-';
            case 2
                sym='--';
            otherwise
        end

        subaxis(2, 2, 1+(InOutInd-1), 's', 0.03,...
            'mr', 0.02, 'mb', 0.07, 'mt', 0.25);
        n=0; h=[];
        for muaInd=1:length(mua_ex)
            for muspInd=1:length(musp_ex)
                n=n+1;
                h(n)=fill(NaN, NaN, cols(n, :));
                plot(n_all,...
                    muaRec(:, InOutInd, muaInd, muspInd, noAss_ind), sym,...
                    'Color', cols(n, :)); hold on;
            end
        end
        xline(n_ass, [sym 'k']);
        h(end+1)=fill(NaN, NaN, [0, 0, 0]);
        h(end+1)=plot(NaN, NaN, '-', 'Color', [0.5, 0.5, 0.5]);
        h(end+1)=plot(NaN, NaN, '--', 'Color', [0.5, 0.5, 0.5]);
        xlim(nLim);
        ylim(muaLim);
        if noAss_ind==2
            ax=gca;
            ax.XTickLabel={};
            switch InOutInd
                case 1
                    ylabel('$\mu_a$ (mm$^{-1}$)');
                    lg=legend(h, legs,...
                        'NumColumns', 3,...
                        'Location', 'northwest');
                    set(lg, 'Position', lg.Position+[-0.15, 0.28, 0, 0]);
                    text(1.03, 0.085, '\textbf{(a)}', 'Interpreter', 'latex');
                case 2
                    ax.YTickLabel={};
                    text(1.03, 0.085, '\textbf{(b)}', 'Interpreter', 'latex');
                otherwise
            end
            hold off;
        end
    
        subaxis(2, 2, 3+(InOutInd-1));
        n=0;
        for muaInd=1:length(mua_ex)
            for muspInd=1:length(musp_ex)
                n=n+1;
                plot(n_all,...
                    muspRec(:, InOutInd, muaInd, muspInd, noAss_ind), sym,...
                    'Color', cols(n, :)); hold on;
            end
        end
        xline(n_ass, [sym 'k']);
        xlim(nLim);
        ylim(muspLim);
        if noAss_ind==2
            ax=gca;
            switch InOutInd
                case 1
                    ylabel('$\mu''_s$ (mm$^{-1}$)');
                    xlabel('$n_{i,true}$');
                    text(1.03, 3.25, '\textbf{(c)}', 'Interpreter', 'latex');
                case 2
                    ax.YTickLabel={};
                    xlabel('$n_{o,true}$');
                    text(1.03, 3.25, '\textbf{(d)}', 'Interpreter', 'latex');
                otherwise
            end
            hold off;
        end
    end

end

h=gcf;
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex');
set(findall(gcf,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex');
set(findall(gcf,'-property','FontSize'),'FontSize',10);
set(findall(gcf,'-property','Linewidth'),'Linewidth',1.5);
set(findall(gcf,'-property','Markersize'),'Markersize',6);
set(gcf, 'Units', 'centimeters', 'Innerposition', [1, 1,...
    [13.85, 16]]);
figName=['wrongIndex_2no'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');