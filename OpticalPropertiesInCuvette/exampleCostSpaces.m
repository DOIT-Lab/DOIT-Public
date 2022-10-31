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
ni=1.3;
no=1;

optProp.nin=ni;
optProp.nout=no;

IerrFrac=0.001;
Perr=0.1*pi/180;

kappa=1.2e-3;

%% Example Cost Spaces
mua_all=linspace(0.0025, 0.0225, 100);
musp_all=linspace(0.25, 2.25, 100);

optsCost=opts;
optsCost.nin=ni;
optsCost.nout=no;
optsCost.L=L;
optsCost.kappa=kappa;

cost=NaN(length(mua_all), length(musp_all),...
    length(mua_ex), length(musp_ex));
for muaInd=1:length(mua_ex)
    for muspInd=1:length(musp_ex)
        
        optProp.mua=mua_ex(muaInd);
        optProp.musp=musp_ex(muspInd);
        T=NaN(size(rDet, 1), 1);
        for i=1:size(rDet, 1)
            T(i)=Tcuv(rSrc, rDet(i, :), L, optProp, opts);
        end
        Terr=(abs(T)*IerrFrac)*exp(1i*Perr);
        
        for i=1:length(mua_all)
            for j=1:length(musp_all)
                cost(i, j, muaInd, muspInd)=...
                    cuvTcost([mua_all(i), musp_all(j)],...
                    T, Terr, rSrc, rDet, optsCost);
            end
        end
    end
end

%% Plot
% cl=[min(cost(:)), max(cost(:))];
cl=[1e-2, 1e4];
cTicks=[1e-1, 1e0, 1e1, 1e2, 1e3];

figure(80); clf; n=0;
for muaInd=1:length(mua_ex)
    for muspInd=1:length(musp_ex)
        n=n+1;
        subaxis(length(mua_ex), length(musp_ex), n,...
            'sv', 0.08, 'sh', 0.03,...
            'mt', 0.08, 'ml', 0.09, 'mr', 0.14, 'mb', 0.08);
        pcolor(mua_all, musp_all, cost(:, :, muaInd, muspInd).'); hold on;
        contour(mua_all, musp_all, cost(:, :, muaInd, muspInd).',...
            cTicks, 'w'); 
        h1=xline(mua_ex(muaInd), '--k');
        yline(musp_ex(muspInd), '--k');
        hold off;
        shading flat;
        set(gca, 'ColorScale', 'log');
        caxis(cl);
        title(sprintf(...
            ['$\\mu_{a,true}=%.3f$ mm$^{-1}$\n'...
            '\\textbf{(%s)} $\\mu''_{s,true}=%.1f$ mm$^{-1}$'],...
            mua_ex(muaInd), char(96+n), musp_ex(muspInd)),...
            'Interpreter', 'latex');
        
        if n==1
            legend(h1, 'True',...
                'Location', 'northeast');
        end
        if muspInd==3
            pos=get(gca, 'Position');
            cb=colorbar;
            set(gca, 'Position', pos);
            set(cb, 'Ticks', cTicks);
            ylabel(cb, '$\chi^2$', 'Interpreter', 'latex');
        end
        if muspInd==1
            ylabel('$\mu''_s$ (mm$^{-1}$)');
        else
            set(gca, 'YTickLabel', {});
        end
        if muaInd==3
            xlabel('$\mu_a$ (mm$^{-1}$)');
        else
            set(gca, 'XTickLabel', {});
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
    [13.85, 14]]);
figName=['costSpaces'];
exportgraphics(h, [figName '.eps'], 'ContentType', 'vector');

%% Functions
function kappaCostVal = kappaCost(kappa, opts)

    opts.kappa=kappa;
    
    mua_all=opts.mua_all;
    musp_all=opts.musp_all;
    Tmeas=opts.Tmeas;
    Terr=opts.Terr;
    rSrc=opts.rSrc;
    rDet=opts.rDet;
    
    cost=NaN(length(opts.mua_all), length(opts.musp_all));
    n_mua=length(opts.mua_all);
    n_musp=length(opts.musp_all);
    parfor i=1:n_mua
        for j=1:n_musp
            cost(i, j)=cuvTcost([mua_all(i), musp_all(j)],...
                Tmeas, Terr, rSrc, rDet, opts);
        end
    end
    costNorm=cost/max(cost(:));
    
    [II, JJ]=meshgrid(linspace(0, 1, length(mua_all)),...
        linspace(0, 1, length(musp_all)));
    [III, JJJ]=meshgrid(...
        linspace(0, 1, 1000),...
        linspace(0, 1, 1000));
    cost_up=interp2(II, JJ, costNorm, III, JJJ, 'spline');

    figure(999); clf;
    pcolor(III(1, :), JJJ(:, 1), cost_up.'); hold on;
    M=contour(III(1, :), JJJ(:, 1), cost_up.',...
        [1, 1]*quantile(cost_up(:), 0.05), 'w'); hold off; shading flat;
    axis equal tight off;
    shp=alphaShape(M(1, 2:M(2, 1)).', M(2, 2:M(2, 1)).');

    kappaCostVal=perimeter(shp)^2/area(shp);

    title(sprintf('$\\kappa=%f$\n$P^2/A=%f$', kappa, kappaCostVal),...
        'Interpreter', 'latex');
    drawnow;
end