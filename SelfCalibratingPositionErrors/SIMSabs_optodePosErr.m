%% Setup
clear; home;
svFig=true;

baseFigName='PHOTC5_optPosErr';

%% Fig 100 - Types of Sets
coors_tmp=load("Coors_mHEX1p2.mat");

DSnms={ ...
    'Linear (LINR)';
    'Asymmetric-Linear (ALIN)';
    'Trapezoidal (TRAP)';
    'Diagonal-Rectangular (DRCT)'
    };
DSnms_short={ ...
    'LINR';
    'ALIN';
    'TRAP';
    'DRCT'
    };
% LINR
rd_all(:, :, 1)=[ ...
    -6, 0, 0; ...
    6, 0, 0;
    ];
rs_all(:, :, 1)=[ ...
    -31, 0, 0; ...
    31, 0, 0;
    ];
% ALIN
rd_all(:, :, 2)=[ ...
    -6, 0, 0; ...
    6, 0, 0;
    ];
rs_all(:, :, 2)=[ ...
    -26, 0, 0; ...
    36, 0, 0;
    ];
% TRAP
rs_tmp=coors_tmp.AllDets;
rd_tmp=coors_tmp.AllSrcs(1:2, :);
rd_all(:, :, 3)=rd_tmp-mean([rs_tmp; rd_tmp]);
rs_all(:, :, 3)=rs_tmp-mean([rs_tmp; rd_tmp]);
rs_all(:, 1, 3)=-rs_all(:, 1, 3);
rd_all(:, 1, 3)=-rd_all(:, 1, 3);
% DRCT
rd_all(:, :, 4)=[ ...
    -18.5, 12.5, 0; ...
    18.5, -12.5, 0;
    ];
rs_all(:, :, 4)=[ ...
    -18.5, -12.5, 0;
    18.5, 12.5, 0; ...
    ];

% FIGURE
h=figure(100); clf;
h.Name=baseFigName;
xl=[-1, 1]* ...
    (max(abs([squeeze(rd_all(:, 1, :)); squeeze(rs_all(:, 1, :))]), ...
    [], 'all')+3);
yl=[-1, 1]* ...
    (max(abs([squeeze(rd_all(:, 2, :)); squeeze(rs_all(:, 2, :))]), ...
    [], 'all')+3);

for i=1:size(rs_all, 3)
    subaxis(2, 2, i);
    plot(rd_all(:, 1, i), rd_all(:, 2, i), 'ob'); hold on;
    plot(rs_all(:, 1, i), rs_all(:, 2, i), 'sr');
    hold off;
    xlim(xl);
    ylim(yl);
    title(sprintf('(%s) %s', ...
        char(i+96), DSnms{i}));
end
drawnow;
subaxis(2, 2, 1, 'sv', 0.07);
set(gca,'XTickLabel',{});
ylabel('$y$ (mm)', 'Interpreter','latex');
legend('Detectors', 'Sources', ...
    'Orientation','horizontal');
axis equal;
subaxis(2, 2, 2);
set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});
axis equal;
subaxis(2, 2, 3);
xlabel('$x$ (mm)', 'Interpreter','latex');
ylabel('$y$ (mm)', 'Interpreter','latex');
axis equal;
subaxis(2, 2, 4);
xlabel('$x$ (mm)', 'Interpreter','latex');
set(gca,'YTickLabel',{});
axis equal;

svFig_params.FigSz_cm=[13.5, 7];
svFig_params.Markersize=5;
saveFig; clear svFig_params;



%% Fig 200 - LINR-1/A-lin/iter
setInd=1;

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

cols=turbo(100);
cols(1, :)=1;
cols(end, :)=1;

h=figure(200); clf;
h.Name=baseFigName;



%######################--- 1 ---###########################################
% Move Src 1
cl_mua=[-1, 1]*20;
cl_musp=[-1, 1]*20;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=(-10:0.1:10);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
mua_RelErr_lin=NaN(length(dx), length(dy));
musp_RelErr_lin=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
        [mua_RelErr_lin(i, j), musp_RelErr_lin(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "lin");
    end
end



subaxis(6, 2, [1, 2], ...
    'sv', 0.04, 'sh', 0, ...
    'mt', 0.01, 'mb', 0.07, 'mr', 0.1, 'ml', 0.07);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rs_all(1, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rs_all(1, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*70);
pos=get(gca, 'Position');
legend('Detectors', 'Sources', ...
    'Orientation','vertical', 'Location','northeast', ...
    'Color', [0.75, 0.75, 0.75], ...
    'EdgeColor',[0.25, 0.25, 0.25]);
set(gca, 'Position', pos);
title(['(a) ' DSnms{setInd} ' Arrangement: Source 1 Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);



subaxis(6, 2, 3);
pcolor(dx, dy, 100*mua_RelErr_lin.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_lin.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
set(gca, 'XTickLabel', {});
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(b) Slopes Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(6, 2, 4);
pcolor(dx, dy, 100*musp_RelErr_lin.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_lin.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
set(gca, 'XTickLabel', {});
set(gca, 'YTickLabel', {});
title('(c) Slopes Method $\mu_s^\prime$ Error', 'Interpreter','latex');





subaxis(6, 2, 5);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(d) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(6, 2, 6);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(e) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');











%######################--- A ---###########################################
% Move Det A
cl_mua=[-1, 1]*40;
cl_musp=[-1, 1]*100;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=(-14:0.1:6);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
mua_RelErr_lin=NaN(length(dx), length(dy));
musp_RelErr_lin=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(1, :)=rd_ass(1, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
        [mua_RelErr_lin(i, j), musp_RelErr_lin(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "lin");
    end
end

subaxis(6, 2, [7, 8]);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rd_all(1, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rd_all(1, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*70);
title(['(f) ' DSnms{setInd} ' Arrangement: Detector A Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);



subaxis(6, 2, 9);
pcolor(dx, dy, 100*mua_RelErr_lin.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_lin.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
set(gca, 'XTickLabel', {});
ylabel('$\Delta y_{A}$ (mm)', ...
    'Interpreter','latex');
title('(g) Slopes Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(6, 2, 10);
pcolor(dx, dy, 100*musp_RelErr_lin.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_lin.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
set(gca, 'XTickLabel', {});
set(gca, 'YTickLabel', {});
title('(h) Slopes Method $\mu_s^\prime$ Error', 'Interpreter','latex');



subaxis(6, 2, 11);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{A}$ (mm)', ...
    'Interpreter','latex');
title('(i) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(6, 2, 12);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(j) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');



svFig_params.FigSz_cm=[13.5, 18];
svFig_params.doPDF=false;
svFig_params.doPDFraster=true;
saveFig; clear svFig_params;

%% TABLE - Absolute systematic error
musp_all=0.5:0.5:1.5; %1/mm
mua_all=0.005:0.005:0.015; %1/mm

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=[]; %1/mm
optProp_act.mua=[]; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

fprintf(['rhos \t mua \t musp \t ' ...
    'mua err frac lin \t musp err frac lin \t ' ...
    'mua err frac iter \t musp err frac iter \n']);
for setInd=1:2 ...size(rs_all, 3)
    rhos=sort([vecnorm(rs_all(:, :, setInd)-rd_all(:, :, setInd), 2, 2);
        vecnorm(rs_all(:, :, setInd)-rd_all(end:-1:1, :, setInd), 2, 2)]);

    rsPen_ass=rs_all(:, :, setInd);
    rd_ass=rd_all(:, :, setInd);

    rsPen_act=rsPen_ass;
    rd_act=rd_ass;

    for muaInd=1:length(mua_all)
        for muspInd=1:length(musp_all)
            
            optProp_act.mua=mua_all(muaInd);
            optProp_act.musp=musp_all(muspInd);
                    
            fprintf('{[}%.0f, %.0f, %.0f, %.0f{]} & ', rhos);
            fprintf('%.3f & %.1f & ', ...
                optProp_act.mua, optProp_act.musp);
            
            [mua_RelErr_lin, musp_RelErr_lin, ...
                ~, ~, ~]=DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, optProp_act, ...
                dmua_act, omega, 'lin');

            [mua_RelErr_iter, musp_RelErr_iter, ...
                ~, ~, ~]=DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, optProp_act, ...
                dmua_act, omega, 'iter');
            
            dec1=-(floor(log10(abs(mua_RelErr_lin*100)))-1);
            dec2=-(floor(log10(abs(musp_RelErr_lin*100)))-1);
            dec3=-(floor(log10(abs(mua_RelErr_iter*100)))-1);
            dec4=-(floor(log10(abs(musp_RelErr_iter*100)))-1);
            formatSpec=[
                '%.' num2str(dec1) 'f & ' ...
                '%.' num2str(dec2) 'f & ' ...
                '%.' num2str(dec3) 'f & ' ...
                '%.' num2str(dec4) 'f ' ...
                '\\\\ \n'];
            fprintf(formatSpec, ...
                mua_RelErr_lin*100, musp_RelErr_lin*100, ...
                mua_RelErr_iter*100, musp_RelErr_iter*100);
        end
    end
end

%% TABLE - Absolute RMSE-orbit
musp_all=0.5:0.5:1.5; %1/mm
mua_all=0.005:0.005:0.015; %1/mm

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=[]; %1/mm
optProp_act.mua=[]; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

rOrb=1; %mm

fprintf(['mua \t musp \t optChar \t' ...
    'lin mua_CV \t lin musp_CV \t ' ...
    'iter mua_CV \t iter musp_CV \n']);
for setInd=1:size(rs_all, 3)
    fprintf('%s \n', DSnms_short{setInd});
    if setInd==2
        optInds=1:2;
    % elseif setInd==4
    %     continue;
    else
        optInds=1;
    end
    
    rsPen_ass=rs_all(:, :, setInd);
    rd_ass=rd_all(:, :, setInd);

    for muaInd=1:length(mua_all)
        for muspInd=1:length(musp_all)
            
            optProp_act.mua=mua_all(muaInd);
            optProp_act.musp=musp_all(muspInd);
    
            for SD='SD'
                if setInd==4 && SD=='D'
                    continue;
                end
                
                for optInd=optInds
                    optChar=sprintf('%s%.0f', SD, optInd);
                    
                    switch SD
                        case 'S'
                            if setInd==1
                                optNm=sprintf('%.0f or %.0f', optInd, optInd+1);
                            elseif setInd==3
                                % optNm=sprintf( ...
                                %     '(%.0f/%.0f)\\textsubscript{TRAP}', ...
                                %     optInd, optInd+1);
                                optNm=sprintf( ...
                                    '%.0f or %.0f', ...
                                    optInd, optInd+1);
                            elseif setInd==4
                                optNm=sprintf('%.0f, %.0f, %s, or %s', ...
                                    optInd, optInd+1, ...
                                    char(optInd+64), char(optInd+64+1));
                            else
                                optNm=sprintf('%.0f', optInd);
                            end
                        case 'D'
                            if setInd==1
                                optNm=sprintf('%s or %s', ...
                                    char(optInd+64), char(optInd+64+1));
                            elseif setInd==3
                                % optNm=sprintf( ...
                                %     ['(%s/%s)\\textsubscript{TRAP}/' ...
                                %     '(1/2/A/B)\\textsubscript{DRCT}'], ...
                                %     char(optInd+64), char(optInd+64+1));
                                optNm=sprintf( ...
                                    '%s or %s', ...
                                    char(optInd+64), char(optInd+64+1));
                            else
                                optNm=sprintf('%s', char(optInd+64));
                            end
                        otherwise
                            error;
                    end
                    
                    fprintf('%.3f & %.1f & ', ...
                        optProp_act.mua, optProp_act.musp);
                    fprintf('%s & ', optNm);
                    
                    if setInd==1
                        [mua_RMSE_O_lin, musp_RMSE_O_lin, mua0, musp0]= ...
                            posOrbitRMSE({optChar}, ...
                            rsPen_ass, rd_ass, rOrb, ...
                            optProp_act, dmua_act, omega, 'lin');
                        mua_CV_O_lin=mua_RMSE_O_lin/mua0;
                        musp_CV_O_lin=musp_RMSE_O_lin/musp0;
                    end

                    [mua_RMSE_O_iter, musp_RMSE_O_iter, mua0, musp0]= ...
                        posOrbitRMSE( ...
                        {optChar}, ...
                        rsPen_ass, rd_ass, rOrb, ...
                        optProp_act, dmua_act, omega, 'iter');
                    mua_CV_O_iter=mua_RMSE_O_iter/mua0;
                    musp_CV_O_iter=musp_RMSE_O_iter/musp0;
                    
                    if setInd==1
                        dec1=-(floor(log10(abs(mua_CV_O_lin*100)))-1);
                        dec2=-(floor(log10(abs(musp_CV_O_lin*100)))-1);
                        dec3=-(floor(log10(abs(mua_CV_O_iter*100)))-1);
                        dec4=-(floor(log10(abs(musp_CV_O_iter*100)))-1);
                        formatSpec=[
                            '%.' num2str(dec1) 'f & ' ...
                            '%.' num2str(dec2) 'f & ' ...
                            '%.' num2str(dec3) 'f & ' ...
                            '%.' num2str(dec4) 'f ' ...
                            '\\\\ \n'];
                        fprintf(formatSpec, ...
                            mua_CV_O_lin*100, musp_CV_O_lin*100, ...
                            mua_CV_O_iter*100, musp_CV_O_iter*100);
                    else
                        dec3=-(floor(log10(abs(mua_CV_O_iter*100)))-1);
                        dec4=-(floor(log10(abs(musp_CV_O_iter*100)))-1);
                        formatSpec=[
                            '%.' num2str(dec3) 'f & ' ...
                            '%.' num2str(dec4) 'f ' ...
                            '\\\\ \n'];
                        fprintf(formatSpec, ...
                            mua_CV_O_iter*100, musp_CV_O_iter*100);
                    end
                end
            end
        end
    end
end





%% Fig 300 - ALIN-1/2/A/B-iter

setInd=2;

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

cols=turbo(100);
cols(1, :)=1;
cols(end, :)=1;

h=figure(300); clf;
h.Name=baseFigName;



%######################--- 1 ---###########################################
% Move Src 1
cl_mua=[-1, 1]*20;
cl_musp=[-1, 1]*20;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=(-10:0.1:10);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
    end
end



subaxis(8, 2, [1, 2], ...
    'sv', 0.04, 'sh', 0, ...
    'mt', 0.01, 'mb', 0.07, 'mr', 0.1, 'ml', 0.07);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rs_all(1, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rs_all(1, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*80);
pos=get(gca, 'Position');
legend('Detectors', 'Sources', ...
    'Orientation','vertical', 'Location','northeast', ...
    'Color', [0.75, 0.75, 0.75], ...
    'EdgeColor',[0.25, 0.25, 0.25]);
set(gca, 'Position', pos);
title(['(a) ' DSnms{setInd} ' Arrangement: Source 1 Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);



subaxis(8, 2, 3);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(b) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(8, 2, 4);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(c) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');











%######################--- A ---###########################################
% Move Det A
cl_mua=[-1, 1]*50;
cl_musp=[-1, 1]*100;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=(-14:0.1:6);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
mua_RelErr_lin=NaN(length(dx), length(dy));
musp_RelErr_lin=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(1, :)=rd_ass(1, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
        [mua_RelErr_lin(i, j), musp_RelErr_lin(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "lin");
    end
end

subaxis(8, 2, [5, 6]);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rd_all(1, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rd_all(1, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*80);
title(['(d) ' DSnms{setInd} ' Arrangement: Detector A Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);



subaxis(8, 2, 7);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{A}$ (mm)', ...
    'Interpreter','latex');
title('(e) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(8, 2, 8);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(f) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');






%######################--- B ---###########################################
% Move Det B
cl_mua=[-1, 1]*50;
cl_musp=[-1, 1]*100;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=(-6:0.1:14);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
mua_RelErr_lin=NaN(length(dx), length(dy));
musp_RelErr_lin=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(2, :)=rd_ass(2, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
        [mua_RelErr_lin(i, j), musp_RelErr_lin(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "lin");
    end
end

subaxis(8, 2, [9, 10]);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rd_all(2, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rd_all(2, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*80);
title(['(g) ' DSnms{setInd} ' Arrangement: Detector B Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);



subaxis(8, 2, 11);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{B}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{B}$ (mm)', ...
    'Interpreter','latex');
title('(h) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(8, 2, 12);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{B}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(i) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');






%######################--- 2 ---###########################################
% Move Det 2
cl_mua=[-1, 1]*10;
cl_musp=[-1, 1]*10;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=(-10:0.1:10);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
mua_RelErr_lin=NaN(length(dx), length(dy));
musp_RelErr_lin=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)

        rsPen_act(2, :)=rsPen_ass(2, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
        [mua_RelErr_lin(i, j), musp_RelErr_lin(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "lin");
    end
end

subaxis(8, 2, [13, 14]);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rs_all(2, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rs_all(2, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*80);
title(['(j) ' DSnms{setInd} ' Arrangement: Source 2 Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);



subaxis(8, 2, 15);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{2}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{2}$ (mm)', ...
    'Interpreter','latex');
title('(k) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(8, 2, 16);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{2}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(l) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');






svFig_params.FigSz_cm=[13.5, 22];
svFig_params.doPDF=false;
svFig_params.doPDFraster=true;
saveFig; clear svFig_params;










%% Fig 400 - TRAP-1/A-iter

setInd=3;

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

cols=turbo(100);
cols(1, :)=1;
cols(end, :)=1;

h=figure(400); clf;
h.Name=baseFigName;



%######################--- 1 ---###########################################
% Move Src 1
cl_mua=[-1, 1]*20;
cl_musp=[-1, 1]*50;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=-10:0.1:10;
dy=-10:0.1:10;

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
    end
end



subaxis(4, 2, [1, 2], ...
    'sv', 0.04, 'sh', 0, ...
    'mt', 0.01, 'mb', 0.07, 'mr', 0.1, 'ml', 0.07);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rs_all(1, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rs_all(1, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*80);
pos=get(gca, 'Position');
legend('Detectors', 'Sources', ...
    'Orientation','vertical', 'Location','northeast', ...
    'Color', [0.75, 0.75, 0.75], ...
    'EdgeColor',[0.25, 0.25, 0.25]);
set(gca, 'Position', pos);
title(['(a) ' DSnms{setInd} ' Arrangement: Source 1 Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);



subaxis(4, 2, 3);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(b) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(4, 2, 4);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(c) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');











%######################--- A ---###########################################
% Move Det A
cl_mua=[-1, 1]*50;
cl_musp=[-1, 1]*100;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=(-11.6:0.1:8.4);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
mua_RelErr_lin=NaN(length(dx), length(dy));
musp_RelErr_lin=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(1, :)=rd_ass(1, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
        [mua_RelErr_lin(i, j), musp_RelErr_lin(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "lin");
    end
end

subaxis(4, 2, [5, 6]);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rd_all(1, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rd_all(1, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*80);
title(['(d) ' DSnms{setInd} ' Arrangement: Detector A Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);



subaxis(4, 2, 7);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{A}$ (mm)', ...
    'Interpreter','latex');
title('(e) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(4, 2, 8);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(f) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');




svFig_params.FigSz_cm=[13.5, 16];
svFig_params.doPDF=false;
svFig_params.doPDFraster=true;
saveFig; clear svFig_params;








%% Fig 500 - DRCT-1-iter

setInd=4;

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

cols=turbo(100);
cols(1, :)=1;
cols(end, :)=1;

h=figure(500); clf;
h.Name=baseFigName;



%######################--- 1 ---###########################################
% Move Src 1
cl_mua=[-1, 1]*50;
cl_musp=[-1, 1]*100;

cbTicks_mua=cl_mua(1):cl_mua(2)/2:cl_mua(2);
cbTicks_musp=cl_musp(1):cl_musp(2)/2:cl_musp(2);

dx=-10:0.1:10;
dy=-10:0.1:10;

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter=NaN(length(dx), length(dy));
musp_RelErr_iter=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        [mua_RelErr_iter(i, j), musp_RelErr_iter(i, j), ~, ~] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
    end
end



subaxis(2, 2, [1, 2], ...
    'sv', 0, 'sh', 0.15, ...
    'mt', 0.2, 'mb', 0.07, 'mr', 0.15, 'ml', 0.10);
plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');
plot( ...
    rs_all(1, 1, setInd)+[                    dx, ...
                     ones(size(dx))*dx(end), ...
                               dx(end:-1:1), ...
                      ones(size(dx))*dx(1)], ...
     rs_all(1, 2, setInd)+[ ones(size(dy))*dy(1), ...
                                         dy, ...
                     ones(size(dy))*dy(end), ...
                              dy(end:-1:1)], ...
                                        '-k' ...
    );
text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)-3, '1', 'Color','r');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)-3, '2', 'Color','r');
text(rd_all(1, 1, setInd), rd_all(1, 2, setInd)-3, 'A', 'Color','b');
text(rd_all(2, 1, setInd), rd_all(2, 2, setInd)-3, 'B', 'Color','b');
hold off;
axis equal tight off;
xlim([-1, 1]*80);
pos=get(gca, 'Position');
legend('Detectors', 'Sources', ...
    'Orientation','vertical', 'Location','northeast', ...
    'Color', [0.75, 0.75, 0.75], ...
    'EdgeColor',[0.25, 0.25, 0.25]);
set(gca, 'Position', pos);
title(['(a) ' DSnms{setInd} ' Arrangement: Source 1 Displacement']);
% set(gca, 'Position', get(gca, 'Position')+[0, -0.01, 0, 0]);



subaxis(2, 2, 3);
pcolor(dx, dy, 100*mua_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_mua;
set(gca, 'Position', pos);
contour(dx, dy, 100*mua_RelErr_iter.', ...
    cbTicks_mua, '--k');
clim(cl_mua);
axis equal tight;
ylabel(cb, '$(\mu_a^{(\times \rho)}-\mu_a)/\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(b) Iterative Method $\mu_a$ Error', 'Interpreter','latex');



subaxis(2, 2, 4);
pcolor(dx, dy, 100*musp_RelErr_iter.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_musp;
set(gca, 'Position', pos);
contour(dx, dy, 100*musp_RelErr_iter.', ...
    cbTicks_musp, '--k');
clim(cl_musp);
axis equal tight;
ylabel(cb, '$(\mu_s^{\prime,(\times \rho)}-\mu_s^\prime)/\mu_s^\prime$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(c) Iterative Method $\mu_s^\prime$ Error', 'Interpreter','latex');





svFig_params.FigSz_cm=[13.5, 11];
svFig_params.doPDF=false;
svFig_params.doPDFraster=true;
saveFig; clear svFig_params;






%% TABLE - Absolute RMSE-orbit for Multiple Optodes
optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

rOrb=1; %mm

fprintf(['mua \t musp \t optChar \t' ...
    'mua_CV \t musp_CV \n']);
for setInd=1:4
    fprintf('%s \n', DSnms_short{setInd});
    switch DSnms_short{setInd}
        case 'LINR'
            optChars_MultOpt={
                {'S1'; 'D1'};
                {'S1'; 'D2'};
                {'S1'; 'S2'};
                {'D1'; 'D2'};
                {'S1'; 'D1'; 'D2'};
                {'S1'; 'S2'; 'D1'};
                {'S1'; 'S2'; 'D1'; 'D2'};
                };
        case 'ALIN'
            optChars_MultOpt={
                {'S1'; 'D1'};
                {'S1'; 'D2'};
                {'S2'; 'D1'};
                {'S2'; 'D2'};
                {'S1'; 'S2'};
                {'D1'; 'D2'};
                {'S1'; 'D1'; 'D2'};
                {'S2'; 'D1'; 'D2'};
                {'S1'; 'S2'; 'D1'};
                {'S1'; 'S2'; 'D2'};
                {'S1'; 'S2'; 'D1'; 'D2'};
                };
        case 'TRAP'
            optChars_MultOpt={
                {'S1'; 'D1'};
                {'S1'; 'D2'};
                {'S1'; 'S2'};
                {'D1'; 'D2'};
                {'S1'; 'D1'; 'D2'};
                {'S1'; 'S2'; 'D1'};
                {'S1'; 'S2'; 'D1'; 'D2'};
                };
        case 'DRCT'
            optChars_MultOpt={
                {'S1'; 'D1'};
                {'S1'; 'D2'};
                {'S1'; 'S2'};
                {'S1'; 'D1'; 'D2'};
                {'S1'; 'S2'; 'D1'; 'D2'};
                };
        otherwise
            error();
    end
    
    rsPen_ass=rs_all(:, :, setInd);
    rd_ass=rd_all(:, :, setInd);
            
    for optSetInd=1:length(optChars_MultOpt)
        optChars=optChars_MultOpt{optSetInd};
        
        optNm=[];
        for i=1:length(optChars)
            switch optChars{i}
                case 'S1'
                    optNm=[optNm, '1'];
                case 'S2'
                    optNm=[optNm, '2'];
                case 'D1'
                    optNm=[optNm, 'A'];
                case 'D2'
                    optNm=[optNm, 'B'];
                otherwise
                    error;
            end
            if i~=length(optChars) && length(optChars)~=2
                optNm=[optNm, ', '];
            elseif i~=length(optChars) 
                optNm=[optNm, ' '];
            end
            if i==(length(optChars)-1)
                optNm=[optNm, 'and '];
            end
        end

        fprintf('& %s & ', optNm);

        [mua_RMSE_O_iter, musp_RMSE_O_iter, mua0, musp0]=...
            posOrbitRMSE( ...
            optChars, ...
            rsPen_ass, rd_ass, rOrb, ...
            optProp_act, dmua_act, omega, 'iter');
        mua_CV_O_iter=mua_RMSE_O_iter/mua0;
        musp_CV_O_iter=musp_RMSE_O_iter/musp0;
        
        dec3=-(floor(log10(abs(mua_CV_O_iter*100)))-1);
        dec4=-(floor(log10(abs(musp_CV_O_iter*100)))-1);
        formatSpec=[
            '%.' num2str(dec3) 'f & ' ...
            '%.' num2str(dec4) 'f ' ...
            '\\\\ \n'];
        fprintf(formatSpec, ...
            mua_CV_O_iter*100, musp_CV_O_iter*100);
    end
end