%% Setup
clear; home;
svFig=true;

baseFigName='PHOTC5_optPosErr';

%% Types of Sets
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



%% Fig 600 - LINR-1/A
setInd=1;

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

cols=parula(100);
cols(1, :)=1;
cols(end, :)=1;

h=figure(600); clf;
h.Name=baseFigName;

%######################--- 1 ---###########################################
% Move Src 1
cl_dmuaDSI=[-1, 1]*1;
cl_dmuaDSP=[-1, 1]*10;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=(-10:0.1:10);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
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
xlim([-1, 1]*70);
pos=get(gca, 'Position');
legend('Detectors', 'Sources', ...
    'Orientation','vertical', 'Location','northeast', ...
    'Color', [0.75, 0.75, 0.75], ...
    'EdgeColor',[0.25, 0.25, 0.25]);
set(gca, 'Position', pos);
title(['(a) ' DSnms{setInd} ' Arrangement: Source 1 Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);





subaxis(4, 2, 3);
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(b) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(4, 2, 4);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(c) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');






%######################--- A ---###########################################
% Move Det A
cl_dmuaDSI=[-1, 1]*3;
cl_dmuaDSP=[-1, 1]*40;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=(-14:0.1:6);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(1, :)=rd_ass(1, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
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
xlim([-1, 1]*70);
title(['(d) ' DSnms{setInd} ' Arrangement: Detector A Displacement']);
set(gca, 'Position', get(gca, 'Position')+[0, -0.02, 0, 0]);





subaxis(4, 2, 7);
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{A}$ (mm)', ...
    'Interpreter','latex');
title('(e) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(4, 2, 8);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(f) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');



svFig_params.FigSz_cm=[13.5, 14];
svFig_params.doPDF=false;
svFig_params.doPDFraster=true;
saveFig; clear svFig_params;










%% Fig 700 - ALIN-1/2/A/B

setInd=2;

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

cols=parula(100);
cols(1, :)=1;
cols(end, :)=1;

h=figure(700); clf;
h.Name=baseFigName;



%######################--- 1 ---###########################################
% Move Src 1
cl_dmuaDSI=[-1, 1]*1;
cl_dmuaDSP=[-1, 1]*20;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=(-10:0.1:10);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
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
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(b) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(8, 2, 4);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(c) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');











%######################--- A ---###########################################
% Move Det A
cl_dmuaDSI=[-1, 1]*3;
cl_dmuaDSP=[-1, 1]*40;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=(-14:0.1:6);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(1, :)=rd_ass(1, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
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
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{A}$ (mm)', ...
    'Interpreter','latex');
title('(e) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(8, 2, 8);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(f) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');






%######################--- B ---###########################################
% Move Det B
cl_dmuaDSI=[-1, 1]*3;
cl_dmuaDSP=[-1, 1]*30;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=(-6:0.1:14);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(2, :)=rd_ass(2, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
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
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{B}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{B}$ (mm)', ...
    'Interpreter','latex');
title('(h) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(8, 2, 12);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{B}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(i) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');






%######################--- 2 ---###########################################
% Move Det 2
cl_dmuaDSI=[-1, 1]*0.5;
cl_dmuaDSP=[-1, 1]*10;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=(-10:0.1:10);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)

        rsPen_act(2, :)=rsPen_ass(2, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
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
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{2}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{2}$ (mm)', ...
    'Interpreter','latex');
title('(k) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(8, 2, 16);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{2}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(l) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');






svFig_params.FigSz_cm=[13.5, 22];
svFig_params.doPDF=false;
svFig_params.doPDFraster=true;
saveFig; clear svFig_params;










%% Fig 800 - TRAP-1/A

setInd=3;

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

cols=parula(100);
cols(1, :)=1;
cols(end, :)=1;

h=figure(800); clf;
h.Name=baseFigName;



%######################--- 1 ---###########################################
% Move Src 1
cl_dmuaDSI=[-1, 1]*2;
cl_dmuaDSP=[-1, 1]*20;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=-10:0.1:10;
dy=-10:0.1:10;

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
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
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(b) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(4, 2, 4);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(c) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');











%######################--- A ---###########################################
% Move Det A
cl_dmuaDSI=[-1, 1]*5;
cl_dmuaDSP=[-1, 1]*30;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=(-11.6:0.1:8.4);
dy=(-10:0.1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
mua_RelErr_lin=NaN(length(dx), length(dy));
musp_RelErr_lin=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(1, :)=rd_ass(1, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
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
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{A}$ (mm)', ...
    'Interpreter','latex');
title('(e) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(4, 2, 8);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{A}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(f) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');




svFig_params.FigSz_cm=[13.5, 16];
svFig_params.doPDF=false;
svFig_params.doPDFraster=true;
saveFig; clear svFig_params;








%% Fig 900 - DRCT-1

setInd=4;

optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

cols=parula(100);
cols(1, :)=1;
cols(end, :)=1;

h=figure(900); clf;
h.Name=baseFigName;



%######################--- 1 ---###########################################
% Move Src 1
cl_dmuaDSI=[-1, 1]*5;
cl_dmuaDSP=[-1, 1]*30;

cbTicks_dmuaDSI=cl_dmuaDSI(1):cl_dmuaDSI(2)/2:cl_dmuaDSI(2);
cbTicks_dmuaDSP=cl_dmuaDSP(1):cl_dmuaDSP(2)/2:cl_dmuaDSP(2);

dx=-10:0.1:10;
dy=-10:0.1:10;

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

dmuaDSI_RelErr=NaN(length(dx), length(dy));
dmuaDSP_RelErr=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        [~, ~, dmuaDSI_RelErr(i, j), dmuaDSP_RelErr(i, j)] = ...
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
pcolor(dx, dy, 100*dmuaDSI_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSI;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSI_RelErr.', ...
    cbTicks_dmuaDSI, '--k');
clim(cl_dmuaDSI);
axis equal tight;
ylabel(cb, '$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
ylabel('$\Delta y_{1}$ (mm)', ...
    'Interpreter','latex');
title('(b) Dual-Slope Intensity $\Delta\mu_a$ Error', 'Interpreter','latex');



subaxis(2, 2, 4);
pcolor(dx, dy, 100*dmuaDSP_RelErr.'); hold on;
plot(0, 0, 'xk');
shading flat;
colormap(cols);
pos=get(gca, 'Position');
cb=colorbar;
cb.Ticks=cbTicks_dmuaDSP;
set(gca, 'Position', pos);
contour(dx, dy, 100*dmuaDSP_RelErr.', ...
    cbTicks_dmuaDSP, '--k');
clim(cl_dmuaDSP);
axis equal tight;
ylabel(cb, '$$(\Delta\mu_a^{(\times \rho)}-\Delta\mu_a)/\Delta\mu_a$ (\%)', ...
    'Interpreter','latex');
xlabel('$\Delta x_{1}$ (mm)', ...
    'Interpreter','latex');
set(gca, 'YTickLabel', {});
title('(c) Dual-Slope Phase $\Delta\mu_a$ Error', 'Interpreter','latex');





svFig_params.FigSz_cm=[13.5, 11];
svFig_params.doPDF=false;
svFig_params.doPDFraster=true;
saveFig; clear svFig_params;






