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



%% Fig 1000


optProp_act.nin=1.4;
optProp_act.nout=1;
optProp_act.musp=1; %1/mm
optProp_act.mua=0.01; %1/mm
dmua_act=0.0001; %1/mm
omega=2*pi*100e6; %rad/s

h=figure(1000); clf;
h.Name=baseFigName;


subaxis(2, 1, 1);

setInd=1;

dx=(-10:1:10);
dy=(-10:1:10);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

[mua_RelErr_iter0, musp_RelErr_iter0, ...
    dmuaDSI_RelErr0, dmuaDSP_RelErr0] = ...
    DSpos2err(...
        rsPen_act, rd_act, rsPen_ass, rd_ass, ...
        optProp_act, dmua_act, omega, "iter");

mua_RelErr_iter_1=NaN(length(dx), length(dy));
musp_RelErr_iter_1=NaN(length(dx), length(dy));
dmuaDSI_RelErr_1=NaN(length(dx), length(dy));
dmuaDSP_RelErr_1=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];
        
        if rsPen_act(1, 1)>0
            continue;
        end
        
        [mua_RelErr_iter_1(i, j), musp_RelErr_iter_1(i, j), ...
            dmuaDSI_RelErr_1(i, j), dmuaDSP_RelErr_1(i, j)] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
    end
end

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter_A=NaN(length(dx), length(dy));
musp_RelErr_iter_A=NaN(length(dx), length(dy));
dmuaDSI_RelErr_A=NaN(length(dx), length(dy));
dmuaDSP_RelErr_A=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(1, :)=rd_ass(1, :)+[dx(i), dy(j), 0];

        if rd_act(1, 1)>0
            continue;
        end
        
        [mua_RelErr_iter_A(i, j), musp_RelErr_iter_A(i, j), ...
            dmuaDSI_RelErr_A(i, j), dmuaDSP_RelErr_A(i, j)] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
    end
end

plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');

text(rs_all(1, 1, setInd)+3, rs_all(1, 2, setInd), '1', 'Color','r');
text(rd_all(1, 1, setInd)-3, rd_all(1, 2, setInd), 'A', 'Color','b');
text(rs_all(2, 1, setInd)-3, rs_all(2, 2, setInd), '2', 'Color','r');
text(rd_all(2, 1, setInd)+3, rd_all(2, 2, setInd), 'B', 'Color','b');

contour(dx+rsPen_ass(1, 1), dy+rsPen_ass(1, 2), 100*mua_RelErr_iter_1.', ...
    [100*mua_RelErr_iter0, 1e6], '-r');
contour(dx+rsPen_ass(1, 1), dy+rsPen_ass(1, 2), 100*musp_RelErr_iter_1.', ...
    [100*musp_RelErr_iter0, 1e6], ':r');
contour(dx+rsPen_ass(1, 1), dy+rsPen_ass(1, 2), 100*dmuaDSI_RelErr_1.', ...
    [100*dmuaDSI_RelErr0, 1e6], '--r');
contour(dx+rsPen_ass(1, 1), dy+rsPen_ass(1, 2), 100*dmuaDSP_RelErr_1.', ...
    [100*dmuaDSP_RelErr0, 1e6], '-.r');

contour(dx+rd_ass(1, 1), dy+rd_ass(1, 2), 100*mua_RelErr_iter_A.', ...
    [100*mua_RelErr_iter0, 1e6], '-b');
contour(dx+rd_ass(1, 1), dy+rd_ass(1, 2), 100*musp_RelErr_iter_A.', ...
    [100*musp_RelErr_iter0, 1e6], ':b');
contour(dx+rd_ass(1, 1), dy+rd_ass(1, 2), 100*dmuaDSI_RelErr_A.', ...
    [100*dmuaDSI_RelErr0, 1e6], '--b');
contour(dx+rd_ass(1, 1), dy+rd_ass(1, 2), 100*dmuaDSP_RelErr_A.', ...
    [100*dmuaDSP_RelErr0, 1e6], '-.b');


contour(-(dx+rsPen_ass(1, 1)), dy+rsPen_ass(1, 2), 100*mua_RelErr_iter_1.', ...
    [100*mua_RelErr_iter0, 1e6], '-r');
contour(-(dx+rsPen_ass(1, 1)), dy+rsPen_ass(1, 2), 100*musp_RelErr_iter_1.', ...
    [100*musp_RelErr_iter0, 1e6], ':r');
contour(-(dx+rsPen_ass(1, 1)), dy+rsPen_ass(1, 2), 100*dmuaDSI_RelErr_1.', ...
    [100*dmuaDSI_RelErr0, 1e6], '--r');
contour(-(dx+rsPen_ass(1, 1)), dy+rsPen_ass(1, 2), 100*dmuaDSP_RelErr_1.', ...
    [100*dmuaDSP_RelErr0, 1e6], '-.r');

contour(-(dx+rd_ass(1, 1)), dy+rd_ass(1, 2), 100*mua_RelErr_iter_A.', ...
    [100*mua_RelErr_iter0, 1e6], '-b');
contour(-(dx+rd_ass(1, 1)), dy+rd_ass(1, 2), 100*musp_RelErr_iter_A.', ...
    [100*musp_RelErr_iter0, 1e6], ':b');
contour(-(dx+rd_ass(1, 1)), dy+rd_ass(1, 2), 100*dmuaDSI_RelErr_A.', ...
    [100*dmuaDSI_RelErr0, 1e6], '--b');
contour(-(dx+rd_ass(1, 1)), dy+rd_ass(1, 2), 100*dmuaDSP_RelErr_A.', ...
    [100*dmuaDSP_RelErr0, 1e6], '-.b');

h1=plot(NaN,NaN,'-k');
h2=plot(NaN,NaN,':k');
h3=plot(NaN,NaN,'--k');
h4=plot(NaN,NaN,'-.k');

hold off;
axis equal tight;
xlim([-35, 35]);
ylim([-5, 5]);
xlabel('$x$ (mm)', 'Interpreter','latex');
ylabel('$y$ (mm)', 'Interpreter','latex');


















subaxis(2, 1, 2);

setInd=3;

dx=(-50:1:50);
dy=(-50:1:50);

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

[mua_RelErr_iter0, musp_RelErr_iter0, ...
    dmuaDSI_RelErr0, dmuaDSP_RelErr0] = ...
    DSpos2err(...
        rsPen_act, rd_act, rsPen_ass, rd_ass, ...
        optProp_act, dmua_act, omega, "iter");

mua_RelErr_iter_1=NaN(length(dx), length(dy));
musp_RelErr_iter_1=NaN(length(dx), length(dy));
dmuaDSI_RelErr_1=NaN(length(dx), length(dy));
dmuaDSP_RelErr_1=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rsPen_act(1, :)=rsPen_ass(1, :)+[dx(i), dy(j), 0];

        if rsPen_act(1, 1)>0
            continue;
        end
        
        [mua_RelErr_iter_1(i, j), musp_RelErr_iter_1(i, j), ...
            dmuaDSI_RelErr_1(i, j), dmuaDSP_RelErr_1(i, j)] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
    end
end

rsPen_ass=rs_all(:, :, setInd);
rd_ass=rd_all(:, :, setInd);
rsPen_act=rsPen_ass;
rd_act=rd_ass;

mua_RelErr_iter_A=NaN(length(dx), length(dy));
musp_RelErr_iter_A=NaN(length(dx), length(dy));
dmuaDSI_RelErr_A=NaN(length(dx), length(dy));
dmuaDSP_RelErr_A=NaN(length(dx), length(dy));
for i=1:length(dx)
    for j=1:length(dy)
        
        rd_act(1, :)=rd_ass(1, :)+[dx(i), dy(j), 0];

        if rd_act(1, 1)>0
            continue;
        end
        
        [mua_RelErr_iter_A(i, j), musp_RelErr_iter_A(i, j), ...
            dmuaDSI_RelErr_A(i, j), dmuaDSP_RelErr_A(i, j)] = ...
            DSpos2err(...
                rsPen_act, rd_act, rsPen_ass, rd_ass, ...
                optProp_act, dmua_act, omega, "iter");
    end
end

plot(rd_all(:, 1, setInd), rd_all(:, 2, setInd), 'ob'); hold on;
plot(rs_all(:, 1, setInd), rs_all(:, 2, setInd), 'sr');

text(rs_all(1, 1, setInd), rs_all(1, 2, setInd)+3, '1', 'Color','r');
text(rd_all(1, 1, setInd)-3, rd_all(1, 2, setInd), 'A', 'Color','b');
text(rs_all(2, 1, setInd), rs_all(2, 2, setInd)+3, '2', 'Color','r');
text(rd_all(2, 1, setInd)+3, rd_all(2, 2, setInd), 'B', 'Color','b');

contour(dx+rsPen_ass(1, 1), dy+rsPen_ass(1, 2), 100*mua_RelErr_iter_1.', ...
    [100*mua_RelErr_iter0, 1e6], '-r');
contour(dx+rsPen_ass(1, 1), dy+rsPen_ass(1, 2), 100*musp_RelErr_iter_1.', ...
    [100*musp_RelErr_iter0, 1e6], ':r');
contour(dx+rsPen_ass(1, 1), dy+rsPen_ass(1, 2), 100*dmuaDSI_RelErr_1.', ...
    [100*dmuaDSI_RelErr0, 1e6], '--r');
contour(dx+rsPen_ass(1, 1), dy+rsPen_ass(1, 2), 100*dmuaDSP_RelErr_1.', ...
    [100*dmuaDSP_RelErr0, 1e6], '-.r');

contour(dx+rd_ass(1, 1), dy+rd_ass(1, 2), 100*mua_RelErr_iter_A.', ...
    [100*mua_RelErr_iter0, 1e6], '-b');
contour(dx+rd_ass(1, 1), dy+rd_ass(1, 2), 100*musp_RelErr_iter_A.', ...
    [100*musp_RelErr_iter0, 1e6], ':b');
contour(dx+rd_ass(1, 1), dy+rd_ass(1, 2), 100*dmuaDSI_RelErr_A.', ...
    [100*dmuaDSI_RelErr0, 1e6], '--b');
contour(dx+rd_ass(1, 1), dy+rd_ass(1, 2), 100*dmuaDSP_RelErr_A.', ...
    [100*dmuaDSP_RelErr0, 1e6], '-.b');



contour(-(dx+rsPen_ass(1, 1)), dy+rsPen_ass(1, 2), 100*mua_RelErr_iter_1.', ...
    [100*mua_RelErr_iter0, 1e6], '-r');
contour(-(dx+rsPen_ass(1, 1)), dy+rsPen_ass(1, 2), 100*musp_RelErr_iter_1.', ...
    [100*musp_RelErr_iter0, 1e6], ':r');
contour(-(dx+rsPen_ass(1, 1)), dy+rsPen_ass(1, 2), 100*dmuaDSI_RelErr_1.', ...
    [100*dmuaDSI_RelErr0, 1e6], '--r');
contour(-(dx+rsPen_ass(1, 1)), dy+rsPen_ass(1, 2), 100*dmuaDSP_RelErr_1.', ...
    [100*dmuaDSP_RelErr0, 1e6], '-.r');

contour(-(dx+rd_ass(1, 1)), dy+rd_ass(1, 2), 100*mua_RelErr_iter_A.', ...
    [100*mua_RelErr_iter0, 1e6], '-b');
contour(-(dx+rd_ass(1, 1)), dy+rd_ass(1, 2), 100*musp_RelErr_iter_A.', ...
    [100*musp_RelErr_iter0, 1e6], ':b');
contour(-(dx+rd_ass(1, 1)), dy+rd_ass(1, 2), 100*dmuaDSI_RelErr_A.', ...
    [100*dmuaDSI_RelErr0, 1e6], '--b');
contour(-(dx+rd_ass(1, 1)), dy+rd_ass(1, 2), 100*dmuaDSP_RelErr_A.', ...
    [100*dmuaDSP_RelErr0, 1e6], '-.b');

hold off;
axis equal tight;
xlim([-30, 30]);
ylim([-15, 15]);
xlabel('$x$ (mm)', 'Interpreter','latex');
ylabel('$y$ (mm)', 'Interpreter','latex');
title(['(b) ' DSnms{3} ' Arrangement']);

sgtitle('Lines of zero change in:');
subaxis(2, 1, 1);
legend([h1, h2, h3, h4], ...
    '$\mu_a^{(\times \rho)}$', '$\mu_s^{\prime,(\times \rho)}$', ...
    'Dual-Slope Intensity''s $\Delta\mu_a$', ...
    'Dual-Slope Phase''s $\Delta\mu_a$', ...
    'Interpreter','latex', ...
    'Location','northoutside', 'Orientation','horizontal', ...
    'NumColumns',2, 'EdgeColor',[1,1,1]);
title(sprintf([ ...
    'from moving a single source (red) or single detector (blue)\n\n' ...
    '(a) ' DSnms{1} ' Arrangement' ...
    ]));



svFig_params.FigSz_cm=[13.5, 12];
svFig_params.doPDF=true;
svFig_params.doPDFraster=false;
saveFig; clear svFig_params;