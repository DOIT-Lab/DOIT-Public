%% Setup
clear; home;
% Giles Blaney Summer 2020

% G. Blaney, A. Sassaroli, and S. Fantini, “Design of a source-detector 
% array for dual-slope diffuse optical imaging,” Review of Scientific 
% Instruments, Submitted.

% NOTE: This code has not yet been optimized. Currently, its heavy
% dependence on loops causes it to run slow. Future revisions will be
% optimized to remove loops. 

rng(100);

col=lines(7);
cols_sd=lines(4);

N=500;

m=1;
k=50;
c=2*sqrt(k*m)*0.25;
a=1e3;
rhomin=20;
rhosnom=[25; 35];
levNom=10;
dt=0.005;
fix1=false;
repelD=true;
levSpg=true;

S0=[randn, randn;...
    randn, randn]*20;
D0=[randn, randn;...
    randn, randn]*10;
nS=size(S0, 1);
nD=size(D0, 1);

DSprs(:, :, 1, 1)=[1, 1; 1, 2];
DSprs(:, :, 2, 1)=[2, 2; 2, 1];

X0=[S0; D0];
V0=zeros(size(X0));

%% Make Pairs
prs=[];
prsLev=[];
LevConPInds=[];
r0=[];
for i=1:size(DSprs, 4)
    for j=1:size(DSprs, 3)
        prs=[prs; [DSprs(1, :, j, i); DSprs(2, :, j, i)]+[0, nS; 0, nS]];
        r0=[r0; rhosnom];
        
        if DSprs(1, 1, j, i)==DSprs(2, 1, j, i) %com src
            prsLev=[prsLev; [DSprs(1, 2, j, i), DSprs(2, 2, j, i)]+nS];
            
            LevConPInds=[LevConPInds; [-1, 0]+size(prs, 1)];
        elseif DSprs(1, 2, j, i)==DSprs(2, 2, j, i) %com det
            prsLev=[prsLev; [DSprs(1, 1, j, i), DSprs(2, 1, j, i)]];
            
            LevConPInds=[LevConPInds; [-1, 0]+size(prs, 1)];
        end
    end
end

for prInd=1:size(prs, 1)
    cols_SD(prInd, :)=mean(cols_sd(prs(prInd, :), :));
end

for prInd=1:size(LevConPInds, 1)
    cols_SS(prInd, :)=mean(cols_SD(LevConPInds(prInd, :), :));
end

%% Run
V=V0;
X=NaN([size(X0), N]);
X(:, :, 1)=X0;
figNum=4;
lv_all=NaN(N, size(prsLev, 1));
r_all=NaN(N, size(prs, 1));
for n=1:(size(X, 3)-1)
    A=zeros(size(X0));
    for prInd=1:size(prs, 1)
        r=vecnorm(X(prs(prInd, 1), :, n)-X(prs(prInd, 2), :, n));
        rhat=(X(prs(prInd, 1), :, n)-X(prs(prInd, 2), :, n))/(r+1e-6);
        delta=r-r0(prInd);
        r_all(n, prInd)=r;
        
        v1=V(prs(prInd, 1), :);
        v2=V(prs(prInd, 2), :);
        
        for dfInd=1:size(X0, 2)
            A(prs(prInd, 1), dfInd)=...
                A(prs(prInd, 1), dfInd)...
                -k*delta*rhat(dfInd)/m...
                -c*v1(dfInd)/m;
            A(prs(prInd, 2), dfInd)=...
                A(prs(prInd, 2), dfInd)...
                +k*delta*rhat(dfInd)/m...
                -c*v2(dfInd)/m;
        end
    end
    
    if levSpg
        for lvPrInd=1:size(prsLev, 1)
            rhat=(X(prsLev(lvPrInd, 1), :, n)-X(prsLev(lvPrInd, 2), :, n))/...
                vecnorm(X(prsLev(lvPrInd, 1), :, n)-...
                X(prsLev(lvPrInd, 2), :, n));
            
            pr1=prs(LevConPInds(lvPrInd, 1), :);
            pr2=prs(LevConPInds(lvPrInd, 2), :);
            r1=vecnorm(X(pr1(1), :, n)-X(pr1(2), :, n));
            r2=vecnorm(X(pr2(1), :, n)-X(pr2(2), :, n));
            
            lv_all(n, lvPrInd)=abs(r1-r2);
            
            delta=lv_all(n, lvPrInd)-levNom;
            
            for dfInd=1:size(X0, 2)
                A(prsLev(lvPrInd, 1), dfInd)=...
                    A(prsLev(lvPrInd, 1), dfInd)...
                    -k*delta*rhat(dfInd)/m;
                
                A(prsLev(lvPrInd, 2), dfInd)=...
                    A(prsLev(lvPrInd, 2), dfInd)...
                    +k*delta*rhat(dfInd)/m;
            end
        end
    end
    
    if repelD
        for sInd=1:nS
            for dInd=1:nD
                r=vecnorm(X(sInd, :, n)-X(dInd+nS, :, n));
                rhat=(X(sInd, :, n)-X(dInd+nS, :, n))/r;

                if r<=rhomin
                    for dfInd=1:size(X0, 2)
                        A(sInd, dfInd)=...
                            A(sInd, dfInd)...
                            +a*(1/r-1/rhomin)*rhat(dfInd)/m;
                        A(dInd+nS, dfInd)=...
                            A(dInd+nS, dfInd)...
                            -a*(1/r-1/rhomin)*rhat(dfInd)/m;
                    end
                end
            end
        end
    end
    
    V=A*dt+V;
    X(:, :, n+1)=V*dt+X(:, :, n);
    if fix1
        X(1, :, n+1)=X0(1, :);
    end
end

%% Plot 
figure(2); clf;
subplot(2, 2, [1, 2]);
plot(squeeze(X(1, 1, 1:n)), squeeze(X(1, 2, 1:n)), '.',...
    'color', cols_sd(1, :)); hold on;
h1=plot(X(1, 1, n+1), X(1, 2, n+1), 's', 'color', cols_sd(1, :));
plot(squeeze(X(2, 1, 1:n)), squeeze(X(2, 2, 1:n)), '.',...
    'color', cols_sd(2, :));
h2=plot(X(2, 1, n+1), X(2, 2, n+1), 's', 'color', cols_sd(2, :));
plot(squeeze(X(3, 1, 1:n)), squeeze(X(3, 2, 1:n)), '.',...
    'color', cols_sd(3, :));
h3=plot(X(3, 1, n+1), X(3, 2, n+1), 'o', 'color', cols_sd(3, :));
plot(squeeze(X(4, 1, 1:n)), squeeze(X(4, 2, 1:n)), '.',...
    'color', cols_sd(4, :));
h4=plot(X(4, 1, n+1), X(4, 2, n+1), 'o', 'color', cols_sd(4, :));
hold off;
axis equal tight;
xl=xlim;
yl=ylim;
xlim(xl*1.1);
ylim(yl*1.1);
xlabel('x (mm)');
ylabel('y (mm)');
lgd=legend([h1, h2, h3, h4], '1', '2', 'A', 'B',...
    'orientation', 'horizontal', 'location', 'best');
lgd.NumColumns=2;
title('Evolution');

subplot(2, 2, 3);
clear h1;
for i=1:4
    h1(i)=plot(r_all(:, i), 'color', cols_SD(i, :)); hold on;
end
plot([0, size(X, 3)], [25, 25], ':k');
plot([0, size(X, 3)], [35, 35], ':k'); hold off;
xlim([0, size(X, 3)]);
xlabel('i_{\Deltat}');
ylabel('\rho (mm)');
lgd=legend(h1, '1A', '1B', '2B', '2A',...
    'orientation', 'horizontal', 'location', 'best');
lgd.NumColumns=2;
title('Single-Distance (SD)');

subplot(2, 2, 4);
clear h1;
for i=1:2
    h1(i)=plot(lv_all(:, i), 'color', cols_SS(i, :)); hold on;
end
plot([0, size(X, 3)], [1, 1]*levNom, ':k'); hold off;
ax=gca;
set(ax, 'YAxisLocation', 'right');
xlim([0, size(X, 3)]);
xlabel('i_{\Deltat}');
ylabel('\Delta\rho (mm)');
legend(h1, '1AB', '2BA',...
    'orientation', 'vertical', 'location', 'best');
title('Single-Slope (SS)');