function h = plot_avg(t,m,s,col,facecol,linestyle)
% upper = m+s;            % Upper envelope
% lower = m-s;            % Lower envelope
% filled=[upper, fliplr(lower)];  % Y of the polygon
% xpoints = t;       
% xpoints = [xpoints, fliplr(xpoints)];   % X of the polygon
% fillhandle = fill(xpoints, filled, col);  % plot the data
% set(fillhandle,'EdgeColor',[1 1 1],'FaceAlpha',0.3,'EdgeAlpha',0);  % set edge color
% hold on
% plot(x, y, 'color',col, 'linewidth', 2);    
% xlim([min(x) max(x)]);

if nargin < 6
    linestyle = '-';
end

%%
alpha = 0.5;
edgeColor=facecol+(1-facecol)*0.55;

t = t(:);
m = m(:);
s = s(:);
if length(m)==1
    m = repmat(m,[length(t) 1]);
    s = repmat(s,[length(t) 1]);
end

t = t.';
m = m.';
s = s.';

uE = m+s;
lE = m-s;
yP = [lE,fliplr(uE)];
tP = [t,fliplr(t)];

hold on
patch(tP,yP,1,'facecolor',facecol,'edgecolor','none','facealpha',0.2);
if sum(col~=[1 1 1])>0
    if col == yell
        h = plot(t,m,'Color',yell-0.11,'linestyle',linestyle,'linewidth',2.5);
    else
        h = plot(t,m,'Color',col,'linestyle',linestyle,'linewidth',2.5);
    end
end
% plot(t,lE,'Color',edgeColor);
% plot(t,uE,'Color',edgeColor);
end