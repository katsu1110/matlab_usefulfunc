close all

%% this function is OLD and not yet updated!!!

[figPars, axPars] = setPlotPars;
% figPos = [10 10 21 29.7]; % this is in cm
figPos = [5 5 21 20]; % this is in cm
figure(figPars);
set(gcf, 'position', figPos, 'paperposition', figPos);

xbegin = 5;
ybegin = 14;
sq = 3;
offset_figlab = 1.5;
y_figlab = 2.65;
figspace_x = 4.5;
figspace_y = 1.5;

% addpath(genpath('Z:\Katsuhisa\code\integrated\CubeHelix-master\CubeHelix-master'))
% map = colormap(parula(12));
map = pink ;

% Figure 3 in the pupil paper 1 =========================
% A2: Kiwi's threshold
ax_new = axes(axPars, 'position', [xbegin+1.25*sq+0.2 ybegin-sq-figspace_y sq/2 1.325*sq]);
% open figures
fig = openfig('Z:\Katsuhisa\pupil_project\Figures\Figure4_Learning\raw_figs\threshold_ki.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
hc = ax_old.Children;
hc.MarkerFaceColor = [1 1 1];
len_ses = length(hc.XData);
copyobj(ax_old.Children, ax_new); delete(fig);
xlim([0 80])
ylim([0.5 len_ses+0.5])
set(gca, 'XTick', [0 40 80])
[a1,a2] = offset_axis([0.05 0.025], axPars);
set(a1,'YTick',[])
set(a1,'YColor','w')
% set(gca,'XTick',[])
% set(gca,'XColor','w')

% axes(axPars,'position',[xbegin+2*figspace_x-2-offset_figlab ybegin+sq-1 1 1])
% title('B','fontsize',15)
% axis off

% A1: colormap of pupilROC in Kiwi
ax_new = axes(axPars, 'position', [xbegin ybegin-sq-figspace_y 1.25*sq 1.325*sq]);
% open figures
fig = openfig('Z:\Katsuhisa\pupil_project\Figures\Figure4_Learning\raw_figs\ROCcolor_ki.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
hc = ax_old.Children;
img = hc(end).CData;
seq = hc(end).XData;
ses = hc(end).YData;
delete(fig);
% copyobj(ax_old.Children, ax_new); delete(fig);
imagesc(seq, ses, img)
colormap(map)
caxis([0.3 0.8])
% set(gca,'box','off')
hold on;
plot([0 0], [0.5 max(ses)+0.5], '-k')
hold on;
plot([1500 1500],[0.5 max(ses)+0.5], '-k')
hold on;
% col = [0.25 0.25 0.25];
% a = arrow([-550 ki_exnum],[-320 ki_exnum],'TipAngle',45,'Length',4,'width',1.5,'edgecolor',col,'facecolor',col);
ylabel('sessions')
xlim(1000* [-0.2 1.7])
set(ax_new, 'XTick', [0 750 1500])
ylim([0.5 max(ses)+0.5])
set(ax_new,'fontsize',8)
set(ax_new,'fontname','helvetica')
set(ax_new,'TickDir','out')
set(ax_new,'box','off')
set(ax_new,'YTick',[max(ses)-200 max(ses)-150 max(ses)-100 max(ses)-50],'YTickLabel',{'200','150','100','50'})

c = colorbar('eastoutside');
cpos = c.Position;
cpos(1) = 1.85*cpos(1);
cpos(3) = 0.9*cpos(3);
c.Position = cpos;
c.AxisLocation = 'in';
c.Box = 'off';
c.FontSize = 8;
c.TickDirection = 'out';
c.Ticks = [0.5 0.7];

[a1,a2] = offset_axis([0.05 0.025], axPars);
set(a1,'TickLength',[0.01 0.01])

axes(axPars,'position',[xbegin-offset_figlab ybegin 1 1])
title('A','fontsize',15)
axis off

% B2: mango's threshold
ax_new = axes(axPars, 'position', [xbegin+8.2+1.25*sq ybegin-sq-figspace_y sq/2 1.325*sq]);
% open figures
fig = openfig('Z:\Katsuhisa\pupil_project\Figures\Figure4_Learning\raw_figs\threshold_ma.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
hc = ax_old.Children;
hc.MarkerFaceColor = [1 1 1];
len_ses = length(hc.XData);
copyobj(ax_old.Children, ax_new); delete(fig);
xlim([30 90])
ylim([0.5 len_ses+0.5])
set(gca, 'XTick', [30 90])
[a1,a2] = offset_axis([0.05 0.05], axPars);
set(a1,'YTick',[])
set(a1,'YColor','w')
% set(gca,'XTick',[])
% set(gca,'XColor','w')

% axes(axPars,'position',[xbegin+2*figspace_x-2-offset_figlab ybegin+sq-1 1 1])
% title('B','fontsize',15)
% axis off

% B1: colormap of pupilROC in Mango
ax_new = axes(axPars, 'position', [xbegin+8 ybegin-sq-figspace_y 1.25*sq 1.325*sq]);
% open figures
fig = openfig('Z:\Katsuhisa\pupil_project\Figures\Figure4_Learning\raw_figs\ROCcolor_ma.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
hc = ax_old.Children;
img = hc(end).CData;
seq = hc(end).XData;
ses = hc(end).YData;
delete(fig);
% copyobj(ax_old.Children, ax_new); delete(fig);
imagesc(seq, ses, img)
colormap(map)
caxis([0.3 0.8])
set(gca,'box','off')
hold on;
plot([0 0], [0.5 max(ses)+0.5], '-k')
hold on;
plot([1500 1500],[0.5 max(ses)+0.5], '-k')
hold on;
% col = [0.25 0.25 0.25];
% a = arrow([-550 ma_exnum],[-320 ma_exnum],'TipAngle',45,'Length',4,'width',1.5,'edgecolor',col,'facecolor',col);
% c = colorbar('northoutside');
% cpos = c.Position;
% cpos(4) = 0.5*cpos(4);
% c.Position = cpos;
% c.Box = 'off';
% c.FontSize = 8;
% c.TickDirection = 'out';
% c.Ticks = [0.2 0.5 0.8];
xlim(1000* [-0.2 1.7])
set(ax_new, 'XTick', [0 750 1500])
ylim([0.5 max(ses)+0.5])
set(ax_new,'fontsize',8)
set(ax_new,'fontname','helvetica')
set(ax_new,'TickDir','out')
set(ax_new,'TickLength',[0.02 0.02])
set(ax_new,'box','off')
set(ax_new,'YTick',[max(ses)-70 max(ses)-40 max(ses)-10],'YTickLabel',{'70','40','10'})
[a1,a2] = offset_axis([0.05 0.05], axPars);
set(a1,'TickLength',[0.01 0.01])


% axes(axPars,'position',[xbegin+8.5-offset_figlab ybegin 1 1])
% title('B','fontsize',15)
% axis off


% B1: scatter
ax_new = axes(axPars, 'position', [xbegin ybegin-sq-figspace_y-6 1.25*sq 1.325*sq]);
% open figures
fig = openfig('Z:\Katsuhisa\pupil_project\Figures\Figure4_Learning\raw_figs\threVSpupil_ki.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
xlim([0.4 0.85])
ylim([10 70])
xlabel('aROC')
mes = sprintf('psychophysical threshold');
ylabel(mes)
set(gca, 'XTick', [0.5 0.85])
set(gca, 'YTick', [10 40 70])
[a1,a2] = offset_axis([0.05 0.025], axPars);


% A1: scatter mango
ax_new = axes(axPars, 'position', [xbegin+8 ybegin-sq-figspace_y-6 1.25*sq 1.325*sq]);
% open figures
fig = openfig('Z:\Katsuhisa\pupil_project\Figures\Figure4_Learning\raw_figs\threVSpupil_ma.fig','invisible');
ax_old = findobj(fig, 'type', 'axes');
copyobj(ax_old.Children, ax_new); delete(fig);
xlim([0.42 0.7])
ylim([30 90])
set(gca, 'XTick', [0.5 0.7])
set(gca, 'YTick', [30 60 90])
[a1,a2] = offset_axis([0.05 0.025], axPars);

% xlabel
axes(axPars, 'position', [xbegin-0.075 ybegin-sq-figspace_y-1 sq sq]);
text(0, 0, 'time from stimulus onset (ms)','fontsize',8)
axis off

axes(axPars, 'position', [xbegin+1.25*sq+0.3 ybegin-sq-figspace_y-1 sq/2 sq]);
text(0, 0, '(% signal)','fontsize',8)
axis off

message = sprintf('psychophysical \n     threshold');
axes(axPars, 'position', [xbegin+1.25*sq+0.1 ybegin sq sq]);
text(0, 0, message,'fontsize',8)
axis off

axes(axPars, 'position', [xbegin+1.2 ybegin sq sq]);
text(0, 0, 'Animal A','fontsize',10)
axis off

axes(axPars, 'position', [xbegin+9.2 ybegin sq sq]);
text(0, 0, 'Animal B','fontsize',10)
axis off


% 
% axes(axPars, 'position', [xbegin+0.65 ybegin+3.4 sq sq]);
% text(0, 0, '(easy vs hard trials)','fontsize',8)
% axis off

% axes(axPars, 'position', [xbegin+8 ybegin-sq-figspace_y-1.3 sq sq]);
% text(0, 0, 'time from stimulus onset (ms)','fontsize',8)
% axis off


% axes(axPars, 'position', [xbegin+3*sq+1.8*figspace_x+1 ybegin-sq-figspace_y-1.7 sq sq])
% text(0,0,'(% signal)', 'fontsize',8)
% axis off



axes(axPars, 'position', [xbegin+6.85 ybegin-0.8 sq sq])
text(0,0,'aROC (easy vs hard trials)', 'rotation', 270, 'fontsize',8)
axis off


% save figure
dpi = 300;
figname = 'Z:\Katsuhisa\pupil_project\Figures\Figure4_Learning\formatted_figs\fig_paper';
savefig(strcat(figname,'.fig'))
print(gcf,'-dpdf',strcat(figname, '.pdf'),sprintf('-r%d',dpi))
printeps(1,figname)

