function [fig_axes] = Plot_sphere()
%% Bloch sphere
fig_axes = figure;
t = -1:0.0001:1;
st = 0*t;

view(145,25);%(azimuthal,elevation)
hold on;
text(0,-0.1,1.3,'\bf{|1\rangle}','Color','black','FontSize',20,'FontName','CMU Serif');
text(0,-0.1,-1.3,'\bf{|0\rangle}','Color','black','FontSize',20,'FontName','CMU Serif');
%
% text(0,-0.1,1.1,'$\hat\mathbf{{z}}$','Color','black','FontSize',25,'Interpreter','latex'); %z-axis label



q01 = quiver3(0,0,0,0,0,1,'Color','blue',...
    'LineStyle','-',...
    'LineWidth',1,...
    'AutoScale','off');
q01.ShowArrowHead = 'off';

q02 = quiver3(0,0,0,0,0,-1,'Color','blue',...
    'LineStyle','-',...
    'LineWidth',1,...
    'AutoScale','off');
q02.ShowArrowHead = 'off';


q1 = quiver3(0,0,0,1,0,0,'Color','blue',...
    'LineStyle','-',...
    'LineWidth',1,...
    'AutoScale','off');
q1.ShowArrowHead = 'off';

q2 = quiver3(0,0,0,-1,0,0,'Color','blue',...
    'LineStyle','-',...
    'LineWidth',1,...
    'AutoScale','off');
q2.ShowArrowHead = 'off';
text(1.3,0,0,'$\bf{\hat{x}}$','Color','black','FontSize',20,'Interpreter','latex','FontName','CMU Serif');

q3 = quiver3(0,0,0,0,-1,0,'Color','blue',...
    'LineStyle','-',...
    'LineWidth',1,...
    'AutoScale','off');
q3.ShowArrowHead = 'off';

q4 = quiver3(0,0,0,0,1,0,'Color','blue',...
    'LineStyle','-',...
    'LineWidth',1,...
    'AutoScale','off');
q4.ShowArrowHead = 'off';
text(0,1.3,0,'$\bf{\hat{y}}$','Color','black','FontSize',20,'Interpreter','latex','FontName','CMU Serif');



%draw the unit circles
py = sqrt( 1 - (t.*t));
ny = -sqrt( 1 - (t.*t));
px = sqrt( 1 - (t.*t));
nx = -sqrt( 1 - (t.*t));
pz = sqrt( 1 - (t.*t));
nz = -sqrt( 1 - (t.*t));

plot3(st,py,t,'Color','[0,0,0]','LineWidth',1.2,'LineStyle','-');
plot3(st,ny,t,'Color','[0,0,0]','LineWidth',1.2,'LineStyle','-');
plot3(px,st,t,'Color','[0,0,0]','LineWidth',1.2,'LineStyle','-');
plot3(nx,st,t,'Color','[0,0,0]','LineWidth',1.2,'LineStyle','-');
plot3(px,t,st,'Color','[0,0,0]','LineWidth',1.2,'LineStyle','-');
plot3(nx,t,st,'Color','[0,0,0]','LineWidth',1.2,'LineStyle','-');
%
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'ZTick',[])
set(gca,'Visible','off')

[a, b, c] = sphere(100);
sp = surfl(a, b, c);
set(sp, 'FaceAlpha', 0.09)
shading flat
end