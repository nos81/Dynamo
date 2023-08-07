%% plot edit

fig_0 = gcf;
axes_0 = axes('Parent',fig_0);
dyn.plot_pop();
set(axes_0,'FontSize',30,'FontName','CMU Serif','Linewidth',3)
ylabel('State probability','FontSize',30,'FontName','CMU Serif')
xlabel('Time [ns]','FontSize',30,'FontName','CMU Serif')
legend('|0>','|1>','FontSize',30,'FontName','CMU Serif', 'Location', 'northwest', ...
    'Orientation','horizontal')
% title('Optmized control pulse','FontSize',30,'FontName','CMU Serif')
ylim([-0.55 0.7])

xlim([time_base(1) time_base(end)])
legend boxoff 
% grid on; grid minor;
