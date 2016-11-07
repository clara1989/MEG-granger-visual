%% Plot the positions of the VE for each visual ROI

figure; hold on;
dataV1 = ft_plot_mesh(sourcemodel_virt.pos(1:49,:),'vertexcolor','r'); 
dataV2 = ft_plot_mesh(sourcemodel_virt.pos(50:91,:),'vertexcolor','b'); 
dataV4 = ft_plot_mesh(sourcemodel_virt.pos(139:177,:),'vertexcolor','g');
dataMT = ft_plot_mesh(sourcemodel_virt.pos(178:184,:),'vertexcolor',[0,0,0]);
dataV7 = ft_plot_mesh(sourcemodel_virt.pos(185:190,:),'vertexcolor',[1,1,0]);
dataPIT = ft_plot_mesh(sourcemodel_virt.pos(191:199,:),'vertexcolor',[1,0,1]);
%legend([data2 data4],'V1','V2')
hold off
legend([dataV1 dataV2 dataV4 dataMT dataV7 dataPIT],'V1','V2','V4','MT','V7','PIT');
ft_plot_mesh(sourcespace,'facecolor','w','edgecolor','none','facealpha',0.1); camlight;
set(gcf,'color','w'); view(0, 0)

