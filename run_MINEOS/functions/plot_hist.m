function [ Ncount, cent ] = plot_hist( subplt,x,Nbins )
% Plot histograms
x_mean = mean(x);
x_std = std(x);
subplot(subplt(1),subplt(2),subplt(3));
[Ncount, cent] = hist(x,Nbins);
hist(x,Nbins,'color',[1 0 0]); hold on;
% plot([x_mean, x_mean],[0 max(Ncount)],'-','color',[0.8500 0.3250 0.0980],'linewidth',2);
% plot([x_mean+1*x_std, x_mean+1*x_std],[0 max(Ncount)],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);
% plot([x_mean-1*x_std, x_mean-1*x_std],[0 max(Ncount)],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);
x_95_50 = prctile(x,[2.5 50 97.5]);
plot([x_95_50(1) x_95_50(1)],[0 max(Ncount)],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);
plot([x_95_50(3) x_95_50(3)],[0 max(Ncount)],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);
plot([x_95_50(2) x_95_50(2)],[0 max(Ncount)],'-','color',[0.8500 0.3250 0.0980],'linewidth',3);
grid on; box on; set(gca,'fontsize',16,'linewidth',2);

end

