% for operation_mode =0:0
%     filename = ['off_power_data_',num2str(operation_mode),'.mat'];
%     load(filename);
%     Ppr_total(operation_mode+1,:) = (power_data(5,:)./2.2e-3).^(1/3);
%     Operating_cost_total(1,operation_mode+1) = sum(power_data(7,:))
% end
% Operating_cost_total(2,:) = (Operating_cost_total(1,:) - Operating_cost_total(1))*100/Operating_cost_total(1);

for operation_mode =0:1:7
    filename = ['off_power_data_',num2str(operation_mode),'.mat'];
    load(filename);
    Ppr_total(operation_mode+1,:) = (power_data(5,:)./2.2e-3).^(1/3);
    Operating_cost_total(1,operation_mode+1) = sum(power_data(7,:))
end
Operating_cost_total(2,:) = (Operating_cost_total(1,:) - Operating_cost_total(1))*100/Operating_cost_total(1);

figure
hold on
xlabel('Time (hours)','FontSize',9); 
ylabel('Speed (kn)','FontSize',9);
% set(gca,'box','on','Ytick',[]);
xlim([0,25]); 
ylim([12.8,14.1]);
set(gca, 'XTick', [1 3 5 7 9 11 13 15 17 19 21 23]) %设置X坐标轴刻度数据点位置  
set(gca,'XTickLabel',{'1','3','5','7','9','11','13','15','17','19','21','23'}) %设置X坐标轴刻度处显示的字符  

hold on
plot([1:24], Ppr_total(1,:),'-','LineWidth',2);
plot([1:24], Ppr_total(2,:),'*-','LineWidth',2);
plot([1:24], Ppr_total(4,:),'o-','LineWidth',2);
plot([1:24], Ppr_total(6,:),'>-','LineWidth',2);
plot([1:24], Ppr_total(8,:),'s-','LineWidth',2);
hold off
legend('Normal wo PPA & ESMC','Normal w PPA wo ESMC','Normal w PPA w ESMC','Fault w PPA wo ESMC','Fault w PPA w ESMC');