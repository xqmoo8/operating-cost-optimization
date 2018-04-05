% for index_distance = 1:6
%     cost_comparison(index_distance) = data_1th(index_distance).optimal_cost_related(2,end);
%     RD(index_distance) = roundn(data_1th(index_distance).RD, -2);
%     LS(index_distance) = roundn(sum(sum(data_1th(index_distance).LS,1), 2), -2);
%     No_sw(index_distance) = roundn(sum(sum(data_1th(index_distance).No_sw,1), 2), -2);
% end
No_test = 3;
time_slot = [6 12 24];
voya_distance = [70 150 180 ];
accele_constraint = [3 4];

varphi = 0.5;

mode_with_all_methods_normal = 3;
mode_3rd = [1 3 4 5 7 8];

for index_mode = 1:1:6
    filename = ['optimal_alg_data_mode_',num2str(mode_3rd(index_mode)-1),'_D_',num2str(voya_distance(3)),'_No.',num2str(No_test),'.mat'];  
    load(filename);
    Ppr(index_mode,:) = data(5,:);
end