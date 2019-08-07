time_slot = [6 12 24];
voya_distance = [70 150 300 ];
optimal = 0;
LNBD = 1; % Low-complexity near-optimal algorithm
% algorithm = [optimal LNBD];

% finished 1th performance analysis: different distance for testing load shedding and reduced distances
No_test = 1;
complete_accelerate = 3;
range_accelerate = 1;
varphi = linspace(0.0, 1.0, 11);
voya_distance_test = linspace(130, 200, 8);
No_distance = length(voya_distance_test);

% operation_mode_input
% 0~3 normal mode; 4~7 fault mode
mode_with_all_methods_fault = 7;
maxi_time_slot = 2;
optimal_alg = 0;
LNBD = 1;
only_operation_cost = zeros(maxi_time_slot, length(varphi)+1);

for index_distance = 3:1:3
%         if index_mode == 1||index_mode == 5
%             figure
%         end

    [data_1th(index_distance).optimal_cost_related, data_1th(index_distance).complexity, data_1th(index_distance).LS, data_1th(index_distance).RD_SW]  ...
        = cost_optimization_for_test_benders( time_slot(2), voya_distance_test(index_distance),complete_accelerate, optimal_alg, mode_with_all_methods_fault, No_test, varphi_Pl, varphi_Ppr);
    
    data_1th(1).only_operation_cost(1,index_distance) = data_1th(index_distance).optimal_cost_related(1,end);
    data_1th(1).only_complexity(1,index_distance) = sum(data_1th(index_distance).complexity(3,:));
    data_1th(1).only_LS(1,index_distance) = sum( sum( data_1th(index_distance).LS, 1), 2);
    data_1th(1).only_RD_SW(index_distance,:) = data_1th(index_distance).RD_SW;

end

save('data_1th.mat','data_1th');