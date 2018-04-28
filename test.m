time_slot = [6 12 24];
voya_distance = [70 150 300 ];

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

% %% 2nd performance analysis: different varphi
% No_test = 2;
% complete_accelerate = 3;
% range_accelerate = 1;
% varphi = linspace(0.1, 0.9, 9);
% 
% % operation_mode_input
% % 0~3 normal mode; 4~7 fault mode
% mode_with_all_methods_fault = 7;
% maxi_time_slot = 2;
% optimal_alg = 0;
% LNBD = 1;
% only_operation_cost = zeros(maxi_time_slot, length(varphi)+1);
% 
% % for index_time_slot = 2:1:maxi_time_slot
% %         if index_mode == 1||index_mode == 5
% %             figure
% %         end
% 
% % based optimal data used for comparsion
% % [data_2nd(1).optimal_cost_related, data_2nd(1).complexity, data_2nd(1).LS, data_2nd(1).RD_SW]  ...
% %     = cost_optimization_for_test_benders( time_slot(2), voya_distance(2),complete_accelerate, optimal_alg, mode_with_all_methods_fault, No_test  );
% % data_2nd(1).only_operation_cost(2,1) = data_2nd(1).optimal_cost_related(1,end);
% 
% for index_varphi = 1:1:9
%     [data_2nd(index_varphi+1).optimal_cost_related, data_2nd(index_varphi+1).complexity, data_2nd(index_varphi+1).LS, data_2nd(index_varphi+1).RD_SW]  ...
%         = cost_optimization_for_test_benders( time_slot(2), voya_distance(2),complete_accelerate, LNBD, mode_with_all_methods_fault, No_test, varphi(index_varphi+1)  );
% 
%     data_2nd(1).only_operation_cost(1,index_varphi+1) = data_2nd(index_distance).optimal_cost_related(1,end);
% end
% 
% data_2nd(1).only_operation_cost(3,:) = (data_2nd(1).only_operation_cost(1,:) - data_2nd(1).only_operation_cost(1,1))*100/data_2nd(1).only_operation_cost(1,1);
% data_2nd(1).only_operation_cost(4,:) = (data_2nd(1).only_operation_cost(2,:) - data_2nd(1).only_operation_cost(1,1))*100/data_2nd(1).only_operation_cost(1,1);
% 
% save('data_2nd.mat','data_2nd');

%% finished 3rd performance analysis: different algorithm and adjustment methods
No_test = 3;
time_slot = [6 12 24];
voya_distance = [70 150 180 ];
accele_constraint = [3 4];

varphi = 0.5;

mode_with_all_methods_normal = 3;
mode_3rd = [1 3 4 5 7 8];
% cost_optimization_for_test_benders( time_slot, voya_distance, accelerate_flag_input, near_opt_optimal_input, operation_mode_input, varphi )
for index_mode = 1:1:6
%     for index_accele = 3:1:4
        [data_3rd(index_mode).optimal_cost_related, data_3rd(index_mode).complexity, data_3rd(index_mode).LS, data_3rd(index_mode).RD_SW]  ...
            = cost_optimization_for_test_benders( time_slot(2), voya_distance(3),accele_constraint(1), optimal, mode_3rd(index_mode)-1, No_test );

        data_3rd(1).only_operation_cost(1,index_mode) = data_3rd(index_mode).optimal_cost_related(1,end);
        data_3rd(1).only_complexity(1,index_mode) = sum(data_3rd(index_mode).complexity(3,:));
%     end
end

data_3rd(1).only_operation_cost(2,:) = (data_3rd(1).only_operation_cost(1,:) - data_3rd(1).only_operation_cost(1,1))*100/data_3rd(1).only_operation_cost(1,1);
data_3rd(1).only_complexity(2,:) = (data_3rd(1).only_complexity(1,:) - data_3rd(1).only_complexity(1,1))*100/data_3rd(1).only_complexity(1,1);

save('data_3rd.mat','data_3rd');

%% 4th performance analysis: different load demand
No_test = 4;
operation_mode_input = linspace(0, 7, 8);
accele_constraint = [3 4];
optimal = 0;
LNBD = 1; % Low-complexity near-optimal algorithm
varphi = 0.5;
algorithm = [optimal LNBD];

for index_time_slot = 1:1:2
    figure
    for index_accelerate = 1:1:4
        [optimal_cost_comparison_4th(index_time_slot,index_accelerate).data] = cost_optimization_for_test_benders( time_slot(index_time_slot), voya_distance(index_time_slot),index_accelerate,0,3 );
        only_operation_cost_4th(index_time_slot,index_accelerate) = optimal_cost_comparison_4th(index_time_slot,index_accelerate).data(1,end);
    end
end

only_operation_cost_4th(3,:) = (only_operation_cost_4th(1,:) - only_operation_cost_4th(1,1))*100/only_operation_cost_4th(1,1);
only_operation_cost_4th(4,:) = (only_operation_cost_4th(2,:) - only_operation_cost_4th(2,1))*100/only_operation_cost_4th(2,1);
save('optimal_cost_comparison_4th.mat','optimal_cost_comparison_4th');
save('only_operation_cost_4th.mat','only_operation_cost_4th');
