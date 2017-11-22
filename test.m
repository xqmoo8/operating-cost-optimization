time_slot = [6 12 24];
voya_distance = [70 150 300 ];

% %% 1th performance analysis: different distance for testing load shedding and reduced distances
% complete_accelerate = 3;
% range_accelerate = 1;
% varphi = linspace(0.0, 1.0, 11);
% voya_distance_test = linspace(130, 200, 8);
% No_distance = length(voya_distance_test);
% 
% % operation_mode_input
% % 0~3 normal mode; 4~7 fault mode
% mode_with_all_methods_fault = 7;
% maxi_time_slot = 2;
% optimal_alg = 0;
% LNBD = 1;
% only_operation_cost_1th = zeros(maxi_time_slot, length(varphi)+1);
% 
% for index_distance = 6:1:No_distance
% %         if index_mode == 1||index_mode == 5
% %             figure
% %         end
% 
%     [data_1th(index_distance).optimal_cost_related, data_1th(index_distance).complexity, data_1th(index_distance).RD, data_1th(index_distance).LS]  ...
%         = cost_optimization_for_test_benders( time_slot(2), voya_distance_test(index_distance),complete_accelerate, optimal_alg, mode_with_all_methods_fault  );
%     only_operation_cost_1th(index_distance,1) = data_1th(index_distance).optimal_cost_related(1,end);
% 
% end
% 
% only_operation_cost_1th(:,3) = (only_operation_cost_1th(:,1) - only_operation_cost_1th(1,1))*100/only_operation_cost_1th(1,1);
% only_operation_cost_1th(:,4) = (only_operation_cost_1th(:,2) - only_operation_cost_1th(2,1))*100/only_operation_cost_1th(2,1);
% 
% save('data_1th.mat','data_1th');


%% 2nd performance analysis: different varphi
complete_accelerate = 3;
range_accelerate = 1;
varphi = linspace(0.0, 1.0, 11);

% operation_mode_input
% 0~3 normal mode; 4~7 fault mode
mode_with_all_methods_fault = 7;
maxi_time_slot = 2;
optimal_alg = 0;
LNBD = 1;
only_operation_cost_2nd = zeros(maxi_time_slot, length(varphi)+1);

for index_time_slot = 2:1:maxi_time_slot
%         if index_mode == 1||index_mode == 5
%             figure
%         end
        
%     [optimal_cost_comparison_2nd(index_time_slot,1).data, complexity_2nd(index_time_slot,1).time_iteration]  ...
%         = cost_optimization_for_test_benders( time_slot(index_time_slot), voya_distance(index_time_slot),complete_accelerate, optimal_alg, mode_with_all_methods_fault  );
%     only_operation_cost_2nd(index_time_slot,1) = optimal_cost_comparison_2nd(index_time_slot,1).data(1,end);
    
    for index_varphi = 1:1:9
        [data_1th(index_varphi).optimal_cost_related, data_1th(index_varphi).complexity, data_1th(index_varphi).RD, data_1th(index_varphi).LS]  ...
            = cost_optimization_for_test_benders( time_slot(index_time_slot), voya_distance(index_time_slot),complete_accelerate, LNBD, mode_with_all_methods_fault, varphi(index_varphi)  );

        only_operation_cost_2nd(index_time_slot,index_varphi+1) = optimal_cost_comparison_2nd(index_time_slot,index_varphi+1).data(1,end);
    end
end

only_operation_cost_2nd(3,:) = (only_operation_cost_2nd(1,:) - only_operation_cost_2nd(1,1))*100/only_operation_cost_2nd(1,1);
only_operation_cost_2nd(4,:) = (only_operation_cost_2nd(2,:) - only_operation_cost_2nd(2,1))*100/only_operation_cost_2nd(2,1);

save('optimal_cost_comparison_2nd.mat','optimal_cost_comparison_2nd');
save('complexity_2nd.mat','complexity_2nd');
save('only_operation_cost_2nd.mat','only_operation_cost_2nd');

%% 3rd performance analysis: different algorithm and adjustment methods
time_slot = [6 12 24];
voya_distance = [70 150 300 ];
accele_constraint = [3 4];

optimal = 0;
LNBD = 1; % Low-complexity near-optimal algorithm
varphi = 0.5;
algorithm = [optimal LNBD];

mode_with_all_methods_normal = 3;
mode_3rd = [1 3 4 5 7 8];
% cost_optimization_for_test_benders( time_slot, voya_distance, accelerate_flag_input, near_opt_optimal_input, operation_mode_input, varphi )
for index_mode = 1:1:6
%     for index_accele = 3:1:4
        [data_3rd(index_mode).optimal_cost_related, data_3rd(index_mode).complexity, data_3rd(index_mode).RD, data_3rd(index_mode).LS]  ...
            = cost_optimization_for_test_benders( time_slot(2), voya_distance(2),accele_constraint(1), optimal, mode_3rd(index_mode)-1, varphi );

        only_operation_cost_3rd(index_mode,1) = data_3rd(index_mode).optimal_cost_related(1,end);
        complexity_3rd(index_mode,1) = sum(data_3rd(index_mode).complexity(3,:));
%     end
end

only_operation_cost_3rd(3,:) = (only_operation_cost_3rd(1,:) - only_operation_cost_3rd(1,1))*100/only_operation_cost_3rd(1,1);
only_operation_cost_3rd(4,:) = (only_operation_cost_3rd(2,:) - only_operation_cost_3rd(2,1))*100/only_operation_cost_3rd(2,1);

save('optimal_cost_comparison_3rd.mat','optimal_cost_comparison_3rd');
save('complexity_3rd.mat','complexity_3rd');
save('only_operation_cost_3rd.mat','only_operation_cost_3rd');

%% 4th performance analysis: different load demand
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
