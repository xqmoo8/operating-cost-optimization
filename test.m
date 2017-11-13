time_slot = [6 12 24];
voya_distance = [70 150 300 ];

%% 1th performance analysis: different varphi
complete_accelerate = 3;
range_accelerate = 1;
varphi = linspace(0.0, 1.0, 11);

% operation_mode_input
% 0~3 normal mode; 4~7 fault mode
mode_with_all_methods_fault = 7;
maxi_time_slot = 1
optimal_alg = 0;
LNBD = 1;
% only_operation_cost_1th = zeros(maxi_time_slot, length(varphi)+1)

for index_time_slot = 1:1:maxi_time_slot
%         if index_mode == 1||index_mode == 5
%             figure
%         end
        
    [optimal_cost_comparison_1th(index_time_slot,1).data, complexity_1th(index_time_slot,1).time_iteration]  ...
        = cost_optimization_for_test_benders( time_slot(index_time_slot), voya_distance(index_time_slot),complete_accelerate, optimal_alg, mode_with_all_methods_fault  );
    only_operation_cost_1th(index_time_slot,1) = optimal_cost_comparison_1th(index_time_slot,1).data(1,end);
    
    for index_varphi = 1:1:9
        [optimal_cost_comparison_1th(index_time_slot,index_varphi+1).data, complexity_1th(index_time_slot,index_varphi+1).time_iteration]  ...
            = cost_optimization_for_test_benders( time_slot(index_time_slot), voya_distance(index_time_slot),complete_accelerate, LNBD, mode_with_all_methods_fault, varphi(index_varphi)  );

        only_operation_cost_1th(index_time_slot,index_varphi+1) = optimal_cost_comparison_1th(index_time_slot,index_varphi+1).data(1,end);
    end
end

only_operation_cost_1th(3,:) = (only_operation_cost_1th(1,:) - only_operation_cost_1th(1,1))*100/only_operation_cost_1th(1,1);
only_operation_cost_1th(4,:) = (only_operation_cost_1th(2,:) - only_operation_cost_1th(2,1))*100/only_operation_cost_1th(2,1);
save('optimal_cost_comparison_1th.mat','optimal_cost_comparison_1th');
save('complexity_1th.mat','complexity_1th');
save('only_operation_cost_1th.mat','only_operation_cost_1th');


%% 2nd performance analysis: different algorithm and adjustment methods
% time_val_comparison(index_time,index_accelerate)
operation_mode_input = linspace(0, 7, 8);
accele_constraint = [3 4];
optimal = 0;
LNBD = 1; % Low-complexity near-optimal algorithm
varphi = 0.5;
algorithm = [optimal LNBD];
% mode = [1 4 5 8];

% mode_1th = 8;
mode_with_all_methods_normal = 3;
mode_2th = [1 3 4 5 7 8];
% cost_optimization_for_test_benders( time_slot, voya_distance, accelerate_flag_input, near_opt_optimal_input, operation_mode_input, varphi )
for index_time_slot = 1:1:1
    for index_mode = 1:1:6
%         if index_mode == 1||index_mode == 5
%             figure
%         end
        for index_accele = 1:1:2
            [optimal_cost_comparison_2nd(index_time_slot,index_mode).data, complexity_2nd(index_time_slot,index_mode).time_iteration]  ...
                = cost_optimization_for_test_benders( time_slot(index_time_slot), voya_distance(index_time_slot),index_accele,0,mode_2th(index_mode)-1 );

            only_operation_cost_2nd(index_time_slot,index_mode) = optimal_cost_comparison_2nd(index_time_slot,index_mode).data(1,end);
        end
    end
end

only_operation_cost_2nd(3,:) = (only_operation_cost_2nd(1,:) - only_operation_cost_2nd(1,1))*100/only_operation_cost_2nd(1,1);
only_operation_cost_2nd(4,:) = (only_operation_cost_2nd(2,:) - only_operation_cost_2nd(2,1))*100/only_operation_cost_2nd(2,1);
save('optimal_cost_comparison_2nd.mat','optimal_cost_comparison_2nd');
save('complexity_2nd.mat','complexity_2nd');
save('only_operation_cost_2nd.mat','only_operation_cost_2nd');

%% 3nd performance analysis: different load demand
operation_mode_input = linspace(0, 7, 8);
accele_constraint = [3 4];
optimal = 0;
LNBD = 1; % Low-complexity near-optimal algorithm
varphi = 0.5;
algorithm = [optimal LNBD];

for index_time_slot = 1:1:2
    figure
    for index_accelerate = 1:1:4
        [optimal_cost_comparison_3rd(index_time_slot,index_accelerate).data] = cost_optimization_for_test_benders( time_slot(index_time_slot), voya_distance(index_time_slot),index_accelerate,0,3 );
        only_operation_cost_3rd(index_time_slot,index_accelerate) = optimal_cost_comparison_3rd(index_time_slot,index_accelerate).data(1,end);
    end
end

only_operation_cost_3rd(3,:) = (only_operation_cost_3rd(1,:) - only_operation_cost_3rd(1,1))*100/only_operation_cost_3rd(1,1);
only_operation_cost_3rd(4,:) = (only_operation_cost_3rd(2,:) - only_operation_cost_3rd(2,1))*100/only_operation_cost_3rd(2,1);
save('optimal_cost_comparison_3rd.mat','optimal_cost_comparison_3rd');
save('only_operation_cost_3rd.mat','only_operation_cost_3rd');
