time_slot = [6 12 24];
voya_distance = [70 150 300 ];
% time_val_comparison(index_time,index_accelerate)

for index_time = 1:1:2
    for index_mode = 1:1:8
        if index_mode == 1||index_mode == 5
            figure
        end

        [optimal_cost_comparison_modes(index_time,index_mode).data] = cost_optimization_for_test_benders( time_slot(index_time), voya_distance(index_time),3,0,index_mode-1 );
        only_operation_cost_modes(index_time,index_mode) = optimal_cost_comparison_modes(index_time,index_mode).data(1,end);
    end
end

only_operation_cost_modes(3,:) = (only_operation_cost_modes(1,:) - only_operation_cost_modes(1,1))*100/only_operation_cost_modes(1,1);
only_operation_cost_modes(4,:) = (only_operation_cost_modes(2,:) - only_operation_cost_modes(2,1))*100/only_operation_cost_modes(2,1);
save('optimal_cost_comparison_modes.mat','optimal_cost_comparison_modes');
save('only_operation_cost_modes.mat','only_operation_cost_modes');


for index_time = 1:1:2
    figure
    for index_accelerate = 1:1:4
        [optimal_cost_comparison_accelerate(index_time,index_accelerate).data] = cost_optimization_for_test_benders( time_slot(index_time), voya_distance(index_time),index_accelerate,0,3 );
        only_operation_cost_accelerate(index_time,index_accelerate) = optimal_cost_comparison_accelerate(index_time,index_accelerate).data(1,end);
    end
end

only_operation_cost_accelerate(3,:) = (only_operation_cost_accelerate(1,:) - only_operation_cost_accelerate(1,1))*100/only_operation_cost_accelerate(1,1);
only_operation_cost_accelerate(4,:) = (only_operation_cost_accelerate(2,:) - only_operation_cost_accelerate(2,1))*100/only_operation_cost_accelerate(2,1);
save('optimal_cost_comparison_accelerate.mat','optimal_cost_comparison_accelerate');
save('only_operation_cost_accelerate.mat','only_operation_cost_accelerate');

B £¨£©