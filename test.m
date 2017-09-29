time_slot = [6 12 24];
voya_distance = [70 150 300 ];
% time_val_comparison(index_time,index_accelerate)

for index_time = 1:1:2
%     for index_voyage = 1:1:3
%         for index_accelerate = 1:1:4
            for index_mode = 1:1:8
                [optimal_cost_comparison(index_time,index_mode).data] = cost_optimization_for_test_benders( time_slot(index_time), voya_distance(index_time),3,0,index_mode-1 );
                only_operation_cost(index_time,index_mode) = optimal_cost_comparison(index_time,index_mode).data(1,end);
            end
%         end
%     end
end

only_operation_cost(3,:) = (only_operation_cost(1,:) - only_operation_cost(1,1))*100/only_operation_cost(1,1);
only_operation_cost(4,:) = (only_operation_cost(2,:) - only_operation_cost(2,1))*100/only_operation_cost(2,1);
save('optimal_cost_comparison.mat','optimal_cost_comparison');