time_slot = [6 12 24];
voya_distance = [70 150 300 ];

for index_time = 1:1:3
%     for index_voyage = 1:1:3
        for index_accelerate = 1:1:4
            [optimal_cost_comparison(index_time,index_accelerate).data ] = cost_optimization_for_test_benders( time_slot(index_time), voya_distance(index_time), index_accelerate );
        end
%     end
end