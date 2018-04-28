time_slot = 10;
voya_distance = 150;
optimal = 0;
LNBD = 1; % Low-complexity near-optimal algorithm

% finished 1th performance analysis: different distance for testing load shedding and reduced distances
No_test = 0;
complete_accelerate = 3;
varphi_Pl = linspace(0.0, 1.0, 11);
varphi_Ppr = linspace(0.0, 1.0, 11);

% operation_mode_input
% 0~3 normal mode; 4~7 fault mode
fault_all_methods = 7;
maxi_time_slot = 2;

% only_operation_cost = zeros(maxi_time_slot, length(varphi)+1);

for index_varphi_Pl = 1:1:11
    for index_varphi_Ppr = 1:1:11
        for index_algorithm = 0:1
            [optimal_cost ] = cost_optimization_for_test_benders( time_slot,  ...
                voya_distance, complete_accelerate, index_algorithm, fault_all_methods, No_test, varphi_Pl(index_varphi_Pl), varphi_Ppr(index_varphi_Ppr));

            data_0th(index_varphi_Pl).compare_cost(index_varphi_Pl, index_algorithm+1) = optimal_cost(2, end);
        end
    end
end

save('data_0th.mat','data_0th');