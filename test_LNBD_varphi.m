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

total_cost_comparison = zeros(12,11);

for index_varphi_Ppr = 1:1:11
    for index_varphi_Pl = 1:1:11
            [suboptimal_cost, cost_for_comparison] = cost_optimization_for_test_benders( time_slot,  ...
                voya_distance, complete_accelerate, LNBD, fault_all_methods, No_test, varphi_Pl(index_varphi_Pl), varphi_Ppr(index_varphi_Ppr));

            total_cost_comparison(index_varphi_Pl + (index_varphi_Ppr-1)*11, 1) = suboptimal_cost(2, end);
            total_cost_comparison(index_varphi_Pl + (index_varphi_Ppr-1)*11, 2:5) = cost_for_comparison;
    end
end

[optimal_cost, cost_for_comparison] = cost_optimization_for_test_benders( time_slot,  ...
    voya_distance, complete_accelerate, optimal, fault_all_methods, No_test);

total_cost_comparison(122, 1) = optimal_cost(2, end);
 total_cost_comparison(122, 2:5) = cost_for_comparison;

save('total_cost_comparison.mat','total_cost_comparison');

total_comparison.three_cost = total_cost_comparison;
total_comparison.cost_LS_RD(:, 1) = total_cost_comparison(:, 2);
total_comparison.cost_LS_RD(:, 2) = total_cost_comparison(:, 3) / OPTIONS.Penalty_L;
total_comparison.cost_LS_RD(:, 3) = total_cost_comparison(:, 4) / (1.02 * OPTIONS.Penalty_D);

index = 0;
for index_varphi_Ppr = 1:1:11
    for index_varphi_Pl = 1:1:11
            if varphi_Pl(index_varphi_Pl) + varphi_Ppr(index_varphi_Ppr) == 1
                index = index+1;
                index_parameter_1(index) = index_varphi_Pl + (index_varphi_Ppr-1)*11;
            end
    end
end

for index = 1:length(index_parameter_1)
    total_comparison.three_cost_lite(index, 1:5) = total_comparison.three_cost(index_parameter_1(index), 1:5);
    total_comparison.cost_LS_RD_lite(index, 1:3) = total_comparison.cost_LS_RD(index_parameter_1(index), 1:3);
end

