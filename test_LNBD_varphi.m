time_slot = 10;
voya_distance = 150;
optimal = 0;
LNBD = 1; % Low-complexity near-optimal algorithm

% finished 1th performance analysis: different distance for testing load shedding and reduced distances
No_test = 0;
full_accelerate = 2;
accelerate = 1;
varphi_Pl = linspace(0.0, 1.0, 11);
varphi_Ppr = linspace(0.0, 1.0, 11);

% operation_mode_input
% 0~3 normal mode; 4~7 fault mode
fault_all_methods = 7;
maxi_time_slot = 2;

% only_operation_cost = zeros(maxi_time_slot, length(varphi)+1);

total_results = zeros(122,11);
% the 4-th is not important.
for index_iteration_D = 0:1
    for index_accelerate = 1:2
        for index_varphi_Ppr = 1:1:11
            for index_varphi_Pl = 1:1:11
                    para_LNBD(1) = varphi_Pl(index_varphi_Pl); 
                    para_LNBD(2) = varphi_Ppr(index_varphi_Ppr);
                    para_LNBD(3) = index_iteration_D;
                    [suboptimal_cost, final_consumed_time] = cost_optimization_for_test_benders( time_slot,  ...
                        voya_distance, index_accelerate, LNBD, fault_all_methods, No_test, para_LNBD);

                    total_results(index_varphi_Pl + (index_varphi_Ppr-1)*11, 1:4) = suboptimal_cost;
                    total_results(index_varphi_Pl + (index_varphi_Ppr-1)*11, 5:6) = final_consumed_time;
            end
        end

        [optimal_cost, final_consumed_time] = cost_optimization_for_test_benders( time_slot,  ...
            voya_distance, full_accelerate, optimal, fault_all_methods, No_test);

        total_results(122, 1:4) = optimal_cost;
        total_results(122, 5:6) = final_consumed_time;


        total_comparison.cost_LS_RD = total_results;

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
            total_comparison.cost_LS_RD_lite(index, 1:6) = total_comparison.cost_LS_RD(index_parameter_1(index), 1:6);
        end
        total_comparison.cost_LS_RD_lite(12, 1:6) = total_comparison.cost_LS_RD(122, 1:6);

        filename = ['total_comparison_Ac.', num2str(index_accelerate), 'iter_D.', num2str(index_iteration_D),  '.mat'];
        save(filename,'total_comparison');
    end
end
