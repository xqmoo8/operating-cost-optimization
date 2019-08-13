time_slot = 10;
voya_distance = 160;
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
index_iteration_D = 0;

% only_operation_cost = zeros(maxi_time_slot, length(varphi)+1);

mode = [2 3 6 7];
total_results = zeros(122,11);
% the 4-th is not important.
for index_mode = 1:1:3
%     for index_iteration_D = 0:0
%     if (index_mode == 4) || (index_mode == 5)
%         continue;
%     end
    for index_accelerate = 1:2
        for index_varphi_Ppr = 3:1:7
            for index_varphi_Pl = 1:1:11
                para_LNBD(1) = varphi_Pl(index_varphi_Pl); 
                para_LNBD(2) = varphi_Ppr(index_varphi_Ppr);
                para_LNBD(3) = 0;
                if para_LNBD(1)+para_LNBD(2) == 1
                    [suboptimal_cost, final_consumed_time, dual_gap, reduced_distance,infeasible_flag] = cost_optimization_for_test_benders( time_slot,  ...
                        voya_distance, index_accelerate, LNBD, mode(index_mode), No_test, para_LNBD);

                                        
                    if infeasible_flag == 0
                        total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + index_varphi_Ppr -2, 1:4) = suboptimal_cost;
                        total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + index_varphi_Ppr -2, 5:6) = final_consumed_time;
                        total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + index_varphi_Ppr -2, 7) = dual_gap;
                        total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + index_varphi_Ppr -2, 8) = reduced_distance;
                    else
                        total_results(index_varphi_Pl + (index_varphi_Ppr-1)*11, 1) = inf;
                    end
                end
            end
        end

        [optimal_cost, final_consumed_time, dual_gap, reduced_distance, infeasible_flag] = cost_optimization_for_test_benders( time_slot,  ...
            voya_distance, index_accelerate, optimal, mode(index_mode), No_test);
                    
        if infeasible_flag == 0
            total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + 6, 1:4) = optimal_cost;
            total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + 6, 5:6) = final_consumed_time;
            total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + 6, 7) = dual_gap;
            total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + 6, 8) = reduced_distance;
        else
            total_results((index_mode - 1)*12 + (index_accelerate-1)*6 + 6, 1:4) = inf;
        end

    end
end

filename = ['comparison_multi_mode_2&3_D.', num2str(voya_distance), '_T.', num2str(time_slot),  '.mat'];
save(filename,'total_results');
