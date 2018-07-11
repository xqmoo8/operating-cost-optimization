%% the main struture of benders decomposition algorithm
function [cost_for_comparison, final_consumed_time, dual_gap, reduced_distance, infeasible_flag_out] ...
         = cost_optimization_for_test_benders( time_slot, voya_distance, accelerate_flag_input, near_opt_optimal_input, operation_mode_input, No_test_in, para_LNBD)
global OPTIONS total_sub total_P total_dual accelerate_flag upper_of_lowerbound with_iteration_D infeasible_flag
%% the adjusted parameters 
dbstop if error

case_numb = 2;

% The number of operation time slot
if ~exist('time_slot', 'var')
    OPTIONS.N_t = 10;
else
    OPTIONS.N_t = time_slot;
end

% The voyage distance
if ~exist('voya_distance', 'var')
    OPTIONS.Distance = 140;
else
    OPTIONS.Distance = voya_distance;
end

% Accelerated options
% 1 no accelerated constraints
% 2 power range and lower bound constraints
% 3 only power range
% 4 only lower bound
if ~exist('accelerate_flag_input', 'var')
    accelerate_flag = 2;
else
    accelerate_flag = accelerate_flag_input;
end

% Algorithm selection
% 0 optimal algorithm
% 1 LNBD
if ~exist('near_opt_optimal_input', 'var')
    near_opt_optimal = 0;
else
    near_opt_optimal = near_opt_optimal_input;
end

% Parameters of generator based on former works
% detail information is given in function initial_parameters().
if ~exist('parameter_test_input', 'var')
    parameter_test = 3;
else
    parameter_test = parameter_test_input;
end

% normal operation mode:
% 0 (Fault wo PPA ESMC)
% 1 (Fault wo PPA w ESMC)
% 2 (Fault w PPA w ESMC)
% 3 (Fault w PPA w ESMC w load shedding & reconfiguration)
% adjustment method 3 is not used in normal mode

% semi-island operation mode:
% 4 (Fault wo PPA ESMC)
% 5 (Fault wo PPA w ESMC)
% 6 (Fault w PPA w ESMC)
% 7 (Fault w PPA w ESMC w load shedding & reconfiguration)

% island operation mode:
% 8 (Fault wo PPA ESMC)
% 9 (Fault wo PPA w ESMC)
% 10 (Fault w PPA w ESMC)
% 11 (Fault w PPA w ESMC w load shedding & reconfiguration)
if ~exist('operation_mode_input', 'var')
    operation_mode = 7;
else
    operation_mode = operation_mode_input;
%     if operation_mode == 3
%         operation_mode = 2;
%     end
end

% Test ID for saving results
if ~exist('No_test_in', 'var')
    No_test = 1;
else
    No_test = No_test_in;
end

% Adjustment parameter of load range in LNBD
% Adjustment parameter of propulsion power in LNBD
% with_iteration_D
if ~exist('para_LNBD', 'var')
    OPTIONS.varphi_Pl = 0.7;
    OPTIONS.varphi_Ppr = 0.7;
    with_iteration_D = 0;
else
    OPTIONS.varphi_Pl = para_LNBD(1);
    OPTIONS.varphi_Ppr = para_LNBD(2);
    with_iteration_D = para_LNBD(3);
end

% number of ESM modules
if case_numb == 1
    OPTIONS.N_g = 2;
    OPTIONS.N_e = 4;
elseif case_numb == 2
    OPTIONS.N_g = 3;
    OPTIONS.N_e = 3;
end

% if ~exist('varphi_Ppr', 'var')
% else
% end

infeasible_flag = 0;
infeasible_flag_out = 0;

%% inintial part
% initial the parameters
initial_parameters(parameter_test, operation_mode);
% initial the load demand
load_demand(operation_mode);
% intial the variables and preallocate an arry to the variables
Max_iteration = 800;
delta_g= ones(OPTIONS.N_g, OPTIONS.N_t);
OPTIONS.initial_redundant_sw = [ones(1, OPTIONS.N_t); zeros(1, OPTIONS.N_t); zeros(1, OPTIONS.N_t); ones(1, OPTIONS.N_t)];
redundant_sw = OPTIONS.initial_redundant_sw ;
objval_upperbound = zeros(1, Max_iteration);
objval_lowerbound = zeros(1, Max_iteration);
best_upperbound = zeros(1, Max_iteration);
best_upperbound(1) = inf;
best_lowerbound = zeros(1, Max_iteration);
best_solution_index = 1;
consumed_time = zeros(3, Max_iteration);
RD = zeros(1, Max_iteration);
min_reduced_distance = +inf;
Ppr_best_solution = zeros(1, OPTIONS.N_t);
% Ppr_max_lowerbound = ones(1, OPTIONS.N_t)* sum(OPTIONS.Pg_Max); 

% initial the objective value and dual variable 
if near_opt_optimal == 0
    total_sub(1).sub_optval = 0;
elseif near_opt_optimal == 1
    total_sub(1).sub_optval = zeros(OPTIONS.N_g, OPTIONS.N_t+1);
end
total_sub(1).delta_g = delta_g;
total_sub(1).redundant_sw = redundant_sw;
total_P(1).Pg = zeros(OPTIONS.N_g, OPTIONS.N_t);
total_P(1).Pb = zeros(OPTIONS.N_e, OPTIONS.N_t);
total_P(1).Ppr = zeros(1, OPTIONS.N_t);
% here load shedding is only supported for two island parts
total_P(1).load_shedding = zeros(2, OPTIONS.N_t);
total_P(1).Pc = zeros(OPTIONS.N_g, OPTIONS.N_t);
total_P(1).Pd = zeros(OPTIONS.N_g, OPTIONS.N_t);
total_dual(1).delta_g = zeros(OPTIONS.N_g, OPTIONS.N_t);
total_dual(1).switch = zeros(4, OPTIONS.N_t); 

% generate lowerbound of benders cut
% delta_bound = [ones(1, OPTIONS.N_t); ones(1, OPTIONS.N_t); ];
% delta_bound = delta_g;
% [objval_upperbound_test, benders_cut_lowerbound] = optimization_subproblem( operation_mode, delta_bound, redundant_sw );
% benders_cut_lowerbound = 0 * benders_cut_lowerbound;
benders_cut_lowerbound = 0;
upper_of_lowerbound = 0;

% only semi-island mode need reconfiguration
if (operation_mode <= 3) || ( (operation_mode <=11) && (operation_mode > 7) )
    OPTIONS.Tmin_sw = OPTIONS.N_t;
else % normal mode and island mode do not need reconfiguration
    OPTIONS.Tmin_sw = OPTIONS.N_t;
end

%% start the iteration of the algorithm
for benders_index = 1:800
    % master_problem optimization based on benders decomposition
    % the last input parameter is used to limit the data dimenson for benders cuts
    tic;
    switch near_opt_optimal
        case 0
            [objval_lowerbound(benders_index), delta_g, redundant_sw, benders_cut ] ...
                = optimization_masterproblem(operation_mode, benders_cut_lowerbound, benders_index-1, Ppr_best_solution);
        case 1
            [objval_lowerbound(benders_index), delta_g, redundant_sw, benders_cut ] ...
                = optimization_masterproblem_t(operation_mode, benders_cut_lowerbound, benders_index-1, Ppr_best_solution);
    end
    
    % storage the computational time of master problem in this iteration
    consumed_time(1, benders_index) = toc;
    
    if infeasible_flag == 1
        infeasible_flag_out = infeasible_flag;
        break;
    end
    
    tic;
    % sub_problem optimization based on benders decomposition
    switch near_opt_optimal
        case 0
            [objval_upperbound(benders_index), sub_optval, sub_P, sub_dual, reduced_distance] ...
                = optimization_subproblem(operation_mode, delta_g, redundant_sw);
        case 1
            [objval_upperbound(benders_index), sub_optval, sub_P, sub_dual, reduced_distance] ...
                = optimization_subproblem_t(operation_mode, delta_g, redundant_sw);
%             [objval_upperbound(benders_index), sub_optval, sub_P, sub_dual, reduced_distance] ...
%                 = optimization_subproblem(operation_mode, delta_g, redundant_sw);
    end
    % storage the computational time of subproblem in this iteration
    consumed_time(2, benders_index) = toc;
    
    if infeasible_flag == 1
        infeasible_flag_out = infeasible_flag;
        break;
    end

    % storage and passing of objective value and dual variable 
    total_sub(benders_index).sub_optval = sub_optval;
    total_sub(benders_index).delta_g = delta_g;
    total_sub(benders_index).redundant_sw = redundant_sw;
    total_P(benders_index).Pg = sub_P.Pg;
    total_P(benders_index).Pb = sub_P.Pb;
    total_P(benders_index).Ppr = sub_P.Ppr;
    total_P(benders_index).load_shedding = sub_P.load_shedding;
    total_P(benders_index).Pc = sub_P.Pc;
    total_P(benders_index).Pd = sub_P.Pd;
    total_dual(benders_index).delta_g = sub_dual.delta_g;
    total_dual(benders_index).switch = sub_dual.switch;
    RD(benders_index) = reduced_distance;
    
    if objval_lowerbound(benders_index) > upper_of_lowerbound
        upper_of_lowerbound = objval_lowerbound(benders_index);
    end
%     if benders_cut > benders_cut_lowerbound
%         benders_cut_lowerbound = benders_cut;
%     end

    % storage the computational time in this iteration
    consumed_time(3, benders_index) = sum(consumed_time(1:2, benders_index), 1);
    % print the iteration number and computation time
    if near_opt_optimal == 0
        fprintf('\nGBD: computation time of MP & SP at iteration %d is %5.2f + %6.2f = %6.2f s \n\n', ...
            benders_index, consumed_time(1, benders_index), consumed_time(2, benders_index), consumed_time(3,benders_index) );
    else
        fprintf('\nLNBD: computation time of MP & SP at iteration %d is %5.2f + %6.2f = %6.2f s \n\n', ...
            benders_index, consumed_time(1, benders_index), consumed_time(2, benders_index), consumed_time(3,benders_index) );
    end

    if benders_index >= 2
        % update upperbound and lowerbound 
        if best_upperbound(benders_index-1) > objval_upperbound(benders_index)
            best_upperbound(benders_index) = objval_upperbound(benders_index);
            best_solution_index = benders_index;
        else
            best_upperbound(benders_index) = best_upperbound(benders_index-1);
        end

        if best_lowerbound(benders_index-1) < objval_lowerbound(benders_index)
            best_lowerbound(benders_index) = objval_lowerbound(benders_index);
        else
            best_lowerbound(benders_index) = best_lowerbound(benders_index-1);
        end
        
%         if reduced_distance < min_reduced_distance
%             min_reduced_distance = reduced_distance;
            Ppr_best_solution = total_P(benders_index).Ppr;
%         end
        
        % performance comparison
        error = ( best_upperbound(benders_index) - best_lowerbound(benders_index) );
        dual_gap = 100*error/best_lowerbound(benders_index);
        disp('upperbound, lowerbound, error, dual_gap');
        disp([best_upperbound(benders_index) best_lowerbound(benders_index) error dual_gap]);
        
        % convergence determination
        if error <= 1e-3
            break;
        end
    end

%     if abs(objval_upperbound(benders_index) - 7.3724e+04) < 10
%         disp('debug');
%     end
end

best_solution.delta_g = total_sub(best_solution_index).delta_g;
best_solution.redundant_sw = total_sub(best_solution_index).redundant_sw;
best_solution.Pg = total_P(best_solution_index).Pg;
best_solution.Pb = total_P(best_solution_index).Pb;
best_solution.Ppr = total_P(best_solution_index).Ppr;
best_solution.load_shedding = total_P(best_solution_index).load_shedding;
best_solution.Pc = total_P(best_solution_index).Pc;
best_solution.Pd = total_P(best_solution_index).Pd;
best_solution.RD = RD(best_solution_index);

optimal_cost = [best_lowerbound(1:benders_index); best_upperbound(1:benders_index); ];
final_consumed_time(1, 1) = sum(consumed_time(3, 1:benders_index));
final_consumed_time(1, 2) = benders_index;
cost_for_comparison = save_optimal_cost_information(No_test, near_opt_optimal, optimal_cost, best_solution, operation_mode, consumed_time(:, 1:benders_index));
% cost_for_comparison = save_optimal_cost_information(No_test, near_opt_optimal, optimal_cost, sub_P, delta_g, redundant_sw, reduced_distance, operation_mode, consumed_time(:, 1:benders_index));

plot_result(optimal_cost);
end

%% Save optimal operation status
function [cost_for_comparison] = save_optimal_cost_information(No_test, near_opt_optimal, optimal_cost, best_solution, operation_mode, consumed_time)
global OPTIONS accelerate_flag
operating_cost = OPTIONS.G(1:OPTIONS.N_g, 1).' * power(best_solution.Pg(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
               + OPTIONS.G(1:OPTIONS.N_g, 2).' * best_solution.Pg(1:OPTIONS.N_g, 1:OPTIONS.N_t) ...
               + OPTIONS.G(1:OPTIONS.N_g, 3).' * best_solution.delta_g ...
               + OPTIONS.Xi_E * OPTIONS.E(1:OPTIONS.N_e, 1).' * power(best_solution.Pb(1:OPTIONS.N_e, 1:OPTIONS.N_t), 2) ...
               + OPTIONS.Xi_E * OPTIONS.E(1:OPTIONS.N_e, 2).' * ones(2, OPTIONS.N_t);
Penalty_cost_LS =  OPTIONS.Penalty_L * sum(best_solution.load_shedding(1:2, 1:OPTIONS.N_t),1);
Penalty_cost_D = 1.02 * OPTIONS.Penalty_D * best_solution.RD;
startup_g =  best_solution.delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) - [zeros(OPTIONS.N_g,1) best_solution.delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t-1)];
startup_cost = OPTIONS.C_ss(1:OPTIONS.N_g,1).'* (round(startup_g) > 0);
operating_cost = operating_cost + startup_cost;

% total cost obtained from the master problem
data.cost(1, 1) = optimal_cost(2, end);
% penalty cost of relaxed power Pc and Pd
data.cost(1, 2) = + sum(10 * OPTIONS.Penalty_D * best_solution.Pc(1, 1:OPTIONS.N_t).' ,1) ...
                  + sum(10 * OPTIONS.Penalty_D * best_solution.Pc(2, 1:OPTIONS.N_t).' ,1) ...
                  + sum(10 * OPTIONS.Penalty_D * best_solution.Pd(1, 1:OPTIONS.N_t).' ,1) ...
                  + sum(10 * OPTIONS.Penalty_D * best_solution.Pd(2, 1:OPTIONS.N_t).' ,1) ;
% operating cost calculated based on the results such as Pg, Ppr and Pe
data.cost(2, 1:OPTIONS.N_t) = operating_cost;
% total operating cost of generators and ESMs
data.cost(3, 1) = sum(operating_cost, 2);
data.cost(4, 1) = sum(operating_cost, 2);
% penalty cost of load shedding
data.cost(3, 2) = sum(Penalty_cost_LS, 2);
data.cost(4, 2) = sum(Penalty_cost_LS, 2)/OPTIONS.Penalty_L;
% penalty cost of reduced distance
data.cost(3, 3) = Penalty_cost_D;
data.cost(4, 3) = best_solution.RD;
% total cost of operating cost, two terms of penalty cost
data.cost(3, 4) = sum(operating_cost) + sum(Penalty_cost_LS) + Penalty_cost_D;
data.cost(4, 4) = sum(operating_cost) + sum(Penalty_cost_LS) + Penalty_cost_D;

% data for parameter passing
cost_for_comparison = data.cost(4, 1:4);

total_switching = sum(abs(best_solution.redundant_sw(1,1) - OPTIONS.initial_redundant_sw(1,1)) ) ...
                + sum(abs(best_solution.redundant_sw(2,1) - OPTIONS.initial_redundant_sw(2,1))) ...
                + sum(abs(best_solution.redundant_sw(1,2:OPTIONS.N_t) - best_solution.redundant_sw(1,1:OPTIONS.N_t-1)) ) ...
                + sum(abs(best_solution.redundant_sw(2,2:OPTIONS.N_t) - best_solution.redundant_sw(2,1:OPTIONS.N_t-1)));

data.Pgc(1:OPTIONS.N_g, 1:OPTIONS.N_t) =  best_solution.Pg;
data.Pgc(OPTIONS.N_g+1:2*OPTIONS.N_g, 1:OPTIONS.N_t) =  best_solution.Pc;
data.Pb(1:OPTIONS.N_e, 1:OPTIONS.N_t) =  best_solution.Pb;
data.Ppr(OPTIONS.N_g+OPTIONS.N_e+1, 1:OPTIONS.N_t) =  best_solution.Ppr;
data.load_shedding(1:2, 1:OPTIONS.N_t) =  best_solution.load_shedding;
if operation_mode <= 3
    data.Pl(1, 1:OPTIONS.N_t) = OPTIONS.P_L_TIME_off;
    data.Pl(3, 1:OPTIONS.N_t) =  sum(best_solution.Pg) + sum(best_solution.Pb, 1) - best_solution.Ppr;
elseif operation_mode <= 11
    data.Pl(1, 1:OPTIONS.N_t) = OPTIONS.island1_load;
    data.Pl(2, 1:OPTIONS.N_t) = OPTIONS.island2_load;
    data.Pl(3, 1:OPTIONS.N_t) =  sum(best_solution.Pg) + sum(best_solution.Pb, 1) - best_solution.Ppr;
end
% data.power(10,:) = OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t);

data.status(1:OPTIONS.N_g, 1:OPTIONS.N_t) = best_solution.delta_g;
data.status(OPTIONS.N_g+1:OPTIONS.N_g+2, 1:OPTIONS.N_t) = [best_solution.redundant_sw(1,1:OPTIONS.N_t); best_solution.redundant_sw(2,1:OPTIONS.N_t)];
data.status(OPTIONS.N_g+3, 1) = total_switching;    

data.distance(1, 1:2) = [OPTIONS.Distance best_solution.RD];

data.consumed_time = consumed_time;
data.consumed_time(4, 1) = sum(consumed_time(3, :));
data.consumed_time(4, 2) = length(consumed_time(3, :));
     
Mg = OPTIONS.M(1:OPTIONS.N_g,1).' * power(best_solution.Pg(1:OPTIONS.N_g,1:OPTIONS.N_t).',2).' ...
   + OPTIONS.M(1:OPTIONS.N_g,2).' * power(best_solution.Pg(1:OPTIONS.N_g,1:OPTIONS.N_t).',1).' ...
   + OPTIONS.M(1:OPTIONS.N_g,3).' * best_solution.delta_g(1:OPTIONS.N_g,1:OPTIONS.N_t);
EEOI_t = Mg./(OPTIONS.M_sl*(best_solution.Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3));

for index_t = find(best_solution.Ppr == 0)
    EEOI_t(index_t) = 0;
end
data.EEOI(1,  1:OPTIONS.N_t) = EEOI_t;

if near_opt_optimal ==0
    filename = ['optimal_mode_',num2str(operation_mode),'_D.',num2str(OPTIONS.Distance),...
                '_T.',num2str(OPTIONS.N_t),'_Ac.',num2str(accelerate_flag),'_No.',num2str(No_test),'.mat'];
elseif near_opt_optimal ==1
    filename = ['LNBD_mode_',num2str(operation_mode),'_D.',num2str(OPTIONS.Distance),'_varphi',num2str(OPTIONS.varphi_Pl),...
                '_T.',num2str(OPTIONS.N_t),'_Ac.',num2str(accelerate_flag),'_No.',num2str(No_test),'.mat'];
end
save(filename,'data');

end

%% Plot the simulation result
function plot_result(result)
hold on
plot(result(1,1:end));
plot(result(2,1:end));
hold off
end

%% Load data generation
function load_demand(operation_mode)
global OPTIONS
OPTIONS.P_L = [4 2]; % P_Generater
OPTIONS.P_L_Scale_t =  [0.5 0.6 0.8 0.9 0.82 0.6 0.47 0.42 0.4 0.47 0.49 0.52 0.49 0.47 0.6 0.8 0.82 0.9 0.9 0.7 0.62 0.5 0.33 0.4 0.5 0.4];
OPTIONS.P_vs_t_Scale = [0.5 0.6 0.8 0.9 0.82 0.6 0.47 0.42 0.4 0.47 0.49 0.52 0.49 0.47 0.6 0.8 0.82 0.9 0.9 0.7 0.62 0.5 0.33 0.4 0.5 0.4];
OPTIONS.P_no_t_Scale = [0.5 0.6 0.8 0.9 0.82 0.6 0.47 0.42 0.4 0.47 0.49 0.52 0.49 0.47 0.6 0.8 0.82 0.9 0.9 0.7 0.62 0.5 0.33 0.4 0.5 0.4];

% total load upper bound 3.6 MW
P_no_base = 1.2; % the scale of non-vital loads are 1.2/3.6 = 0.33
P_vs_base = ones(1,OPTIONS.zone).' * (sum(OPTIONS.P_L) - 1.2)/OPTIONS.zone;

if operation_mode <= 3
    P_vs = P_vs_base * OPTIONS.P_vs_t_Scale;
    P_no = P_no_base * OPTIONS.P_no_t_Scale;
    P_total_vs = sum(P_vs);
    
    P_total_time = P_total_vs + P_no;
    P_total_average = sum(P_total_time(1:OPTIONS.N_t))/OPTIONS.N_t;
    
    OPTIONS.P_total_vs = P_total_vs(1:OPTIONS.N_t);
    OPTIONS.P_no = P_no(1:OPTIONS.N_t);
    OPTIONS.P_L_TIME_off = P_total_time(1:OPTIONS.N_t);
    OPTIONS.P_L_TIME_off_avg = P_total_average;
    OPTIONS.Delta_PL =  OPTIONS.P_L_TIME_off - OPTIONS.P_L_TIME_off_avg ;
    
elseif operation_mode <= 7
    % the default fault is at 5-6 bus and 7-8 bus. Thus,the scale factors of
    % vital and semi-vital loads in different island parts are 3/6 and 1/6;
    % The scale factors of non-vital loads are (5+9-6)/12 and 1-(5+9-6)/12
    % (0.75 and 0.25)
    P_vs_island1 = sum(P_vs_base(1:3)) * OPTIONS.P_vs_t_Scale;
    P_vs_island2 = sum(P_vs_base(6)) * OPTIONS.P_vs_t_Scale;
    P_no_island1 = 0.75 * P_no_base * OPTIONS.P_no_t_Scale;
    P_no_island2 = 0.25 * P_no_base * OPTIONS.P_no_t_Scale;
    P_coupled_load = (P_vs_base(4:5)) * OPTIONS.P_vs_t_Scale;
    
    OPTIONS.Coupled_load = P_coupled_load(1:2,1:OPTIONS.N_t);
    OPTIONS.island1_load = P_vs_island1(1:OPTIONS.N_t) + P_no_island1(1:OPTIONS.N_t);
    OPTIONS.island1_load_average = sum(OPTIONS.island1_load)/OPTIONS.N_t;
    OPTIONS.island1_non_load = P_no_island1(1:OPTIONS.N_t);
    % not the real max and min, which can be considered as the minimum value of upperbound
    % , and the maximum value of lowbound.
    OPTIONS.island1_max = P_vs_island1(1:OPTIONS.N_t) ;
    OPTIONS.island1_min = OPTIONS.island1_load(1:OPTIONS.N_t) + sum(P_coupled_load(1:2, 1:OPTIONS.N_t), 1);
    OPTIONS.Delta_PL_island1 = OPTIONS.island1_load - OPTIONS.island1_load_average;
    
    OPTIONS.island2_load = P_vs_island2(1:OPTIONS.N_t) + P_no_island2(1:OPTIONS.N_t);
    OPTIONS.island2_load_average = sum(OPTIONS.island2_load)/OPTIONS.N_t;
    OPTIONS.island2_non_load = P_no_island2(1:OPTIONS.N_t);
    % not the real max, it's the maximum of the supply power of vital and semi loads
    OPTIONS.island2_max = P_vs_island2(1:OPTIONS.N_t) ;
    OPTIONS.island2_min = OPTIONS.island2_load(1:OPTIONS.N_t) + sum(P_coupled_load(1:2, 1:OPTIONS.N_t), 1);
    OPTIONS.Delta_PL_island2 = OPTIONS.island2_load - OPTIONS.island2_load_average;
    
%     P_vs = P_vs_base * OPTIONS.P_vs_t_Scale;
%     P_no = P_no_base * OPTIONS.P_no_t_Scale;
%     P_total_vs = sum(P_vs);
%     
%     P_total_time = P_total_vs + P_no;
%     OPTIONS.P_L_TIME_off = P_total_time(1:OPTIONS.N_t);
elseif operation_mode <= 11
    % the default fault is at 5-6 bus and 11-12 bus. Thus,the scale factors of
    % vital and semi-vital loads in different island parts are 5/6 and 1/6;
    % The scale factors of non-vital loads are (5+5)/12 and 1-(5+5)/12
    % (5/6 and 1/6)
    P_vs_island1 = sum(P_vs_base(1:5)) * OPTIONS.P_vs_t_Scale;
    P_vs_island2 = sum(P_vs_base(6)) * OPTIONS.P_vs_t_Scale;
    P_no_island1 = 5/6 * P_no_base * OPTIONS.P_no_t_Scale;
    P_no_island2 = 1/6 * P_no_base * OPTIONS.P_no_t_Scale;
    
    OPTIONS.island1_load = P_vs_island1(1:OPTIONS.N_t) + P_no_island1(1:OPTIONS.N_t);
    OPTIONS.island1_load_average = sum(OPTIONS.island1_load)/OPTIONS.N_t;
    OPTIONS.island1_non_load = P_no_island1(1:OPTIONS.N_t);
    % not the real max and min, which can be considered as the minimum value of upperbound
    % , and the maximum value of lowbound.
    OPTIONS.island1_max = OPTIONS.island1_load;
    OPTIONS.island1_min = P_vs_island1(1:OPTIONS.N_t);
    OPTIONS.Delta_PL_island1 = OPTIONS.island1_load - OPTIONS.island1_load_average;
    
    OPTIONS.island2_load = P_vs_island2(1:OPTIONS.N_t) + P_no_island2(1:OPTIONS.N_t);
    OPTIONS.island2_load_average = sum(OPTIONS.island2_load)/OPTIONS.N_t;
    OPTIONS.island2_non_load = P_no_island2(1:OPTIONS.N_t);
    % not the real max, it's the maximum of the supply power of vital and semi loads
    OPTIONS.island2_max = OPTIONS.island2_load;
    OPTIONS.island2_min = P_vs_island2(1:OPTIONS.N_t);
    OPTIONS.Delta_PL_island2 = OPTIONS.island2_load - OPTIONS.island2_load_average;
    
%     P_vs = P_vs_base * OPTIONS.P_vs_t_Scale;
%     P_no = P_no_base * OPTIONS.P_no_t_Scale;
%     P_total_vs = sum(P_vs);
%     
%     P_total_time = P_total_vs + P_no;
%     OPTIONS.P_L_TIME_off = P_total_time(1:OPTIONS.N_t);
end

end

%% Inital the parameters for optimization
function [] = initial_parameters(parameter_test, operation_mode)
global OPTIONS
OPTIONS.velocity = [25 0];
% % number of ESM modules
% OPTIONS.N_e = 4;
% % number of generators
% OPTIONS.N_g = 2;
% number of time slots
% OPTIONS.N_t = 12;
% minmum time between startup and shuntdown
OPTIONS.Tmin_g_st = 1;
% the minimum shutdown time must be greater than 1 (not 0)
OPTIONS.Tmin_g_sd = 1;
% minmum time between switch change
OPTIONS.Tmin_sw = 1;
% time of faults happen
OPTIONS.Fault_time = 4;
% the number of zones in shipboard power system
OPTIONS.zone = 6;

OPTIONS.velocity_avg = OPTIONS.Distance/OPTIONS.N_t;
OPTIONS.P_pr_avg = (OPTIONS.velocity_avg).^3*2.2e-3;
OPTIONS.Pg_Max(1,1) = 8;
OPTIONS.Pg_Min(1,1) = 1;
OPTIONS.Pg_Max(1,2) = 4;
OPTIONS.Pg_Min(1,2) = 0.5;
OPTIONS.Pg_Max(1,3) = 4;
OPTIONS.Pg_Min(1,3) = 0.5;

OPTIONS.Xi_E = 1;
OPTIONS.Ppr_Max = (OPTIONS.velocity(1)).^3*2.2e-3;
% in the normal mode, there are two scenarios (4 ESMs in case 1, 3 ESMs in case 2)
% if operation_mode <= 3
if OPTIONS.N_e == 4
    OPTIONS.Pb_Max(1) = 0.5;
    OPTIONS.Pb_Min(1) = -0.5;
    OPTIONS.E_Max = [1 1 1 1];
    OPTIONS.E_Min = [0.2 0.2 0.2 0.2];
    OPTIONS.E(1:OPTIONS.N_e,1) = 1 * ones(OPTIONS.N_e,1);
    OPTIONS.E(1:OPTIONS.N_e,2) = 0.5 * ones(OPTIONS.N_e,1);
%     OPTIONS.E(1:2,1:2) = [1 0.5; 1 0.5];
elseif OPTIONS.N_e == 3
    OPTIONS.Pb_Max(1) = 0.5;
    OPTIONS.Pb_Min(1) = -0.5;
    OPTIONS.E_Max = [1 1 1];
    OPTIONS.E_Min = [0.2 0.2 0.2];
    OPTIONS.E(1:OPTIONS.N_e,1) = 1 * ones(OPTIONS.N_e,1);
    OPTIONS.E(1:OPTIONS.N_e,2) = 0.5 * ones(OPTIONS.N_e,1);
%     OPTIONS.E(1:2,1:2) = [1 0.5; 1 0.5];
end
% % in the operation mode 4~7, 
% elseif operation_mode <= 7 
%     OPTIONS.Pb_Max(1) = 1;
%     OPTIONS.Pb_Min(1) = -1;
%     OPTIONS.E_Max = [3 1];
%     OPTIONS.E_Min = [0.2 0.2];
% elseif operation_mode <= 11
%     OPTIONS.Pb_Max(1) = 1;
%     OPTIONS.Pb_Min(1) = -1;
%     OPTIONS.E_Max = [3 1];
%     OPTIONS.E_Min = [0.2 0.2];
% end
% a_lc and c_lc


% generator function parameters
switch parameter_test
    case 1 % 
        OPTIONS.G(1,1:3) = [10 19 170];
        OPTIONS.G(3,1:3) = [28.5 28 190];
        OPTIONS.G(2,1:3) = [10 19 170];
        
    case 2 % similar minimum SFC (160), large SFC deviation (250 and 150) between lower bound Pmin and efficient point Peff 
        OPTIONS.G(1,1:3) = [12 60 350];
        OPTIONS.G(3,1:3) = [30 12 170];
        OPTIONS.G(2,1:3) = [12 60 350];
        
    % Kanellos: OPMS in IEEE Transactions on Sustainable Energy, Electric power systems research, and Inventions     
    case 3 
        OPTIONS.G(1,1:3) = [5.4 61.5 390];
        OPTIONS.G(2,1:3) = [13.1 12 430];
        OPTIONS.G(3,1:3) = [5.4 63 400];
        OPTIONS.G(4,1:3) = [13.5 10 450];
        
    % Ce Shang: Economic and environmental generation and voyage scheduling of all electric ships in IEEE Transactions on Power Systems
    case 4 
        OPTIONS.G(1,1:3) = [29.16 2185.65 300];
        OPTIONS.G(2,1:3) = [18.83 2330.10 291];
        OPTIONS.G(3,1:3) = [9.37 622.50 210.00];
        OPTIONS.G(4,1:3) = [9.89 618.75 204.00];
        OPTIONS.G(5,1:3) = [26.86 886.50 135.0];
        OPTIONS.G(6,1:3) = [13.67 532.80 111.0];
end


% GHG function parameters
switch parameter_test
    case 1 % 
        OPTIONS.M(1,1:3) = [10 19 170];
        OPTIONS.M(2,1:3) = [28.5 28 190];
        OPTIONS.M(3,1:3) = [10 19 170];
        OPTIONS.M_sl = 30;
        
    case 2 % similar minimum SFC (160), large SFC deviation (250 and 150) between lower bound Pmin and efficient point Peff 
        OPTIONS.M(1,1:3) = [12 60 350];
        OPTIONS.M(2,1:3) = [30 12 170];
        OPTIONS.M(3,1:3) = [12 60 350];
        OPTIONS.M_sl = 30;
        
    % Kanellos: OPMS in IEEE Transactions on Sustainable Energy, Electric power systems research, and Inventions 
    case 3 
        OPTIONS.M(1,1:3) = 3.2*[5.4 61.5 390];
        OPTIONS.M(2,1:3) = 2.5*[13.1 12 430];
        OPTIONS.M(3,1:3) = 3.2*[5.4 63 400];
        OPTIONS.M(4,1:3) = 2.5*[13.5 10 450];
        OPTIONS.M_sl = 30;
        
    % Ce Shang: Economic and environmental generation and voyage scheduling of all electric ships in IEEE Transactions on Power Systems
    case 4 
        OPTIONS.M(1,1:3) = [8383.50 -2070.36 386.10];
        OPTIONS.M(2,1:3) = [956.25 -613.12 362.43];
        OPTIONS.M(3,1:3) = [8383.50 -2070.36 386.10];
        OPTIONS.M(4,1:3) = [956.25 613.12 362.43];
        OPTIONS.M(5,1:3) = [850.50 442.39 124.74];
        OPTIONS.M(6,1:3) = [241.90 235.94 124.74];
        OPTIONS.M_sl = 80;
end

% OPTIONS.C_ss = [90; 45];
OPTIONS.C_ss = 0 * ones(OPTIONS.N_g, 1);
OPTIONS.R_G = OPTIONS.Pg_Max*0.75;
% OPTIONS.R_G = [2 2];
OPTIONS.error = 1e-3;
OPTIONS.Penalty = max([2*OPTIONS.G(1,1)*OPTIONS.Pg_Max(1) + OPTIONS.G(1,2), 2*OPTIONS.G(2,1)*OPTIONS.Pg_Max(2) + OPTIONS.G(2,2)]);
A_max = max(OPTIONS.G(1:2,1));
B_max = max(OPTIONS.G(1:2,2));
Css_max = max(OPTIONS.G(1:2,3));
P_G_max = max(OPTIONS.Pg_Max);
% P_E_max = max(OPTIONS.Pb_Max);

MIN_PRECISE = 0.01;

OPTIONS.Penalty_L = max(2*A_max*P_G_max + B_max, 2*OPTIONS.Xi_E*OPTIONS.E(1,1)*OPTIONS.Pb_Max);
OPTIONS.Penalty_D = 3*(2.2e-3)^(1/3)*OPTIONS.Penalty_L/(OPTIONS.Ppr_Max)^(1/3-1);
OPTIONS.Penalty_xi = max(2*A_max*P_G_max + B_max + (Css_max*(OPTIONS.N_t-1))/MIN_PRECISE, 2*OPTIONS.Xi_E*OPTIONS.E(1,1)*OPTIONS.Pb_Max + OPTIONS.E(1,2)/MIN_PRECISE);

% OPTIONS.varphi = 2/3;
end

%% Optimal algorithm: solve the subproblem which is used to calculate the optimal power of generators and ESMs
function [objval_upperbound, sub_optval, sub_P, sub_dual, reduced_distance] = optimization_subproblem( operation_mode, delta_g, redundant_sw )
global OPTIONS infeasible_flag
% cvx_begin
% cvx_solver SeDuMi
cvx_begin quiet
    variable Ppr(1, OPTIONS.N_t) nonnegative
    variable Pb(OPTIONS.N_e, OPTIONS.N_t)
    variable E(OPTIONS.N_e, OPTIONS.N_t) nonnegative
    variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    variable redundant_sw_s(2, OPTIONS.N_t)
    variable sub_delta_g(OPTIONS.N_g, OPTIONS.N_t)
    variable load_shedding(2, OPTIONS.N_t) nonnegative % it denotes the load shedding amount
    variable reduced_distance nonnegative
    variable Pc(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    variable Pd(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    dual variable dual_delta_g1
    dual variable dual_delta_g2
    dual variable dual_Sp
    dual variable dual_Ss

    minimize( sum( OPTIONS.G(1:OPTIONS.N_g, 1).' * power(Pg(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2), 2) ...
            + sum( OPTIONS.G(1:OPTIONS.N_g, 2).' * power(Pg(1:OPTIONS.N_g, 1:OPTIONS.N_t), 1), 2) ...
            + OPTIONS.Xi_E * sum(OPTIONS.E(1:OPTIONS.N_e, 1).'* power(Pb(1:OPTIONS.N_e, 1:OPTIONS.N_t),2) ,2) ...
            + OPTIONS.Xi_E * sum(OPTIONS.E(1:OPTIONS.N_e, 2).'* ones(OPTIONS.N_e, OPTIONS.N_t) ,2) ...
            + sum(OPTIONS.Penalty_L * sum(load_shedding(1:2, 1:OPTIONS.N_t),2) ,1) ...
            + 1.02*OPTIONS.Penalty_D * reduced_distance ...
            + sum(100 * OPTIONS.Penalty_D * Pc(1, 1:OPTIONS.N_t).' ,1) ...
            + sum(100 * OPTIONS.Penalty_D * Pc(2, 1:OPTIONS.N_t).' ,1) ...
            + sum(100 * OPTIONS.Penalty_D * Pd(1, 1:OPTIONS.N_t).' ,1) ...
            + sum(100 * OPTIONS.Penalty_D * Pd(2, 1:OPTIONS.N_t).' ,1) )

    subject to
        % the range constraints of all the variables
        Pg(1, 1:OPTIONS.N_t) <= sub_delta_g(1, 1:OPTIONS.N_t) * OPTIONS.Pg_Max(1) + Pc(1, 1:OPTIONS.N_t)
        Pg(1, 1:OPTIONS.N_t) >= sub_delta_g(1, 1:OPTIONS.N_t) * OPTIONS.Pg_Min(1) - Pd(1, 1:OPTIONS.N_t)
        Pg(2, 1:OPTIONS.N_t) <= sub_delta_g(2, 1:OPTIONS.N_t) * OPTIONS.Pg_Max(2) + Pc(2, 1:OPTIONS.N_t)
        Pg(2, 1:OPTIONS.N_t) >= sub_delta_g(2, 1:OPTIONS.N_t) * OPTIONS.Pg_Min(2) - Pd(2, 1:OPTIONS.N_t)

        dual_delta_g1 : sub_delta_g(1, 1:OPTIONS.N_t) == delta_g(1, 1:OPTIONS.N_t)
        dual_delta_g2 : sub_delta_g(2, 1:OPTIONS.N_t) == delta_g(2, 1:OPTIONS.N_t)

        % ramping rate power of generators
        Pg(1, 1) <= OPTIONS.R_G(1) + Pc(1, 1)
        Pg(2, 1) <= OPTIONS.R_G(2) + Pc(2, 1)
        for index_t = 2:OPTIONS.N_t
            Pg(1, index_t) - Pg(1, index_t-1) <= OPTIONS.R_G(1)
            Pg(1, index_t) - Pg(1, index_t-1) >= -OPTIONS.R_G(1)
            Pg(2, index_t) - Pg(2, index_t-1) <= OPTIONS.R_G(2)
            Pg(2, index_t) - Pg(2, index_t-1) >= -OPTIONS.R_G(2)
        end
        Pg(1, OPTIONS.N_t) <= OPTIONS.R_G(1) + Pc(1, OPTIONS.N_t)
        Pg(2, OPTIONS.N_t) <= OPTIONS.R_G(2) + Pc(2, OPTIONS.N_t)

        % propulsion power limitation
        Ppr(1, 1:OPTIONS.N_t) <= OPTIONS.Ppr_Max * ones(1, OPTIONS.N_t)

        % ESM limitation
        for index_e = 1:OPTIONS.N_e
            Pb(index_e, 1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1, OPTIONS.N_t)
            Pb(index_e, 1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1, OPTIONS.N_t)
        end
%         Pb(2, 1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1, OPTIONS.N_t)
%         Pb(2, 1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1, OPTIONS.N_t)

        % charging and discharging
        for index_e = 1:OPTIONS.N_e
            E(index_e, 1:OPTIONS.N_t) <= OPTIONS.E_Max(1) * ones(1, OPTIONS.N_t)
            E(index_e, 1:OPTIONS.N_t) >= OPTIONS.E_Min(1) * ones(1, OPTIONS.N_t)
        end
%         E(2, 1:OPTIONS.N_t) <= OPTIONS.E_Max(2) * ones(1, OPTIONS.N_t)
%         E(2, 1:OPTIONS.N_t) >= OPTIONS.E_Min(2) * ones(1, OPTIONS.N_t)

        % ESM output power and the capacity constraints
        % When Pb is positive, ESM is in discharge status,
        % Thus, Pg + Pg is used in power balance comstraint. 
        for index_e = 1:OPTIONS.N_e
            E(index_e, 1) == OPTIONS.E_Max(index_e) - Pb(index_e,1)
%             E(2, 1) == OPTIONS.E_Max(2) - Pb(2,1)
            for t_index = 2:OPTIONS.N_t
                E(index_e, t_index) == E(index_e, t_index-1) - Pb(index_e, t_index)
%                 E(2, t_index) == E(2, t_index-1) - Pb(2, t_index)
            end
        end

        % voyage planning and ESS usage
        switch mod(operation_mode, 4)
            case 0  % only generator scheduling
                Ppr(1, 1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
                Pb(1:OPTIONS.N_e, 1:OPTIONS.N_t) == 0;
            case 1  % generator scheduling & ESMC
                Ppr(1, 1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
            case 2  % generator scheduling & ESMC & PPA
                sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance; 
            case 3  % generator scheduling & ESMC & PPA, load shedding is modelled in the code of next part
                sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance; 
        end
        
        % system power balance & load shedding range in three modes
        if operation_mode <= 3
            % normal operation without load shedding
            for t_index = 1:OPTIONS.N_t
                 OPTIONS.P_L_TIME_off(1, t_index) - sum(load_shedding(1:2, t_index)) + Ppr(1, t_index) == sum( Pg(1:OPTIONS.N_g, t_index) ) + sum(Pb(1:OPTIONS.N_e, t_index) )
            end
            % load shedding range
            if operation_mode == 3
                load_shedding(1, 1:OPTIONS.N_t) + load_shedding(2, 1:OPTIONS.N_t) <= OPTIONS.P_no
            else
                load_shedding(1, 1:OPTIONS.N_t) == zeros(1, OPTIONS.N_t)
                load_shedding(2, 1:OPTIONS.N_t) == zeros(1, OPTIONS.N_t)
            end
        elseif operation_mode <= 7
            % semi-island mode operation with load shedding and
            % reconfiguration (only in mode 7)
            for t_index = 1:OPTIONS.N_t
                redundant_sw_s(1:2, t_index).' * OPTIONS.Coupled_load(:, t_index) + OPTIONS.island1_load(1, t_index) ...
                    - load_shedding(1, t_index) + Ppr(1, t_index) == (Pg(1:OPTIONS.N_g-1, t_index)) + sum(Pb(1:OPTIONS.N_e-1, t_index)) 
                (ones(OPTIONS.N_g,1) - redundant_sw_s(1:OPTIONS.N_g, t_index)).' * OPTIONS.Coupled_load(:, t_index) ...
                    + OPTIONS.island2_load(1, t_index) - load_shedding(2, t_index) == (Pg(OPTIONS.N_g, t_index)) + Pb(2, t_index)
            end
            
            if operation_mode == 7
            % load shedding amount is limited by the demand of non-vital loads
                load_shedding(1, 1:OPTIONS.N_t) <= OPTIONS.island1_non_load(1, 1:OPTIONS.N_t)
                load_shedding(2, 1:OPTIONS.N_t) <= OPTIONS.island2_non_load(1, 1:OPTIONS.N_t)
            else
                load_shedding(1, 1:OPTIONS.N_t) == zeros(1, OPTIONS.N_t)
                load_shedding(2, 1:OPTIONS.N_t) == zeros(1, OPTIONS.N_t)
            end
        elseif operation_mode <= 11 
            % island fault mode operation with load shedding only in mode
            % 11
            for t_index = 1:OPTIONS.N_t
                OPTIONS.island1_load(1, t_index) - load_shedding(1, t_index) + Ppr(1, t_index) == (Pg(1:OPTIONS.N_g-1, t_index)) + sum(Pb(1:OPTIONS.N_e-1, t_index) )
                OPTIONS.island2_load(1, t_index) - load_shedding(2, t_index) == (Pg(OPTIONS.N_g, t_index)) + Pb(OPTIONS.N_e, t_index)
            end
            
            if operation_mode ==11
                % load shedding amount is limited by the demand of non-vital loads
                load_shedding(1, 1:OPTIONS.N_t) <= OPTIONS.island1_non_load(1, 1:OPTIONS.N_t)
                load_shedding(2, 1:OPTIONS.N_t) <= OPTIONS.island2_non_load(1, 1:OPTIONS.N_t)
            else
                load_shedding(1, 1:OPTIONS.N_t) == zeros(1, OPTIONS.N_t)
                load_shedding(2, 1:OPTIONS.N_t) == zeros(1, OPTIONS.N_t)
            end
        end
        
        % the difference of redundant switch is shown in master problem,
        % here only obtain the dual value
        dual_Sp : redundant_sw_s(1:2, :) == redundant_sw(1:2, :)
%         dual_Ss : redundant_sw_s(3:4, :) == redundant_sw(3:4, :)

cvx_end
sub_optval = cvx_optval;

% sub_P.Pg = round(Pg, 5);
% sub_P.Pb = round(Pb, 5);
% sub_P.Pc = round(Pc, 5);
% sub_P.Pd = round(Pd, 5);
% sub_P.Ppr = round(Ppr, 5);
% sub_P.load_shedding = round(load_shedding, 5);

sub_P.Pg = Pg;
sub_P.Pb = Pb;
sub_P.Pc = Pc;
sub_P.Pd = Pd;
sub_P.Ppr = Ppr;
sub_P.load_shedding = load_shedding;

if strcmp(cvx_status, 'Infeasible')
    disp('stop');
    infeasible_flag = 1;
end

%% startup cost: first startup cost and other startup cost
startup_g = ( (sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) - [zeros(OPTIONS.N_g, 1) sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t-1)] ) >=1);
startup_cost = sum(OPTIONS.C_ss(1:OPTIONS.N_g, 1).'* startup_g ,2);
% shutdown = ( (sub_delta_g(1:OPTIONS.N_g,1:OPTIONS.N_t) - [ sub_delta_g(1:OPTIONS.N_g, 2:OPTIONS.N_t) zeros(OPTIONS.N_g,1) ] ) >=1);
% shutdown_cost = OPTIONS.weight(1)*sum(OPTIONS.C_ss(1:OPTIONS.N_g,1).'* shutdown ,2);
one_item_cost = sum(OPTIONS.G(1:OPTIONS.N_g, 3).' * sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) ,2) ;
objval_upperbound = sub_optval + sum(startup_cost) + one_item_cost;

sub_dual.switch(1:2, :) = dual_Sp;
% sub_dual.switch(3:4, :) = dual_Ss;
sub_dual.delta_g(1, :) = dual_delta_g1;
sub_dual.delta_g(2, :) = dual_delta_g2;

if abs(objval_upperbound - 14899.4643324579) <= 2
    disp(objval_upperbound);
end

% disp('upperbound');
end

%% Optimal algorithm: solve the master problem which is used to determine the redundant switches and generator state
function [cvx_optval, master_delta_g, master_redundant_switch, benders_cut ] = optimization_masterproblem( operation_mode, benders_cut_lowerbound, Max_benders_iteration, Ppr_best_solution )
global OPTIONS accelerate_flag  total_P upper_of_lowerbound total_sub total_dual infeasible_flag
% global total_sub total_dual

if Max_benders_iteration == 0
    Max_benders_iteration = 1;
end

if ~exist('benders_cut', 'var') 
    benders_cut = 0;
end

if ~exist('Pd', 'var')
    Pd = 0;
end

% cvx_begin
cvx_solver Mosek
cvx_begin quiet
    variable master_delta_g(OPTIONS.N_g, OPTIONS.N_t) binary
    variable startup_g(OPTIONS.N_g, OPTIONS.N_t) binary
%     variable shutdown_g(OPTIONS.N_g, OPTIONS.N_t) binary
    variable master_redundant_switch(2, OPTIONS.N_t) binary
    variable switch_change(2, OPTIONS.N_t) binary
    variable benders_cut
    minimize( sum(OPTIONS.G(1:OPTIONS.N_g, 3).' * master_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
            + sum(OPTIONS.C_ss(1:OPTIONS.N_g, 1).' * startup_g(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
            + benders_cut )

    subject to
        % startup detect
        startup_g(1:OPTIONS.N_g, 1) >= master_delta_g(1:OPTIONS.N_g, 1) - zeros(OPTIONS.N_g, 1)
        for t_index = 2 : OPTIONS.N_t
            startup_g(1:OPTIONS.N_g, t_index) >= master_delta_g(1:OPTIONS.N_g, t_index) - master_delta_g(1:OPTIONS.N_g, t_index-1)
        end
        
        %  generator operation time (rising edge)
        for t_index = 1 : OPTIONS.N_t-OPTIONS.Tmin_g_st+1
            sum(master_delta_g(1:OPTIONS.N_g, t_index:t_index + OPTIONS.Tmin_g_st-1), 2) >= OPTIONS.Tmin_g_st * startup_g(1:OPTIONS.N_g, t_index)
        end
        % edge condition in OPTIONS.N_t - OPTIONS.Tmin_g_st+2 to OPTIONS.N_t
        for t_index = OPTIONS.N_t-OPTIONS.Tmin_g_st+2 : OPTIONS.N_t
            sum(master_delta_g(1:OPTIONS.N_g, t_index:OPTIONS.N_t), 2) >= OPTIONS.Tmin_g_st * startup_g(1:OPTIONS.N_g, t_index)
        end
        
%         % shutdown detect
%         for t_index = 1 : OPTIONS.N_t-1
%             shutdown_g(1:OPTIONS.N_g, t_index) >= master_delta_g(1:OPTIONS.N_g, t_index) - master_delta_g(1:OPTIONS.N_g, t_index+1)
%         end
%         shutdown_g(1:OPTIONS.N_g, OPTIONS.N_t) >= master_delta_g(1:OPTIONS.N_g, OPTIONS.N_t) - zeros(OPTIONS.N_g, 1)
        
%         %  generator shutdown time (falling edge)
%         for t_index = 1 : OPTIONS.N_t-OPTIONS.Tmin_g_sd
%             sum(ones(OPTIONS.N_g, OPTIONS.Tmin_g_sd) - master_delta_g(1:OPTIONS.N_g, t_index:t_index + OPTIONS.Tmin_g_sd-1), 2) ...
%                 >= OPTIONS.Tmin_g_sd * shutdown_g(1:OPTIONS.N_g, t_index)
%         end
%         % edge condition in OPTIONS.N_t - OPTIONS.Tmin_g_sd+1 to OPTIONS.N_t
%         for t_index = OPTIONS.N_t-OPTIONS.Tmin_g_sd+1 : OPTIONS.N_t
%             sum(ones(OPTIONS.N_g, OPTIONS.Tmin_g_sd) - master_delta_g(1:OPTIONS.N_g, t_index:OPTIONS.N_t), 2) ...
%                 >= OPTIONS.Tmin_g_sd * shutdown_g(1:OPTIONS.N_g, t_index)
%         end

        %  Redundent_switch
%         master_redundant_switch(1, :) + master_redundant_switch(2, :) == ones(1, OPTIONS.N_t)
%         master_redundant_switch(3, :) + master_redundant_switch(4, :) == ones(1, OPTIONS.N_t)
        
        % in the normal mode, redundant switches are determined.
        if operation_mode <= 3 
            if operation_mode == 3
                
                % switch change detect
                switch_change(1, 1) >= (master_redundant_switch(1, 1) - 1)
                switch_change(2, 1) >= (master_redundant_switch(2, 1) - 0)
                switch_change(1, 2:OPTIONS.N_t) >= (master_redundant_switch(1, 2:OPTIONS.N_t) - master_redundant_switch(1, 1:OPTIONS.N_t-1))
                switch_change(2, 2:OPTIONS.N_t) >= (master_redundant_switch(2, 2:OPTIONS.N_t) - master_redundant_switch(2, 1:OPTIONS.N_t-1))

                %  Redundent switch time
                for t_index = 1 : OPTIONS.N_t-OPTIONS.Tmin_sw+1
                    sum(master_redundant_switch(1, t_index:t_index + OPTIONS.Tmin_sw-1), 2) >= OPTIONS.Tmin_sw * switch_change(1, t_index)
                    sum(master_redundant_switch(2, t_index:t_index + OPTIONS.Tmin_sw-1), 2) >= OPTIONS.Tmin_sw * switch_change(2, t_index)
                end
                for t_index = OPTIONS.N_t-OPTIONS.Tmin_sw+2 : OPTIONS.N_t
                    sum(master_redundant_switch(1, t_index:OPTIONS.N_t), 2) >= OPTIONS.Tmin_sw * switch_change(1, t_index)
                    sum(master_redundant_switch(2, t_index:OPTIONS.N_t), 2) >= OPTIONS.Tmin_sw * switch_change(2, t_index)
                end
            else
                master_redundant_switch(1, 1:OPTIONS.N_t) == total_sub(1).redundant_sw(1, 1:OPTIONS.N_t)
                master_redundant_switch(2, 1:OPTIONS.N_t) == total_sub(1).redundant_sw(2, 1:OPTIONS.N_t)
            end
%             master_redundant_switch(3, 1:OPTIONS.N_t) == total_sub(1).redundant_sw(3, 1:OPTIONS.N_t)
%             master_redundant_switch(4, 1:OPTIONS.N_t) == total_sub(1).redundant_sw(4, 1:OPTIONS.N_t)
            
            % benders cuts
            for index_benders = 1:Max_benders_iteration
                benders_cut >=  total_sub(index_benders).sub_optval ...
                                + total_dual(index_benders).delta_g(1, 1:OPTIONS.N_t) ...
                                  * (master_delta_g(1, :).' - total_sub(index_benders).delta_g(1, 1:OPTIONS.N_t).') ...
                                + total_dual(index_benders).delta_g(2, 1:OPTIONS.N_t) ...
                                  * (master_delta_g(2, :).' - total_sub(index_benders).delta_g(2, 1:OPTIONS.N_t).');
            end
        % in the fault mode, redundant switches need to be reconfigurated.
        elseif operation_mode <= 7
            if operation_mode == 7
                % switch change detect
                switch_change(1, 1) >= (master_redundant_switch(1, 1) - 1)
                switch_change(2, 1) >= (master_redundant_switch(2, 1) - 0)
                switch_change(1, 2:OPTIONS.N_t) >= (master_redundant_switch(1, 2:OPTIONS.N_t) - master_redundant_switch(1, 1:OPTIONS.N_t-1))
                switch_change(2, 2:OPTIONS.N_t) >= (master_redundant_switch(2, 2:OPTIONS.N_t) - master_redundant_switch(2, 1:OPTIONS.N_t-1))

                %  Redundent switch time
                for t_index = 1 : OPTIONS.N_t-OPTIONS.Tmin_sw+1
                    sum(master_redundant_switch(1, t_index:t_index + OPTIONS.Tmin_sw-1), 2) >= OPTIONS.Tmin_sw * switch_change(1, t_index)
                    sum(master_redundant_switch(2, t_index:t_index + OPTIONS.Tmin_sw-1), 2) >= OPTIONS.Tmin_sw * switch_change(2, t_index)
                end
                for t_index = OPTIONS.N_t-OPTIONS.Tmin_sw+2 : OPTIONS.N_t
                    sum(master_redundant_switch(1, t_index:OPTIONS.N_t), 2) >= OPTIONS.Tmin_sw * switch_change(1, t_index)
                    sum(master_redundant_switch(2, t_index:OPTIONS.N_t), 2) >= OPTIONS.Tmin_sw * switch_change(2, t_index)
                end
            else
                master_redundant_switch(1, 1:OPTIONS.N_t) == total_sub(1).redundant_sw(1, 1:OPTIONS.N_t)
                master_redundant_switch(2, 1:OPTIONS.N_t) == total_sub(1).redundant_sw(2, 1:OPTIONS.N_t)
            end
            
            % benders cuts
            for index_benders = 1:Max_benders_iteration
                benders_cut >=  total_sub(index_benders).sub_optval ...
                                + total_dual(index_benders).delta_g(1, 1:OPTIONS.N_t) ...
                                  * (master_delta_g(1, :).' - total_sub(index_benders).delta_g(1, 1:OPTIONS.N_t).') ...
                                + total_dual(index_benders).delta_g(2, 1:OPTIONS.N_t) ...
                                  * (master_delta_g(2, :).' - total_sub(index_benders).delta_g(2, 1:OPTIONS.N_t).') ... 
                                + total_dual(index_benders).switch(1, 1:OPTIONS.N_t) ...
                                  * (master_redundant_switch(1, :).' - total_sub(index_benders).redundant_sw(1, 1:OPTIONS.N_t).')...
                                + total_dual(index_benders).switch(2, 1:OPTIONS.N_t) ...
                                  * ((master_redundant_switch(2, :)).' - total_sub(index_benders).redundant_sw(2, 1:OPTIONS.N_t).');
%                                 + total_dual(index_benders).switch(3, 1:OPTIONS.N_t) ...
%                                   *(master_redundant_switch(3, :).' - total_sub(index_benders).redundant_sw(3, 1:OPTIONS.N_t).')...
%                                 + total_dual(index_benders).switch(4, 1:OPTIONS.N_t) ...
%                                   *((master_redundant_switch(4, :)).' - total_sub(index_benders).redundant_sw(4, 1:OPTIONS.N_t).')
            end        % in the fault mode, redundant switches need to be reconfigurated.
        elseif operation_mode <= 11
%             switch_change(1, 1:OPTIONS.N_t-1) >= (master_redundant_switch(1, 2:OPTIONS.N_t) - master_redundant_switch(1, 1:OPTIONS.N_t-1))
%             switch_change(2, 1:OPTIONS.N_t-1) >= (master_redundant_switch(2, 2:OPTIONS.N_t) - master_redundant_switch(2, 1:OPTIONS.N_t-1))
            
            % benders cuts. In this case, there are three generators.
            for index_benders = 1:Max_benders_iteration
                benders_cut >=  total_sub(index_benders).sub_optval ...
                                + total_dual(index_benders).delta_g(1, 1:OPTIONS.N_t) ...
                                  * (master_delta_g(1, :).' - total_sub(index_benders).delta_g(1, 1:OPTIONS.N_t).') ...
                                + total_dual(index_benders).delta_g(2, 1:OPTIONS.N_t) ...
                                  * (master_delta_g(2, :).' - total_sub(index_benders).delta_g(2, 1:OPTIONS.N_t).') ...
                                + total_dual(index_benders).delta_g(3, 1:OPTIONS.N_t) ...
                                  * (master_delta_g(3, :).' - total_sub(index_benders).delta_g(3, 1:OPTIONS.N_t).');
            end
        end

        switch accelerate_flag
            case 1
            case 2
                % low bound of master problem
                sum(OPTIONS.G(1:OPTIONS.N_g, 3).' * master_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
                + sum(OPTIONS.C_ss(1:OPTIONS.N_g, 1).' * startup_g(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
                + benders_cut >= upper_of_lowerbound
                % speedup constraints: power range & lower bound
                if operation_mode <=3 % normal mode
%                     master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) ...
                    OPTIONS.Pg_Max(1, 1:OPTIONS.N_g)*master_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) ...
                        +  OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.P_L_TIME_off + Ppr_best_solution
%                             +  OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.P_L_TIME_off + total_P(Max_benders_iteration).Ppr
                        
%                     master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) ...
                    OPTIONS.Pg_Min(1, 1:OPTIONS.N_g)*master_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) ...
                        +  OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.P_total_vs + 0 % Ppr
                elseif operation_mode <= 7 % semi-island fault mode
                    if Max_benders_iteration >= 1
                        master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) ...
                            >= OPTIONS.island1_max + Ppr_best_solution
%                             >= OPTIONS.island1_max + total_P(Max_benders_iteration).Ppr
                        master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island1_min + 0 % Ppr
                        master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.island2_max
                        master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island2_min + 0 % Ppr
                    end
                    % ppr in island 1, there must have a worked generator.
%                     master_delta_g(1, 1:OPTIONS.N_t) >= master_delta_g(2, 1:OPTIONS.N_t)
                elseif operation_mode <= 11 % island fault mode
                    if Max_benders_iteration >= 1
%                         master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) ...
                        OPTIONS.Pg_Max(1, 1:2)*master_delta_g(1:2, 1:OPTIONS.N_t) >= OPTIONS.island1_max + Ppr_best_solution
%                             >= OPTIONS.island1_max + total_P(Max_benders_iteration).Ppr
                        OPTIONS.Pg_Min(1, 1:2)*master_delta_g(1:2, 1:OPTIONS.N_t) <= OPTIONS.island1_min
%                         master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island1_min + 0 % Ppr
                        master_delta_g(3, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(3) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.island2_max
                        master_delta_g(3, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(3) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island2_min + 0 % Ppr
                    end
                    % ppr in island 1, there must have a worked generator.
%                     sum(master_delta_g(1:2, 1:OPTIONS.N_t), 1)/2 >= master_delta_g(3, 1:OPTIONS.N_t)
                end
                
                % speedup constraint: lower bound
%                 benders_cut >= benders_cut_lowerbound
            case 3
                % speedup constraints: power range
                if operation_mode <=3 % normal mode
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) ...
                        +  OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.P_L_TIME_off + Ppr_best_solution
%                         +  OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.P_L_TIME_off + total_P(Max_benders_iteration).Ppr

                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) ...
                        +  OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.P_total_vs + 0 % Ppr
                elseif operation_mode <= 7 % fault mode
                    if Max_benders_iteration >= 1
                        master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) ...
                            >= OPTIONS.island1_max + Ppr_best_solution
%                             >= OPTIONS.island1_max + total_P(Max_benders_iteration).Ppr
                        
                        master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island1_min + 0 % Ppr
                        master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.island2_max
                        master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island2_min + 0 % Ppr
                    end
                end
            case 4
                % speedup constraint: lower bound
                benders_cut >= benders_cut_lowerbound
%             otherwise % non-constriant
        end
        
cvx_end
        
if strcmp(cvx_status, 'Infeasible')
    disp('Master Infeasible');
    infeasible_flag = 1;
end
if strcmp(cvx_status, 'Failed')
    disp('Master Failed');
    infeasible_flag = 1;
end

end

%% LNBD algorithm: sovlve the subproblem t which is used to calculate the optimal power of generators and ESMs at time t  
function [objval_upperbound, sub_optval, sub_P, sub_dual, reduced_distance] = optimization_subproblem_t( operation_mode, delta_g, redundant_sw )
% in this function, the return variable: suboptval is different with that
% in the function of optimal algorithm.

global OPTIONS with_iteration_D infeasible_flag
% XXX_on denotes the variable that calculated in this iteration.
% inital the rest parameters
Pg = zeros(OPTIONS.N_g, OPTIONS.N_t);
Pb = zeros(OPTIONS.N_e, OPTIONS.N_t);
E = zeros(OPTIONS.N_e, OPTIONS.N_t);
Pc = zeros(OPTIONS.N_g, OPTIONS.N_t);
Pd = zeros(OPTIONS.N_g, OPTIONS.N_t);
Ppr = zeros(1, OPTIONS.N_t);
load_shedding = zeros(2, OPTIONS.N_t);

dual_Sp = zeros(2, OPTIONS.N_t);
dual_Ss = zeros(2, OPTIONS.N_t);
dual_delta_g1 = zeros(1, OPTIONS.N_t);
dual_delta_g2 = zeros(1, OPTIONS.N_t);

Rest_Distance = OPTIONS.Distance;
Distance_obj_on = zeros(1, OPTIONS.N_t);
Rest_ppr_avg = OPTIONS.P_pr_avg;
D_d = zeros(1, OPTIONS.N_t);
sub_optval_on = 0;
sub_optval_t = zeros(1, OPTIONS.N_t+1);
sub_optval = zeros(2, OPTIONS.N_t+1);

if operation_mode <= 3
    p_ESM_avg = (sum(OPTIONS.E_Max) - sum(OPTIONS.E_Min))/OPTIONS.N_t;
elseif operation_mode <= 11
    p_ESM_avg(1, 1) = (OPTIONS.E_Max(1) - OPTIONS.E_Min(1))/OPTIONS.N_t;
    p_ESM_avg(2, 1) = (OPTIONS.E_Max(2) - OPTIONS.E_Min(2))/OPTIONS.N_t;
end

% update the rest parameters: rest distance, rest average ppr, and the
% power deviation of load at next time slot.
for index_loop = 1:10
    sub_optval_on = 0;
    for index_time = 1:OPTIONS.N_t
    %     tic;
        switch index_time
            % all the load demand must be meet at the last time slot (N_t). 
            % So varphi = 1
            case OPTIONS.N_t
                sub_varphi_Pl = OPTIONS.varphi_Pl;
                sub_varphi_Ppr = 0.0;

            % the load demand at other time slot is meet proportionally.
            otherwise
                sub_varphi_Pl = OPTIONS.varphi_Pl;
                sub_varphi_Ppr = OPTIONS.varphi_Ppr;
        end

        if index_time == 1
            E_last_time = OPTIONS.E_Max.';
        else
            E_last_time = E(1:OPTIONS.N_e, index_time-1);
        end

        if index_loop == 1
            [upper_bound_Pb, Distance_slot_obj] = generate_upperbound(operation_mode, Rest_ppr_avg, sub_varphi_Ppr, ...
                             sub_varphi_Pl, E_last_time, index_time, p_ESM_avg, delta_g);
        else
            [upper_bound_Pb] = generate_upperbound(operation_mode, Rest_ppr_avg, sub_varphi_Ppr, ...
                             sub_varphi_Pl, E_last_time, index_time, p_ESM_avg, delta_g);
            Distance_slot_obj = Distance_obj_on(index_time);
        end

    %     tic;
        % cvx_solver SeDuMi
        cvx_begin quiet
            variable Pg_on(OPTIONS.N_g) nonnegative
            variable Pb_on(2,1)
            variable E_on(2,1) nonnegative
            variable Ppr_on(1) nonnegative
            variable sub_delta_g(OPTIONS.N_g)
            variable redundant_sw_s(4,1)
            variable load_shedding_on(2, 1) nonnegative
            variable reduced_distance_on nonnegative
            variable Pc_on(OPTIONS.N_g) nonnegative
            variable Pd_on(OPTIONS.N_g) nonnegative
            variable PcR(OPTIONS.N_g) nonnegative
            variable PdR(OPTIONS.N_g) nonnegative
            dual variable temp_dual_delta_g1
            dual variable temp_dual_delta_g2
            dual variable temp_dual_Sp
            dual variable temp_dual_Ss
            minimize( OPTIONS.G(1:OPTIONS.N_g, 1).' * power(Pg_on(1:OPTIONS.N_g, 1), 2) ...
                    + OPTIONS.G(1:OPTIONS.N_g, 2).' * power(Pg_on(1:OPTIONS.N_g, 1), 1) ...
                    + OPTIONS.Xi_E * OPTIONS.E(1:OPTIONS.N_e, 1).' * power(Pb_on(1:OPTIONS.N_e, 1), 2) ...
                    + OPTIONS.Xi_E * sum(OPTIONS.E(1:OPTIONS.N_e, 2), 1)...
                    + sum(OPTIONS.Penalty_L * sum(load_shedding_on(1:OPTIONS.N_g, 1), 2) ,1) ...
                    + 1.02 * OPTIONS.Penalty_D * reduced_distance_on ...
                    + 10 * OPTIONS.Penalty_D * Pc_on(1) ...
                    + 10 * OPTIONS.Penalty_D * Pc_on(2) ...
                    + 10 * OPTIONS.Penalty_D * Pd_on(1) ...
                    + 10 * OPTIONS.Penalty_D * Pd_on(2) )

            subject to
                % the range constraints of all the variables
                Pg_on(1, 1) <= sub_delta_g(1, 1) * OPTIONS.Pg_Max(1) + Pc_on(1)
                Pg_on(1, 1) >= sub_delta_g(1, 1) * OPTIONS.Pg_Min(1) - Pd_on(1)
                Pg_on(2, 1) <= sub_delta_g(2, 1) * OPTIONS.Pg_Max(2) + Pc_on(2)
                Pg_on(2, 1) >= sub_delta_g(2, 1) * OPTIONS.Pg_Min(2) - Pd_on(2)

                temp_dual_delta_g1 : sub_delta_g(1,1) == delta_g(1,index_time)
                temp_dual_delta_g2 : sub_delta_g(2,1) == delta_g(2,index_time)

                % ramping rate power of generators
                if index_time == 1
                    Pg_on(1, 1) <= OPTIONS.R_G(1) %+ PcR(1)
                    Pg_on(2, 1) <= OPTIONS.R_G(2) %+ PcR(2)
                elseif index_time == OPTIONS.N_t
                    Pg_on(1, 1) - Pg(1, index_time-1) >= -OPTIONS.R_G(1) %- PdR(1)
                    Pg_on(2, 1) - Pg(2, index_time-1) >= -OPTIONS.R_G(2) %- PdR(2)
                    0 - Pg_on(1, 1) >= -OPTIONS.R_G(1) %- PcR(1)
                    0 - Pg_on(2, 1) >= -OPTIONS.R_G(2) %- PcR(2)
                else
                    Pg_on(1, 1) - Pg(1, index_time-1) <= OPTIONS.R_G(1) %+ PcR(1)
                    Pg_on(1, 1) - Pg(1, index_time-1) >= -OPTIONS.R_G(1) %- PdR(1) 
                    Pg_on(2, 1) - Pg(2, index_time-1) <= OPTIONS.R_G(2) %+ PcR(2)
                    Pg_on(2, 1) - Pg(2, index_time-1) >= -OPTIONS.R_G(2) %- PdR(2)
                end

                % propulsion power limitation
                Ppr_on(1,1) <= OPTIONS.Ppr_Max

                % ESM limitation constraints
                % the basic constraints of ESMs
                Pb_on(1,1) <= OPTIONS.Pb_Max
                Pb_on(1,1) >= OPTIONS.Pb_Min 
                Pb_on(2,1) <= OPTIONS.Pb_Max
                Pb_on(2,1) >= OPTIONS.Pb_Min
                % the upper bound in suboptimal algorithm
                if operation_mode <= 3
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
    %                 Pb_on(1,1) + Pb_on(2,1) <= upper_bound_ESM_P
                    Pb_on(1,1) + Pb_on(2,1) <= upper_bound_Pb
                elseif operation_mode <= 7
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
    %                 Pb_on(1,1) <= upper_bound_ESM1P
    %                 Pb_on(2,1) <= upper_bound_ESM2P
                    Pb_on(1,1) <= upper_bound_Pb(1)
                    Pb_on(2,1) <= upper_bound_Pb(2)
                end

                % charging and discharging
                E_on(1, 1) <= OPTIONS.E_Max(1)
                E_on(2, 1) <= OPTIONS.E_Max(2)
                E_on(1, 1) >= OPTIONS.E_Min(1)
                E_on(2, 1) >= OPTIONS.E_Min(2)

                % ESM output power and the capacity constraints
                if index_time == 1
                    E_on(1,1) == OPTIONS.E_Max(1) - Pb_on(1,1)
                    E_on(2,1) == OPTIONS.E_Max(2) - Pb_on(2,1)
                else
                    E_on(1,1) == E(1, index_time-1) - Pb_on(1,1)
                    E_on(2,1) == E(2, index_time-1) - Pb_on(2,1)
                end

                % system power balance & load shedding
                if operation_mode <= 3
                    OPTIONS.P_L_TIME_off(1, index_time) - sum(load_shedding_on(1:2, 1)) + Ppr_on(1, 1) ...
                         == sum( Pg_on(1:OPTIONS.N_g, 1) ) + sum(Pb_on(1:OPTIONS.N_e, 1))

                    % load shedding range
                    if operation_mode == 3
                        sum(load_shedding_on(1:2, 1)) <= OPTIONS.P_no(index_time)
                    else
                        load_shedding_on(1) == zeros(1, 1)
                        load_shedding_on(2) == zeros(1, 1)
                    end
                elseif operation_mode <= 7
                    for t_index = index_time
                        redundant_sw_s(1:2, 1).' * OPTIONS.Coupled_load(:, index_time) + OPTIONS.island1_load(1, t_index) ...
                            - load_shedding_on(1, 1) + Ppr_on(1, 1) == (Pg_on(1, 1)) + Pb_on(1, 1) 
                        redundant_sw_s(3:4, 1).'* OPTIONS.Coupled_load(:, index_time) + OPTIONS.island2_load(1, t_index) ...
                            - load_shedding_on(2, 1) == (Pg_on(2, 1)) + Pb_on(2, 1)
                    end

                    % load shedding range
                    if operation_mode == 7 
                        load_shedding_on(1) <= OPTIONS.island1_non_load(1, index_time)
                        load_shedding_on(2) <= OPTIONS.island2_non_load(1, index_time)
                    else
                        load_shedding_on(1) == zeros(1, 1)
                        load_shedding_on(2) == zeros(1, 1)
                    end
                elseif operation_mode <= 11
                    for t_index = index_time
                        OPTIONS.island1_load(1, t_index) - load_shedding_on(1, 1) + Ppr_on(1, 1) == (Pg_on(1, 1)) + Pb_on(1, 1) 
                        OPTIONS.island2_load(1, t_index) - load_shedding_on(2, 1) == (Pg_on(2, 1)) + Pb_on(2, 1)
                    end

                    % load shedding range
                    if operation_mode == 11 
                        load_shedding_on(1) <= OPTIONS.island1_non_load(1, index_time)
                        load_shedding_on(2) <= OPTIONS.island2_non_load(1, index_time)
                    else
                        load_shedding_on(1) == zeros(1, 1)
                        load_shedding_on(2) == zeros(1, 1)
                    end
                end

                temp_dual_Sp : redundant_sw_s(1:2, 1) == redundant_sw(1:2, index_time)
                temp_dual_Ss : redundant_sw_s(3:4, 1) == redundant_sw(3:4, index_time)

                % voyage planning
                switch mod(operation_mode, 4)
                    case 0 % Only GS
                        Ppr_on(1, 1) == OPTIONS.P_pr_avg;
                        Pb_on(1:2, 1) == 0;
                    case 1 % GS & ESSC
                        Ppr_on(1, 1) == OPTIONS.P_pr_avg;
                    case 2 % GS & ESSC & PPA
                        (Ppr_on(1, 1)./2.2e-3).^(1/3) >= Distance_slot_obj - reduced_distance_on; 
                    case 3 % GS & ESSC & PPA (load shedding )
                        (Ppr_on(1, 1)./2.2e-3).^(1/3) >= Distance_slot_obj - reduced_distance_on; 
                end

        cvx_end
    %     test_time = toc;
    %     fprintf('\nSP part 2 at %d is %5.2f s \n\n', index_time, test_time);

        if strcmp(cvx_status, 'Infeasible')
            disp('Sub Infeasible');
            infeasible_flag = 1;
            break;
        end

        Pg(1:OPTIONS.N_g, index_time) = Pg_on(1:OPTIONS.N_g, 1);
        Pb(1:OPTIONS.N_e, index_time) = Pb_on(1:OPTIONS.N_g, 1);
        E(1:OPTIONS.N_e, index_time) = E_on(1:OPTIONS.N_g, 1);
        Pc(1:OPTIONS.N_g, index_time) = Pc_on(1:OPTIONS.N_g, 1);
        Pd(1:OPTIONS.N_g, index_time) = Pd_on(1:OPTIONS.N_g, 1);
        Ppr(1, index_time) = Ppr_on(1, 1);
        load_shedding(1:2, index_time) = load_shedding_on(1:2, 1);

        % with_iteration_D = 0: D(t) is only adjusted by load change and
        % varphi_Pl. sub_optval equals to cvx_optval
        % with_iteration_D = 1: D(t) is iterated at each time slot based on the
        % load change, varphi_Pl and D_d of last time slot. 

        if with_iteration_D == 0
    %         Rest_Distance =  Rest_Distance - Distance_slot_obj;
    %         Rest_velocity_avg = Rest_Distance/(OPTIONS.N_t-index_time);
            Rest_ppr_avg = OPTIONS.P_pr_avg;
            sub_optval_t(index_time) = cvx_optval;
            sub_optval_on = sub_optval_on + cvx_optval;
            D_d(index_time) = reduced_distance_on;
            Distance_obj_on(index_time) = Distance_slot_obj - reduced_distance_on;
            if index_time == OPTIONS.N_t
                Rest_Distance = sum(D_d);

            end
        elseif with_iteration_D == 1
            Rest_Distance =  Rest_Distance - Distance_slot_obj + reduced_distance_on;
            if index_time < OPTIONS.N_t
                Rest_velocity_avg = Rest_Distance/(OPTIONS.N_t-index_time);
            end 
            Rest_ppr_avg = (Rest_velocity_avg).^3*2.2e-3;

            % the penalty cost of reduced distance is calculated at the last
            % iteration, thus it must subtract the penalty cost in the total cost
            %  before the last iteration.
            if index_time < OPTIONS.N_t
                sub_optval_t(index_time) = cvx_optval - 1.02 * OPTIONS.Penalty_D * reduced_distance_on;
                sub_optval_on = sub_optval_on + sub_optval_t(index_time);
            else
                sub_optval_on = sub_optval_on + cvx_optval;
                sub_optval_t(index_time) = cvx_optval - 1.02 * OPTIONS.Penalty_D * reduced_distance_on;
                sub_optval_t(index_time+1) = 1.02 * OPTIONS.Penalty_D * reduced_distance_on;
            end
        elseif with_iteration_D == 2
        end

        dual_Sp(1:2, index_time) = temp_dual_Sp;
        dual_Ss(1:2, index_time) = temp_dual_Ss;
        dual_delta_g1(1, index_time) = temp_dual_delta_g1;
        dual_delta_g2(1, index_time) = temp_dual_delta_g2;
    end
    if length(find(D_d <= 0.01)) >= 1
        Distance_obj_on = Distance_obj_on + Rest_Distance/OPTIONS.N_t;
    else
        break;
    end
    if infeasible_flag == 1
        break;
    end
end

reduced_distance = Rest_Distance;
sub_delta_g(1:2, 1:OPTIONS.N_t) = delta_g(1:2, 1:OPTIONS.N_t);
sub_optval(1, 1) = sub_optval_on;
sub_optval(2, :) = sub_optval_t;

%% startup cost: first startup cost and other startup cost
startup_g = ( (sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) - [zeros(OPTIONS.N_g, 1) sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t-1)] ) >=1);
startup_cost = sum(OPTIONS.C_ss(1:OPTIONS.N_g, 1).' * startup_g ,2);
one_item_cost = sum(OPTIONS.G(1:OPTIONS.N_g, 3).' * sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) ,2) ;
objval_upperbound = sub_optval_on + sum(startup_cost) + one_item_cost;

sub_P.Pg = Pg;
sub_P.Pb = Pb;
sub_P.Pc = Pc;
sub_P.Pd = Pd;
sub_P.Ppr = Ppr;
sub_P.load_shedding = load_shedding;

sub_dual.switch(1:2, 1:OPTIONS.N_t) = dual_Sp;
sub_dual.switch(3:4, 1:OPTIONS.N_t) = dual_Ss;
sub_dual.delta_g(1, 1:OPTIONS.N_t) = dual_delta_g1;
sub_dual.delta_g(2, 1:OPTIONS.N_t) = dual_delta_g2;

if abs(objval_upperbound - 18085) < 10
    disp(objval_upperbound);
end

end

%% Subfunction of LNBD algorithm: generate the upper bound of output power of ESMs based on suboptimal decomposition strategy
function [upper_bound_Pb, Distance_slot_obj] = generate_upperbound(operation_mode, Rest_ppr_avg, sub_varphi_Ppr, sub_varphi_Pl, E_last_time, index_time, p_ESM_avg, delta_g)
global OPTIONS 
    % upper bound in different operation mode 
    if operation_mode <= 3
        Distance_slot_obj = (( Rest_ppr_avg - sub_varphi_Ppr * OPTIONS.Delta_PL(index_time) )/2.2e-3)^(1/3);
        
        % ESM output power and the capacity constraints
        if index_time == 1
            rest_pmax_ESM = sum(OPTIONS.E_Max) - sum(OPTIONS.E_Min);
            rest_pmin_ESM = sum(OPTIONS.E_Max) - sum(OPTIONS.E_Max);
        else
            rest_pmax_ESM = sum(E_last_time, 1) - sum(OPTIONS.E_Min);
            rest_pmin_ESM = sum(E_last_time, 1) - sum(OPTIONS.E_Max);
        end
        if sum(delta_g(1:OPTIONS.N_g, index_time)) > 0
            if OPTIONS.Delta_PL(index_time) >= 0
                upper_bound_ESM_P = roundn(min(p_ESM_avg + sub_varphi_Pl * OPTIONS.Delta_PL(index_time), rest_pmax_ESM), -2);
            else
                upper_bound_ESM_P = roundn(max(p_ESM_avg + sub_varphi_Pl * OPTIONS.Delta_PL(index_time), rest_pmin_ESM), -2);
            end
        else
            upper_bound_ESM_P = sum(OPTIONS.E_Max) - sum(OPTIONS.E_Min);
        end
        upper_bound_Pb = upper_bound_ESM_P;
    elseif operation_mode <= 7
        % propulsion power modules locate in island part 1, 
        % power adjustment only uses OPTIONS.Delta_PL_island1
        Distance_slot_obj = ((Rest_ppr_avg - sub_varphi_Ppr * OPTIONS.Delta_PL_island1(index_time) )/2.2e-3 )^(1/3);
        
        % Maximum charge and discharge power of ESM in different islands
        % The initial capacity of ESM is the maximum value (OPTIONS.E_Max)
        if index_time == 1
            rest_pmax_ESM1 = OPTIONS.E_Max(1) - OPTIONS.E_Min(1);
            rest_pmax_ESM2 = OPTIONS.E_Max(2) - OPTIONS.E_Min(2);
            rest_pmin_ESM1 = OPTIONS.E_Max(1) - OPTIONS.E_Max(1);
            rest_pmin_ESM2 = OPTIONS.E_Max(2) - OPTIONS.E_Max(2);
        else
            rest_pmax_ESM1 = E_last_time(1, 1) - OPTIONS.E_Min(1);
            rest_pmax_ESM2 = E_last_time(2, 1) - OPTIONS.E_Min(2);
            rest_pmin_ESM1 = E_last_time(1, 1) - OPTIONS.E_Max(1);
            rest_pmin_ESM2 = E_last_time(2, 1) - OPTIONS.E_Max(2);
        end
        % charge and discharge bound in island 1
        if sum(delta_g(1,index_time))>=1
            if OPTIONS.Delta_PL_island1(index_time) >= 0
                upper_bound_ESM1P = roundn(min(p_ESM_avg(1, 1) + ...
                    sub_varphi_Pl * OPTIONS.Delta_PL_island1(index_time), rest_pmax_ESM1), -2);
            else
                upper_bound_ESM1P = roundn(max(p_ESM_avg(1, 1) + ...
                    sub_varphi_Pl * OPTIONS.Delta_PL_island1(index_time), rest_pmin_ESM1), -2);
            end
        else
            upper_bound_ESM1P = OPTIONS.E_Max(1) - OPTIONS.E_Min(1);
        end
        % charge and discharge bound in i
        if sum(delta_g(2,index_time))>=1
            if OPTIONS.Delta_PL_island2(index_time) >= 0
                upper_bound_ESM2P = roundn(min(p_ESM_avg(2, 1) + ...
                    sub_varphi_Pl * OPTIONS.Delta_PL_island2(index_time), rest_pmax_ESM2), -2);
            else
                upper_bound_ESM2P = roundn(max(p_ESM_avg(2, 1) + ...
                    sub_varphi_Pl * OPTIONS.Delta_PL_island2(index_time), rest_pmin_ESM2), -2);
            end
        else
            upper_bound_ESM2P = OPTIONS.E_Max(2) - OPTIONS.E_Min(2);
        end
        upper_bound_Pb(1) = upper_bound_ESM1P;
        upper_bound_Pb(2) = upper_bound_ESM2P;
    end
    
end

%% LNBD algorithm: the master problem which is used to determine the redundant switches and generator state
function [master_optval, master_delta_g, master_redundant_switch, benders_cut ] = optimization_masterproblem_t( operation_mode, benders_cut_lowerbound, Max_benders_iteration, Ppr_best_solution)
global OPTIONS accelerate_flag  total_P upper_of_lowerbound total_sub total_dual with_iteration_D infeasible_flag
% global total_sub total_dual

if Max_benders_iteration == 0
    Max_benders_iteration = 1;
end

if ~exist('benders_cut', 'var') 
    benders_cut = 0;
end

if ~exist('Pd', 'var')
    Pd = 0;
end

master_delta_g = zeros(OPTIONS.N_g, OPTIONS.N_t);
master_redundant_switch = zeros(OPTIONS.N_g, OPTIONS.N_t);
sub_optval_on = 0;

for index_time  = 1:OPTIONS.N_t
    % cvx_begin
    cvx_solver Mosek
    cvx_begin quiet
        variable delta_g_on(OPTIONS.N_g, 1) binary
        variable startup_g_on(OPTIONS.N_g, 1) binary
        variable redundant_switch_on(4, 1) binary
    %     variable shutdown_g(OPTIONS.N_g, 1) binary
%         variable switch_change(4, 1) binary
        variable benders_cut
        minimize( OPTIONS.G(1:OPTIONS.N_g, 3).' * delta_g_on(1:OPTIONS.N_g, 1) ...
                + OPTIONS.C_ss(1:OPTIONS.N_g, 1).' * startup_g_on(1:OPTIONS.N_g, 1) ...
                + benders_cut )

        subject to
            % startup detect
            if index_time == 1
                startup_g_on(1:OPTIONS.N_g, 1) >= delta_g_on(1:OPTIONS.N_g, 1) - zeros(OPTIONS.N_g, 1);
            else
                startup_g_on(1:OPTIONS.N_g, 1) >= delta_g_on(1:OPTIONS.N_g, 1) - master_delta_g(1:OPTIONS.N_g, index_time-1);
            end
        
            %  Redundent_switch
            redundant_switch_on(1, 1) + redundant_switch_on(2, 1) == 1;
            redundant_switch_on(3, 1) + redundant_switch_on(4, 1) == 1;

            % in the normal mode, redundant switches are determined.
            if operation_mode <= 3
                redundant_switch_on(1, 1) == total_sub(1).redundant_sw(1, index_time);
                redundant_switch_on(2, 1) == total_sub(1).redundant_sw(2, index_time);
                redundant_switch_on(3, 1) == total_sub(1).redundant_sw(3, index_time);
                redundant_switch_on(4, 1) == total_sub(1).redundant_sw(4, index_time);

                % benders cuts
                for index_benders = 1:Max_benders_iteration
                    benders_cut >=  total_sub(index_benders).sub_optval(2, index_time) ...
                                    + total_dual(index_benders).delta_g(1, index_time) ...
                                      * (delta_g_on(1, 1).' - total_sub(index_benders).delta_g(1, index_time).') ...
                                    + total_dual(index_benders).delta_g(2, index_time) ...
                                      * (delta_g_on(2, 1).' - total_sub(index_benders).delta_g(2, index_time).');
                end
            % in the fault mode, redundant switches need to be reconfigurated.
            elseif operation_mode <= 7
                % benders cuts
                for index_benders = 1:Max_benders_iteration
                    benders_cut >=  total_sub(index_benders).sub_optval(2, index_time) ...
                                    + total_dual(index_benders).delta_g(1, index_time) ...
                                      * (delta_g_on(1, 1).' - total_sub(index_benders).delta_g(1, index_time).') ...
                                    + total_dual(index_benders).delta_g(2, index_time) ...
                                      * (delta_g_on(2, 1).' - total_sub(index_benders).delta_g(2, index_time).') ... 
                                    + total_dual(index_benders).switch(1, index_time) ...
                                      * (redundant_switch_on(1, 1).' - total_sub(index_benders).redundant_sw(1, index_time).')...
                                    + total_dual(index_benders).switch(2, index_time) ...
                                      * (redundant_switch_on(2, 1).' - total_sub(index_benders).redundant_sw(2, index_time).')... 
                                    + total_dual(index_benders).switch(3, index_time) ...
                                      * (redundant_switch_on(3, 1).' - total_sub(index_benders).redundant_sw(3, index_time).')...
                                    + total_dual(index_benders).switch(4, index_time) ...
                                      * (redundant_switch_on(4, 1).' - total_sub(index_benders).redundant_sw(4, index_time).');
                end
            elseif operation_mode <= 11
            end

            switch accelerate_flag
                case 1
                case 2
                    % low bound of master problem
                    OPTIONS.G(1:OPTIONS.N_g, 3).' * delta_g_on(1:OPTIONS.N_g, 1) ...
                    + OPTIONS.C_ss(1:OPTIONS.N_g, 1).' * startup_g_on(1:OPTIONS.N_g, 1) ...
                    + benders_cut >= upper_of_lowerbound;
                    % speedup constraints: power range & lower bound
                    if operation_mode <=3 % normal mode
                        delta_g_on(1, 1)*OPTIONS.Pg_Max(1) + delta_g_on(2, 1)*OPTIONS.Pg_Max(2) ...
                            +  OPTIONS.Pb_Max >= OPTIONS.P_L_TIME_off(index_time) + Ppr_best_solution;
%                             +  OPTIONS.Pb_Max >= OPTIONS.P_L_TIME_off(index_time) + total_P(Max_benders_iteration).Ppr(index_time);
                        
                        delta_g_on(1, 1)*OPTIONS.Pg_Min(1) + delta_g_on(2, 1)*OPTIONS.Pg_Min(2) ...
                            +  OPTIONS.Pb_Min <= OPTIONS.P_total_vs(index_time) + 0; % Ppr
                    elseif operation_mode <= 7 % fault mode
                        if Max_benders_iteration >= 1
                            delta_g_on(1, 1)*OPTIONS.Pg_Max(1) + OPTIONS.Pb_Max ...
                                >= OPTIONS.island1_max(index_time) + Ppr_best_solution;
%                                 >= OPTIONS.island1_max(index_time) + total_P(Max_benders_iteration).Ppr(index_time);
                            delta_g_on(1, 1)*OPTIONS.Pg_Min(1) + OPTIONS.Pb_Min <= OPTIONS.island1_min(index_time) + 0; % Ppr
                            delta_g_on(2, 1)*OPTIONS.Pg_Max(2) + OPTIONS.Pb_Max >= OPTIONS.island2_max(index_time);
                            delta_g_on(2, 1)*OPTIONS.Pg_Min(2) + OPTIONS.Pb_Min <= OPTIONS.island2_min(index_time) + 0; % Ppr
                        end
                    end
                    % speedup constraint: lower bound
%                     benders_cut >= benders_cut_lowerbound;
                case 3
                    % speedup constraints: power range
                    if operation_mode <=3 % normal mode
                        delta_g_on(1, 1)*OPTIONS.Pg_Max(1) + delta_g_on(2, 1)*OPTIONS.Pg_Max(2) ...
                            +  OPTIONS.Pb_Max >= OPTIONS.P_L_TIME_off(index_time) + total_P(Max_benders_iteration).Ppr(index_time);
                        delta_g_on(1, 1)*OPTIONS.Pg_Min(1) + delta_g_on(2, 1)*OPTIONS.Pg_Min(2) ...
                            +  OPTIONS.Pb_Min <= OPTIONS.P_total_vs(index_time) + 0; % Ppr
                    elseif operation_mode <= 7 % fault mode
                        if Max_benders_iteration >= 1
                            delta_g_on(1, 1)*OPTIONS.Pg_Max(1) + OPTIONS.Pb_Max ...
                                >= OPTIONS.island1_max + total_P(Max_benders_iteration).Ppr;
                            delta_g_on(1, 1)*OPTIONS.Pg_Min(1) + OPTIONS.Pb_Min <= OPTIONS.island1_min(index_time) + 0; % Ppr
                            delta_g_on(2, 1)*OPTIONS.Pg_Max(2) + OPTIONS.Pb_Max >= OPTIONS.island2_max(index_time);
                            delta_g_on(2, 1)*OPTIONS.Pg_Min(2) + OPTIONS.Pb_Min <= OPTIONS.island2_min(index_time) + 0; % Ppr
                        end
                    end
                case 4
                    % speedup constraint: lower bound
                    benders_cut >= benders_cut_lowerbound;
    %             otherwise % non-constriant
            end

    cvx_end
    
    master_delta_g(1:OPTIONS.N_g, index_time) = delta_g_on(1:OPTIONS.N_g, 1);
    master_redundant_switch(1:4, index_time) = redundant_switch_on(1:4, 1);
    sub_optval_on = sub_optval_on + cvx_optval;
end

% the related description is given in the function of subproblem
if with_iteration_D == 0
    master_optval = sub_optval_on;
elseif with_iteration_D == 1
    master_optval = sub_optval_on + total_sub(Max_benders_iteration).sub_optval(2, OPTIONS.N_t);
end
        
if strcmp(cvx_status, 'Infeasible')
    disp('Master Infeasible');
    infeasible_flag = 1;
end
if strcmp(cvx_status, 'Failed')
    disp('Master Failed');
    infeasible_flag = 1;
end

% avoid the decimals affect the calculation of startup cost 
redundant_sw = round(master_redundant_switch);
delta_g = round(master_delta_g);
end