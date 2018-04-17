%% the main struture of benders decomposition algorithm
function [optimal_cost, time_val_problem ] ...
         = cost_optimization_for_test_benders( time_slot, voya_distance, accelerate_flag_input, near_opt_optimal_input, operation_mode_input, No_test_in, varphi )
global OPTIONS total_sub total_P total_dual near_opt_optimal accelerate_flag subobj_upperbound No_test parameter_test
% decomposition_flag
%% the adjusted parameters 
dbstop if error
if ~exist('time_slot', 'var')
    OPTIONS.N_t = 12;
else
    OPTIONS.N_t = time_slot;
end

if ~exist('voya_distance', 'var')
    OPTIONS.Distance = 150;
else
    OPTIONS.Distance = voya_distance;
end

if ~exist('accelerate_flag_input', 'var')
    accelerate_flag = 4;
else
    accelerate_flag = accelerate_flag_input;
end

if ~exist('near_opt_optimal_input', 'var')
    near_opt_optimal = 0;
else
    near_opt_optimal = near_opt_optimal_input;
end

if ~exist('parameter_test_input', 'var')
    parameter_test = 3;
else
    parameter_test = parameter_test_input;
end

% normal operation mode:
% 0 (normal wo PPA ESMC) 
% 1 (normal w PPA wo ESMC) 
% 2 (normal wo PPA w ESMC) 
% 3 (normal w PPA w ESMC)

% semi-island operation mode:
% 4 (Fault wo PPA ESMC)
% 5 (Fault w PPA wo ESMC)
% 6 (Fault wo PPA w ESMC)
% 7 (Fault w PPA w ESMC)

% island operation mode:
% 8 (Fault wo PPA ESMC)
% 9 (Fault w PPA wo ESMC)
% 10 (Fault wo PPA w ESMC)
% 11 (Fault w PPA w ESMC)
if ~exist('operation_mode_input', 'var')
    operation_mode = 7;
else
    operation_mode = operation_mode_input;
end

if ~exist('No_test_in', 'var')
    No_test = 0;
else
    No_test = No_test_in;
end

if ~exist('varphi', 'var')
    OPTIONS.varphi = 0.5;
else
    OPTIONS.varphi = varphi;
end

% decomposition_flag = 1;

%% initial the parameters
initial_parameters();
load_demand(operation_mode);

% preallocating an arry to the variables
Max_iteration = 800;
delta_g= zeros(OPTIONS.N_g, OPTIONS.N_t);
OPTIONS.initial_redundent_sw = [ones(1, OPTIONS.N_t); zeros(1, OPTIONS.N_t); zeros(1, OPTIONS.N_t); ones(1, OPTIONS.N_t)];
redundent_sw = OPTIONS.initial_redundent_sw ;
objval_upperbound = zeros(1, Max_iteration);
objval_lowerbound = zeros(1, Max_iteration);
best_upperbound = zeros(1, Max_iteration);
best_upperbound(1) = inf;
best_lowerbound = zeros(1, Max_iteration);
time_val_problem = zeros(3, Max_iteration);

delta_bound = [ones(1, OPTIONS.N_t); ones(1, OPTIONS.N_t); ];
% delta_bound = delta_g;
[objval_upperbound_test, benders_cut_lowerbound] = optimization_subproblem( operation_mode, delta_bound, redundent_sw );
benders_cut_lowerbound = 0 * benders_cut_lowerbound;
% benders_cut_lowerbound = 0;
subobj_upperbound = 0;

%% start the iteration of the algorithm
for benders_index = 1: 800
    tic;
    % sub_problem optimization based on benders decomposition
    [objval_upperbound(benders_index), sub_optval, sub_P, sub_dual, reduced_distance] ...
        = optimization_subproblem( operation_mode, delta_g, redundent_sw);
    
    % storage and passing of objective value and dual variable 
    total_sub(benders_index).sub_optval = sub_optval;
    total_sub(benders_index).delta_g = delta_g;
    total_sub(benders_index).redundent_sw = redundent_sw;
    
    total_P(benders_index).Pg = sub_P.Pg;
    total_P(benders_index).Pb = sub_P.Pb;
    total_P(benders_index).Ppr = sub_P.Ppr;
    total_P(benders_index).load_shedding = sub_P.load_shedding;
    total_P(benders_index).Pc = sub_P.Pc;
    total_P(benders_index).Pd = sub_P.Pd;
    
    total_dual(benders_index).delta_g = sub_dual.delta_g;
    total_dual(benders_index).switch = sub_dual.switch;
    
%     time_val_sub(benders_index) = toc;
    time_val_problem(1, benders_index) = toc;
    
    tic;
    % master_problem optimization based on benders decomposition
    [objval_lowerbound(benders_index), delta_g, redundent_sw, benders_cut ] ...
        = optimization_masterproblem( operation_mode, total_sub, total_dual, benders_cut_lowerbound );
    
    if objval_lowerbound(benders_index) > subobj_upperbound
        subobj_upperbound = objval_lowerbound(benders_index);
    end
%     if benders_cut > benders_cut_lowerbound
%         benders_cut_lowerbound = benders_cut;
%     end

%     time_val_master(benders_index) = toc;
    time_val_problem(2, benders_index) = toc;
    time_val_problem(3,benders_index) = sum(time_val_problem(1:2, benders_index), 1);
    fprintf('\nComputation time of iteration %d is %5.2f + %6.2f = %6.2f s \n\n', ...
        benders_index, time_val_problem(1, benders_index), time_val_problem(2, benders_index), time_val_problem(3,benders_index) );
    
    if benders_index >= 2
        % performance comparison
        error = ( objval_upperbound(benders_index) - objval_lowerbound(benders_index-1) );
%         error(benders_index) = ( objval_upperbound(benders_index) - objval_lowerbound(benders_index) );

        dual_gap = 100*error/objval_lowerbound(benders_index);
        disp('upperbound, lowerbound, error, dual_gap');
        disp([objval_upperbound(benders_index) objval_lowerbound(benders_index) error dual_gap]);
        
        if best_upperbound(benders_index-1) > objval_upperbound(benders_index)
            best_upperbound(benders_index) = objval_upperbound(benders_index);
        else
            best_upperbound(benders_index) = best_upperbound(benders_index-1);
        end

        if best_lowerbound(benders_index-1) < objval_lowerbound(benders_index)
            best_lowerbound(benders_index) = objval_lowerbound(benders_index);
        else
            best_lowerbound(benders_index) = best_lowerbound(benders_index-1);
        end
        
        % convergence determination
        if error <= 1e-3
            break;
        end
    end
    
%     if abs(objval_upperbound(benders_index) - 7.3724e+04) < 10
%         disp(1);
%     end

end

optimal_cost = [best_lowerbound(1:benders_index-1); best_upperbound(2:benders_index); ];
save_optimal_cost_information(optimal_cost, sub_P, delta_g, redundent_sw, reduced_distance, operation_mode );

plot_result(optimal_cost);

end

%% save optimal operation status
function save_optimal_cost_information(optimal_cost, sub_P, delta_g, redundent_sw, reduced_distance, operation_mode )
global OPTIONS  No_test near_opt_optimal accelerate_flag
operating_cost = OPTIONS.G(1:OPTIONS.N_g, 1).' * power(sub_P.Pg(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
               + OPTIONS.G(1:OPTIONS.N_g, 2).' * sub_P.Pg(1:OPTIONS.N_g, 1:OPTIONS.N_t) ...
               + OPTIONS.G(1:OPTIONS.N_g, 3).' * delta_g ...
               + OPTIONS.Xi_E * OPTIONS.E(1:OPTIONS.N_e, 1).' * power(sub_P.Pb(1:OPTIONS.N_e, 1:OPTIONS.N_t), 2) ...
               + OPTIONS.Xi_E * OPTIONS.E(1:OPTIONS.N_e, 2).' * ones(2, OPTIONS.N_t);
Penalty_cost_LS =  OPTIONS.Penalty_L * sum(sub_P.load_shedding(1:2, 1:OPTIONS.N_t),1);
startup =  delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) - [zeros(OPTIONS.N_g,1) delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t-1)];
startup_cost = OPTIONS.C_ss(1:OPTIONS.N_g,1).'* startup;
operating_cost = operating_cost + startup_cost + Penalty_cost_LS;

total_switching = sum( abs(redundent_sw(1,1) - OPTIONS.initial_redundent_sw(1,1)) ) + sum(abs(redundent_sw(3,1) - OPTIONS.initial_redundent_sw(3,1))) ...
        + sum( abs(redundent_sw(1,2:OPTIONS.N_t) - redundent_sw(1,1:OPTIONS.N_t-1)) ) + sum(abs(redundent_sw(3,2:OPTIONS.N_t) - redundent_sw(3,1:OPTIONS.N_t-1)));

data.power(1:OPTIONS.N_g, 1:OPTIONS.N_t) =  sub_P.Pg;
data.power(1:OPTIONS.N_e, 1:OPTIONS.N_t) =  sub_P.Pb;
data.power(5, 1:OPTIONS.N_t) =  sub_P.Ppr;
data.power(6:7, 1:OPTIONS.N_t) =  sub_P.load_shedding;
data.power(7, 1:OPTIONS.N_t) =  sum(sub_P.Pg + sub_P.Pb, 1) - sub_P.Ppr;
data.power(8:8+OPTIONS.N_g-1, 1:OPTIONS.N_t) =  sub_P.Pc;

data.status(1:OPTIONS.N_g, 1:OPTIONS.N_t) = delta_g;
data.status(OPTIONS.N_g+1:OPTIONS.N_g+2, 1:OPTIONS.N_t) = [redundent_sw(1,1:OPTIONS.N_t); redundent_sw(3,1:OPTIONS.N_t)];
data.status(OPTIONS.N_g+3, 1) = total_switching;

% data(6,:) = OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t);
data.cost(1,:) = operating_cost;
data.cost(2,1) = optimal_cost(1, end);

data.distance(1, 1:2) = [OPTIONS.Distance reduced_distance];

if near_opt_optimal ==0
    filename = ['optimal_data_mode_',num2str(operation_mode),'_D.',num2str(OPTIONS.Distance),...
                '_T.',num2str(OPTIONS.N_t),'_Ac.',num2str(accelerate_flag),'_No.',num2str(No_test),'.mat'];
elseif near_opt_optimal ==1
    filename = ['LNBD_data_mode_',num2str(operation_mode),'_D.',num2str(OPTIONS.Distance),...
                '_T.',num2str(OPTIONS.N_t),'_Ac.',num2str(accelerate_flag),'_No.',num2str(No_test),'.mat'];
end

save(filename,'data');

end

%% plot the simulation result
function plot_result(result)
hold on
plot(result(1,1:end));
plot(result(2,1:end));
hold off
end

%% load data generation
function load_demand(operation_mode)
global OPTIONS
OPTIONS.P_L = [2.7 0.9]; % P_Generater
OPTIONS.P_L_Scale_t = [0.5 0.6 0.8 0.9 0.82 0.6 0.4 0.35 0.25 0.33 0.4 0.5 0.4 0.3 0.6 0.8 0.82 0.9 0.9 0.7 0.62 0.5 0.33 0.4 0.5 0.4];
OPTIONS.P_vs_t_Scale = [0.5 0.6 0.8 0.9 0.82 0.6 0.4 0.35 0.25 0.33 0.4 0.5 0.4 0.3 0.6 0.8 0.82 0.9 0.9 0.7 0.62 0.5 0.33 0.4 0.5 0.4];
OPTIONS.P_no_t_Scale = [0.5 0.6 0.8 0.9 0.82 0.6 0.4 0.35 0.25 0.33 0.4 0.5 0.4 0.3 0.6 0.8 0.82 0.9 0.9 0.7 0.62 0.5 0.33 0.4 0.5 0.4];

% total load upper bound 3.6 MW
P_no_base = 1.2; % the scale of non-vital loads are 1.2/3.6 = 0.33
P_vs_base = ones(1,OPTIONS.zone).' * (2.7 + 0.9 - 1.2)/OPTIONS.zone;

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
    % not the real max, it's the maximum of the supply power of vital and
    % semi loads, which can be considered as the minimum value of upperbound 
    OPTIONS.island1_max = OPTIONS.island1_load + sum(P_coupled_load(1:2, 1:OPTIONS.N_t), 1);
    OPTIONS.island1_min = P_vs_island1(1:OPTIONS.N_t);
    
    OPTIONS.island2_load = P_vs_island2(1:OPTIONS.N_t) + P_no_island2(1:OPTIONS.N_t);
    OPTIONS.island2_load_average = sum(OPTIONS.island2_load)/OPTIONS.N_t;
    OPTIONS.island2_non_load = P_no_island2(1:OPTIONS.N_t);
    % not the real max, it's the maximum of the supply power of vital and semi loads
    OPTIONS.island2_max = OPTIONS.island2_load + sum(P_coupled_load(1:2, 1:OPTIONS.N_t), 1);
    OPTIONS.island2_min = P_vs_island2(1:OPTIONS.N_t);
    
    P_vs = P_vs_base * OPTIONS.P_vs_t_Scale;
    P_no = P_no_base * OPTIONS.P_no_t_Scale;
    P_total_vs = sum(P_vs);
    
    P_total_time = P_total_vs + P_no;
    OPTIONS.P_L_TIME_off = P_total_time(1:OPTIONS.N_t);
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
    % not the real max, it's the maximum of the supply power of vital and
    % semi loads, which can be considered as the minimum value of upperbound 
    OPTIONS.island1_max = OPTIONS.island1_load;
    OPTIONS.island1_min = P_vs_island1(1:OPTIONS.N_t);
    
    OPTIONS.island2_load = P_vs_island2(1:OPTIONS.N_t) + P_no_island2(1:OPTIONS.N_t);
    OPTIONS.island2_load_average = sum(OPTIONS.island2_load)/OPTIONS.N_t;
    OPTIONS.island2_non_load = P_no_island2(1:OPTIONS.N_t);
    % not the real max, it's the maximum of the supply power of vital and semi loads
    OPTIONS.island2_max = OPTIONS.island2_load;
    OPTIONS.island2_min = P_vs_island2(1:OPTIONS.N_t);
    
    P_vs = P_vs_base * OPTIONS.P_vs_t_Scale;
    P_no = P_no_base * OPTIONS.P_no_t_Scale;
    P_total_vs = sum(P_vs);
    
    P_total_time = P_total_vs + P_no;
    OPTIONS.P_L_TIME_off = P_total_time(1:OPTIONS.N_t);
end

end

%% inital the parameters for optimization
function [] = initial_parameters()
global OPTIONS
OPTIONS.velocity = [25 0];
% number of ESM modules
OPTIONS.N_e = 2;
% number of generators
OPTIONS.N_g = 2;
% number of time slots
% OPTIONS.N_t = 12;
% minmum time between startup and shuntdown
OPTIONS.Tmin_g = 2;
% minmum time between switch change
OPTIONS.Tmin_sw = 2;
% time of faults happen
OPTIONS.Fault_time = 4;
% the number of zones in shipboard power system
OPTIONS.zone = 6;

OPTIONS.velocity_avg = OPTIONS.Distance/OPTIONS.N_t;
OPTIONS.P_pr_avg = (OPTIONS.velocity_avg).^3*2.2e-3;
OPTIONS.Pg_Max(1) = 8;
OPTIONS.Pg_Min(1) = 1;
OPTIONS.Pg_Max(2) = 4;
OPTIONS.Pg_Min(2) = 0.0;

OPTIONS.Xi_E = 1;
OPTIONS.Ppr_Max = (OPTIONS.velocity(1)).^3*2.2e-3;
OPTIONS.Pb_Max(1) = 1;
OPTIONS.Pb_Min(1) = -1;
OPTIONS.E_Max = [2 2];
OPTIONS.E_Min = [0.2 0.2];

% generator function parameters
OPTIONS.G(1,1:3) = [10 30 220];
OPTIONS.G(2,1:3) = [14.5 12 170];
OPTIONS.E(1:2,1:2) = [1 0.5; 1 0.5];
OPTIONS.C_ss = [90; 45];
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

% load_information = [OPTIONS.P_L_TIME_off; OPTIONS.P_L_TIME_on; OPTIONS.Coupled_load];
% save('load_information');
end

%% the subproblem which is used to calculate the optimal power of generators and ESMs
function [objval_upperbound, sub_optval, sub_P, sub_dual, reduced_distance] = optimization_subproblem( operation_mode, delta_g, redundent_sw )
global OPTIONS
%% optimal algorithm
%  cvx_begin
% cvx_solver SeDuMi
cvx_begin quiet
    variable Ppr(1,OPTIONS.N_t) nonnegative
    variable Pb(2,OPTIONS.N_t)
    variable E(2,OPTIONS.N_t) nonnegative
    variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    variable redundent_sw_s(4,OPTIONS.N_t)
    variable sub_delta_g(OPTIONS.N_g,OPTIONS.N_t)
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
            + sum(10 * OPTIONS.Penalty_D * Pc(1, 1:OPTIONS.N_t).' ,1) ...
            + sum(10 * OPTIONS.Penalty_D * Pc(2, 1:OPTIONS.N_t).' ,1) ...
            + sum(10 * OPTIONS.Penalty_D * Pd(1, 1:OPTIONS.N_t).' ,1) ...
            + sum(10 * OPTIONS.Penalty_D * Pd(2, 1:OPTIONS.N_t).' ,1) )

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
        Pb(1, 1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1, OPTIONS.N_t)
        Pb(1, 1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1, OPTIONS.N_t)
        Pb(2, 1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1, OPTIONS.N_t)
        Pb(2, 1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1, OPTIONS.N_t)

        % charging and discharging
        E(1, 1:OPTIONS.N_t) <= OPTIONS.E_Max(1) * ones(1, OPTIONS.N_t)
        E(1, 1:OPTIONS.N_t) >= OPTIONS.E_Min(1) * ones(1, OPTIONS.N_t)
        E(2, 1:OPTIONS.N_t) <= OPTIONS.E_Max(2) * ones(1, OPTIONS.N_t)
        E(2, 1:OPTIONS.N_t) >= OPTIONS.E_Min(2) * ones(1, OPTIONS.N_t)

        % ESM output power and the capacity constraints
        E(1, 1) == OPTIONS.E_Max(1) - Pb(1,1)
        E(2, 1) == OPTIONS.E_Max(2) - Pb(2,1)
        for t_index = 2:OPTIONS.N_t
            E(1, t_index) == E(1, t_index-1) - Pb(1, t_index)
            E(2, t_index) == E(2, t_index-1) - Pb(2, t_index)
        end

        % system power balance & load shedding range
        if operation_mode <= 3 % normal operation without load shedding
            for t_index = 1:OPTIONS.N_t
                 OPTIONS.P_L_TIME_off(1, t_index) - sum(load_shedding(1:2, t_index)) + Ppr(1, t_index) == sum( Pg(1:OPTIONS.N_g, t_index) ) + sum(Pb(1:OPTIONS.N_e, t_index) )
            end

            % load shedding range
            load_shedding(1, 1:OPTIONS.N_t) == zeros(1, OPTIONS.N_t)
            load_shedding(2, 1:OPTIONS.N_t) == zeros(1, OPTIONS.N_t)
        elseif operation_mode <= 7 % semi-island mode operation with load shedding and reconfiguration
            for t_index = 1:OPTIONS.N_t
                redundent_sw_s(1:2, t_index).'*OPTIONS.Coupled_load(:, t_index) + OPTIONS.island1_load(1, t_index) - load_shedding(1, t_index) + Ppr(1, t_index) == (Pg(1, t_index)) + Pb(1, t_index) 
                redundent_sw_s(3:4, t_index).'*OPTIONS.Coupled_load(:, t_index) + OPTIONS.island2_load(1, t_index) - load_shedding(2, t_index) == (Pg(2, t_index)) + Pb(2, t_index)
            end

            % load shedding amount is limited by the demand of non-vital loads
            load_shedding(1, 1:OPTIONS.N_t) <= OPTIONS.island1_non_load(1, 1:OPTIONS.N_t)
            load_shedding(2, 1:OPTIONS.N_t) <= OPTIONS.island2_non_load(1, 1:OPTIONS.N_t)
        elseif operation_mode <= 11 % island fault mode operation with load shedding and reconfiguration
            for t_index = 1:OPTIONS.N_t
                OPTIONS.island1_load(1, t_index) - load_shedding(1, t_index) + Ppr(1, t_index) == (Pg(1, t_index)) + Pb(1, t_index) 
                OPTIONS.island2_load(1, t_index) - load_shedding(2, t_index) == (Pg(2, t_index)) + Pb(2, t_index)
            end

            % load shedding amount is limited by the demand of non-vital loads
            load_shedding(1, 1:OPTIONS.N_t) <= OPTIONS.island1_non_load(1, 1:OPTIONS.N_t)
            load_shedding(2, 1:OPTIONS.N_t) <= OPTIONS.island2_non_load(1, 1:OPTIONS.N_t)
        end

        dual_Sp : redundent_sw_s(1:2, :) == redundent_sw(1:2, :)
        dual_Ss : redundent_sw_s(3:4, :) == redundent_sw(3:4, :)

        % voyage planning
        switch mod(operation_mode, 4)
            case 0  % only generator scheduling
                Ppr(1, 1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
                Pb(1:2, 1:OPTIONS.N_t) == 0;
            case 1  % generator scheduling & PPA
                sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance;
                Pb(1:2,1:OPTIONS.N_t) == 0;
            case 2  % generator scheduling & ESMC
                Ppr(1, 1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
            case 3  % generator scheduling & ESMC & PPA
                sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance; 
%             case 4  % only generator scheduling
%                 Ppr(1, 1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
%                 Pb(1:2, 1:OPTIONS.N_t) == 0;
%             case 5  % generator scheduling & PPA
%                 sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance;
%                 Pb(1:2, 1:OPTIONS.N_t) == 0;
%             case 6  % generator scheduling & ESMC
%                 Ppr(1, 1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
%             case 7  % generator scheduling & ESMC & PPA
%                 sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance; 
        end

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
end

%% startup cost: first startup cost and other startup cost
startup = ( (sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) - [zeros(OPTIONS.N_g, 1) sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t-1)] ) >=1);
% shutdown = ( (sub_delta_g(1:OPTIONS.N_g,1:OPTIONS.N_t) - [ sub_delta_g(1:OPTIONS.N_g, 2:OPTIONS.N_t) zeros(OPTIONS.N_g,1) ] ) >=1);
startup_cost = sum(OPTIONS.C_ss(1:OPTIONS.N_g, 1).'* startup ,2);
% shutdown_cost = OPTIONS.weight(1)*sum(OPTIONS.C_ss(1:OPTIONS.N_g,1).'* shutdown ,2);
one_item_cost = sum(OPTIONS.G(1:OPTIONS.N_g, 3).' * sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) ,2) ;

objval_upperbound = sub_optval + startup_cost + one_item_cost;

sub_dual.switch(1:2, :) = dual_Sp;
sub_dual.switch(3:4, :) = dual_Ss;
sub_dual.delta_g(1, :) = dual_delta_g1;
sub_dual.delta_g(2, :) = dual_delta_g2;

if objval_upperbound == 5.5164e+04
    disp(objval_upperbound);
end

% disp('upperbound');
end


%% the master problem which is used to determine the redundent switches and ??
function [cvx_optval, master_delta_g, master_redundent_switch, benders_cut ] = optimization_masterproblem( operation_mode, total_sub, total_dual, benders_cut_lowerbound )
global OPTIONS accelerate_flag  total_P subobj_upperbound

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
    variable startup(OPTIONS.N_g, OPTIONS.N_t) binary
    variable master_redundent_switch(4, OPTIONS.N_t) binary
    variable switch_change(4, OPTIONS.N_t-1) binary
    variable benders_cut
    minimize( sum(OPTIONS.G(1:OPTIONS.N_g, 3).' * master_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
            + sum(OPTIONS.C_ss(1:OPTIONS.N_g, 1).' * startup(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
            + benders_cut )

    subject to
        % startup detect
        startup(1:OPTIONS.N_g, 1) >= master_delta_g(1:OPTIONS.N_g, 1) -0
        for t_index = 2 : OPTIONS.N_t
            startup(1:OPTIONS.N_g, t_index) >= master_delta_g(1:OPTIONS.N_g, t_index) - master_delta_g(1:OPTIONS.N_g, t_index-1)
        end
        
        %  generator switch time
        sum(master_delta_g(1:OPTIONS.N_g, 1:OPTIONS.Tmin_g), 2) >= OPTIONS.Tmin_g * startup(1:OPTIONS.N_g, 1)
        for t_index = 2:1:OPTIONS.N_t - OPTIONS.Tmin_g+1
            sum(master_delta_g(1:OPTIONS.N_g, t_index:t_index + OPTIONS.Tmin_g-1), 2) >= OPTIONS.Tmin_g * startup(1:OPTIONS.N_g, t_index)
        end

        %  Redundent_switch
        master_redundent_switch(1, :) + master_redundent_switch(2, :) == ones(1, OPTIONS.N_t)
        master_redundent_switch(3, :) + master_redundent_switch(4, :) == ones(1, OPTIONS.N_t)
        
        if operation_mode < 3
            master_redundent_switch(1, :) == total_sub(1).redundent_sw(1, :)
            master_redundent_switch(2, :) == total_sub(1).redundent_sw(3, :)
            master_redundent_switch(3, :) == total_sub(1).redundent_sw(3, :)
            master_redundent_switch(4, :) == total_sub(1).redundent_sw(4, :)
        end
        
        % switch change detect
        switch_change(1, 1:OPTIONS.N_t-1) >= (master_redundent_switch(1, 2:OPTIONS.N_t) - master_redundent_switch(1, 1:OPTIONS.N_t-1))
        switch_change(2, 1:OPTIONS.N_t-1) >= (master_redundent_switch(2, 2:OPTIONS.N_t) - master_redundent_switch(2, 1:OPTIONS.N_t-1))
        switch_change(3, 1:OPTIONS.N_t-1) >= (master_redundent_switch(3, 2:OPTIONS.N_t) - master_redundent_switch(3, 1:OPTIONS.N_t-1))
        switch_change(4, 1:OPTIONS.N_t-1) >= (master_redundent_switch(4, 2:OPTIONS.N_t) - master_redundent_switch(4, 1:OPTIONS.N_t-1))
        
        %  Redundent switch time
        sum(master_redundent_switch(1, 1:OPTIONS.Tmin_sw), 2) >= OPTIONS.Tmin_sw * switch_change(1:OPTIONS.N_g, 1)
        sum(master_redundent_switch(2, 1:OPTIONS.Tmin_sw), 2) >= OPTIONS.Tmin_sw * switch_change(2:OPTIONS.N_g, 1)
        for t_index = 2 : OPTIONS.N_t-OPTIONS.Tmin_sw+1
            sum(master_redundent_switch(1, t_index:t_index + OPTIONS.Tmin_sw-1), 2) >= OPTIONS.Tmin_sw * switch_change(1, t_index)
            sum(master_redundent_switch(2, t_index:t_index + OPTIONS.Tmin_sw-1), 2) >= OPTIONS.Tmin_sw * switch_change(2, t_index)
        end

         % benders cuts
%         for index_benders = max([size(total_sub,2)-2 1] ):size(total_sub,2)
        for index_benders = 1:size(total_sub, 2)
            benders_cut >=  total_sub(index_benders).sub_optval ...
                            + total_dual(index_benders).delta_g(1, :)*(master_delta_g(1, :).' - total_sub(index_benders).delta_g(1, :).')...
                            + total_dual(index_benders).delta_g(2, :)*(master_delta_g(2, :).' - total_sub(index_benders).delta_g(2, :).')... 
                            + total_dual(index_benders).switch(1, :)*(master_redundent_switch(1, :).' - total_sub(index_benders).redundent_sw(1, :).')...
                            + total_dual(index_benders).switch(2, :)*((master_redundent_switch(2, :)).' - total_sub(index_benders).redundent_sw(2, :).')... 
                            + total_dual(index_benders).switch(3, :)*(master_redundent_switch(3, :).' - total_sub(index_benders).redundent_sw(3, :).')...
                            + total_dual(index_benders).switch(4, :)*((master_redundent_switch(4, :)).' - total_sub(index_benders).redundent_sw(4, :).');
%                             + total_dual(index_benders).delta_g(1,:)*master_delta_g(1,:).' - total_dual(index_benders).delta_g(1,:)*total_sub(index_benders).delta_g(1,:).'...
%                             + total_dual(index_benders).delta_g(2,:)*master_delta_g(2,:).' - total_dual(index_benders).delta_g(2,:)*total_sub(index_benders).delta_g(2,:).' ...
%                             + total_dual(index_benders).switch(1,:)*master_redundent_switch(1,:).' - total_dual(index_benders).switch(1,:)*total_sub(index_benders).redundent_sw(1,:).'...
%                             + total_dual(index_benders).switch(2,:)*master_redundent_switch(2,:).' - total_dual(index_benders).switch(2,:)*total_sub(index_benders).redundent_sw(2,:).';
        end
 
        switch accelerate_flag
            case 1 
                % only power range
                % speedup constraints: power range
                if operation_mode <= 3 % normal mode
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) +  OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.P_L_TIME_off + total_P(end).Ppr
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) +  OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.P_total_vs + 0 % Ppr
                elseif operation_mode <= 7 % fault mode
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.island1_max + total_P(end).Ppr
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island1_min + 0 % Ppr
                    master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.island2_max
                    master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island2_min + 0 % Ppr
                end    
            case 2
                % only lower bound
                % speedup constraint: lower bound
                benders_cut >= benders_cut_lowerbound
            case 3 
                % have two constraints
                % speedup constraints: power range
                if operation_mode <=3 % normal mode
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) +  OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.P_L_TIME_off + total_P(end).Ppr
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) +  OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.P_total_vs + 0 % Ppr
%                 % speedup constraint: lower bound
%                 benders_cut >= benders_cut_lowerbound
                elseif operation_mode <= 7 % fault mode
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.island1_max + total_P(end).Ppr
                    master_delta_g(1, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island1_min + 0 % Ppr
                    master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) + OPTIONS.Pb_Max * ones(1, OPTIONS.N_t) >= OPTIONS.island2_max
                    master_delta_g(2, 1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) + OPTIONS.Pb_Min * ones(1, OPTIONS.N_t) <= OPTIONS.island2_min + 0 % Ppr
                end

                % speedup constraint: lower bound
                benders_cut >= benders_cut_lowerbound
            otherwise % non-constriant
        end
        
%         sum(OPTIONS.G(1,3) * master_delta_g(1,1:OPTIONS.N_t) ,2) ...
%         + sum(OPTIONS.G(2,3) * master_delta_g(2,1:OPTIONS.N_t) ,2) ...
%         + sum(OPTIONS.C_ss(1:2,1).'* startup(1:2, 1:OPTIONS.N_t) ,2) ...
%         + benders_cut <= subobj_upperbound
cvx_end
        
if strcmp(cvx_status, 'Infeasible')
    disp('stop');
end

redundent_sw = master_redundent_switch;

end


%% the subproblem t which is used to calculate the optimal power of generators and ESMs at time t  
function [objval_upperbound, sub_optval, sub_P, sub_dual, reduced_distance] = optimization_subproblem_t( operation_mode, delta_g, redundent_sw )
global OPTIONS near_opt_optimal

%% near_opt algorithm
% inital the rest parameters
Rest_Distance = OPTIONS.Distance;
Rest_ppr_avg = OPTIONS.P_pr_avg;
sub_optval_on = 0;

% update the rest parameters: rest distance, rest average ppr, and the
% power deviation of load at next time slot. 
for index_time = 1:OPTIONS.N_t
    switch index_time
        case OPTIONS.N_t
            varphi_sub = 1;
        otherwise
            varphi_sub = OPTIONS.varphi;
    end

    if operation_mode <= 3
        Delta_PL =  OPTIONS.P_L_TIME_off(index_time) - OPTIONS.P_L_TIME_off_avg ;
        Distance_slot_obj = ((Rest_ppr_avg - (1 - varphi_sub)*Delta_PL)/2.2e-3)^(1/3);

        % ESM output power and the capacity constraints
        if index_time == 1
            rest_pmax_ESM = OPTIONS.E_Max(1) + OPTIONS.E_Max(2) - OPTIONS.E_Min(1) - OPTIONS.E_Min(2);
            rest_pmin_ESM = OPTIONS.E_Max(1) + OPTIONS.E_Max(2)  -OPTIONS.E_Max(1) - OPTIONS.E_Max(2);
        else
            rest_pmax_ESM = E(1,index_time-1) + E(2,index_time-1) - OPTIONS.E_Min(1) - OPTIONS.E_Min(2) ;
            rest_pmin_ESM = E(1,index_time-1) + E(2,index_time-1) - OPTIONS.E_Max(1) - OPTIONS.E_Max(2) ;
        end
        upper_bound_ESM_P = roundn( min( varphi_sub * Delta_PL, rest_pmax_ESM ), -4);
        upper_bound_ESM_N = roundn( max( varphi_sub * Delta_PL, rest_pmin_ESM ), -4);
    elseif operation_mode <= 7
        Delta_PL_island1 = OPTIONS.island1_load(index_time) - OPTIONS.island1_load_average;
        Delta_PL_island2 = OPTIONS.island2_load(index_time) - OPTIONS.island2_load_average;
        Distance_slot_obj = ((Rest_ppr_avg - (1 - varphi_sub)*Delta_PL_island1)/2.2e-3)^(1/3);

        % ESM output power and the capacity constraints
        if index_time == 1
            rest_pmax_ESM1 = OPTIONS.E_Max(1) - OPTIONS.E_Min(1) ;
            rest_pmax_ESM2 = OPTIONS.E_Max(2) - OPTIONS.E_Min(2) ;
            rest_pmin_ESM1 = OPTIONS.E_Max(1) - OPTIONS.E_Max(1) ;
            rest_pmin_ESM2 = OPTIONS.E_Max(2) - OPTIONS.E_Max(2) ;
        else
            rest_pmax_ESM1 = E(1,index_time-1) - OPTIONS.E_Min(1) ;
            rest_pmax_ESM2 = E(2,index_time-1) - OPTIONS.E_Min(2) ;
            rest_pmin_ESM1 = E(1,index_time-1) - OPTIONS.E_Max(1) ;
            rest_pmin_ESM2 = E(2,index_time-1) - OPTIONS.E_Max(2) ;
        end
        if sum(delta_g(1,index_time))>=1
            upper_bound_ESM1P = roundn( min( varphi_sub * Delta_PL_island1, rest_pmax_ESM1 ), -4);
            upper_bound_ESM1N = roundn( max( varphi_sub * Delta_PL_island1, rest_pmin_ESM1 ), -4);
%                 upper_bound_ESM2P = roundn( min( varphi_sub * Delta_PL_island2, rest_pmax_ESM2 ), -4);
%                 upper_bound_ESM2N = roundn( max( varphi_sub * Delta_PL_island2, rest_pmin_ESM2 ), -4);
        else
            upper_bound_ESM1P = roundn( rest_pmax_ESM1, -4);
            upper_bound_ESM1N = roundn( rest_pmin_ESM1, -4);
%                 upper_bound_ESM2P = roundn( min( varphi_sub * Delta_PL_island2, rest_pmax_ESM2 ), -4);
%                 upper_bound_ESM2N = roundn( max( varphi_sub * Delta_PL_island2, rest_pmin_ESM2 ), -4);
        end
        if sum(delta_g(2,index_time))>=1
%                 upper_bound_ESM1P = roundn( min( varphi_sub * Delta_PL_island1, rest_pmax_ESM1 ), -4);
%                 upper_bound_ESM1N = roundn( max( varphi_sub * Delta_PL_island1, rest_pmin_ESM1 ), -4);
            upper_bound_ESM2P = roundn( min( varphi_sub * Delta_PL_island2, rest_pmax_ESM2 ), -4);
            upper_bound_ESM2N = roundn( max( varphi_sub * Delta_PL_island2, rest_pmin_ESM2 ), -4);
        else
%                 upper_bound_ESM1P = roundn( min( varphi_sub * Delta_PL_island1, rest_pmax_ESM1 ), -4);
%                 upper_bound_ESM1N = roundn( max( varphi_sub * Delta_PL_island1, rest_pmin_ESM1 ), -4);
            upper_bound_ESM2P = roundn( rest_pmax_ESM2, -4);
            upper_bound_ESM2N = roundn( rest_pmin_ESM2, -4);
        end
    end       

    % cvx_begin
    % cvx_solver SeDuMi
    cvx_begin quiet
        variable Pg_on(OPTIONS.N_g) nonnegative
        variable Pb_on(2,1)
        variable E_on(2,1) nonnegative
        variable Ppr_on(1) nonnegative
        variable sub_delta_g(OPTIONS.N_g)
        variable redundent_sw_s(4,1)
        variable load_shedding_on(2, 1) nonnegative
        variable reduced_distance_on nonnegative
        variable Pc(OPTIONS.N_g) nonnegative
        variable Pd(OPTIONS.N_g) nonnegative
        dual variable temp_dual_delta_g1
        dual variable temp_dual_delta_g2
        dual variable temp_dual_Sp
        dual variable temp_dual_Ss
        minimize( sum(OPTIONS.G(1,1) * power(Pg_on(1,1).',2) ,1) ...
                + sum(OPTIONS.G(2,1) * power(Pg_on(2,1).',2) ,1) ...
                + sum(OPTIONS.G(1,2) * power(Pg_on(1,1).',1) ,1) ...
                + sum(OPTIONS.G(2,2) * power(Pg_on(2,1).',1) ,1) ...
                + OPTIONS.Xi_E * sum(OPTIONS.E(1:2,1).'* power(Pb_on(1:OPTIONS.N_e,1),2) ,2) ...
                + OPTIONS.Xi_E * sum([1 1] * ones(OPTIONS.N_e,1) ,2)...
                + sum(OPTIONS.Penalty_L * sum(load_shedding_on(1:2,1),2) ,1) ...
                + 1.02 * OPTIONS.Penalty_D * reduced_distance_on ...
                + sum(10 * OPTIONS.Penalty_D * Pc(1).' ,1) ...
                + sum(10 * OPTIONS.Penalty_D * Pc(2).' ,1) ...
                + sum(10 * OPTIONS.Penalty_D * Pd(1).' ,1) ...
                + sum(10 * OPTIONS.Penalty_D * Pd(2).' ,1))

        subject to
            % the range constraints of all the variables
            Pg_on(1,1) <= sub_delta_g(1,1) * OPTIONS.Pg_Max(1) + Pc(1)
            Pg_on(1,1) >= sub_delta_g(1,1) * OPTIONS.Pg_Min(1) - Pd(1)
            Pg_on(2,1) <= sub_delta_g(2,1) * OPTIONS.Pg_Max(2) + Pc(2)
            Pg_on(2,1) >= sub_delta_g(2,1) * OPTIONS.Pg_Min(2) - Pd(2)

            temp_dual_delta_g1 : sub_delta_g(1,1) == delta_g(1,index_time)
            temp_dual_delta_g2 : sub_delta_g(2,1) == delta_g(2,index_time)

            % ramping rate power of generators
            switch index_time
                case 1
                    Pg_on(1,1) <= OPTIONS.R_G(1) 
                    Pg_on(2,1) <= OPTIONS.R_G(2) 
%                     case OPTIONS.N_t
%                         Pg_on(1,1) <= OPTIONS.R_G(1) 
%                         Pg_on(2,1) <= OPTIONS.R_G(2)
                otherwise
                    Pg_on(1,1) - Pg(1,index_time-1) <= OPTIONS.R_G(1) 
                    Pg_on(1,1) - Pg(1,index_time-1) >= -OPTIONS.R_G(1) 
                    Pg_on(2,1) - Pg(2,index_time-1) <= OPTIONS.R_G(2) 
                    Pg_on(2,1) - Pg(2,index_time-1) >= -OPTIONS.R_G(2) 
            end

            % propulsion power limitation
            Ppr_on(1,1) <= OPTIONS.Ppr_Max

            % ESM limitation constraints
            if operation_mode <= 3
                if Delta_PL >= 0
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
                    Pb_on(1,1) + Pb_on(2,1) <= 2 * OPTIONS.Pb_Max
                    Pb_on(1,1) + Pb_on(2,1) <= upper_bound_ESM_P
                    Pb_on(1,1) + Pb_on(2,1) >= 0
                elseif Delta_PL < 0
                    Pb_on(1,1) + Pb_on(2,1) <= 0
                    Pb_on(1,1) + Pb_on(2,1) <= upper_bound_ESM_N
                    Pb_on(1,1) + Pb_on(2,1) >= 2 * OPTIONS.Pb_Min
                end
            elseif operation_mode <= 7
                if Delta_PL_island1 >=0
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
                    Pb_on(1,1) <= OPTIONS.Pb_Max
%                         Pb_on(1,1) <= upper_bound_ESM1P
                    Pb_on(1,1) >= 0
                elseif Delta_PL_island1 <0
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
                    Pb_on(1,1) <= 0
%                         Pb_on(1,1) <= upper_bound_ESM1N
                    Pb_on(1,1) >= OPTIONS.Pb_Min 
                end

                if Delta_PL_island2 >=0 && delta_g(2,index_time) == 0 
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
                    Pb_on(2,1) <= OPTIONS.Pb_Max
%                         Pb_on(2,1) <= upper_bound_ESM2P
                    Pb_on(2,1) >= 0
                elseif Delta_PL_island2 >=0 && delta_g(2,index_time) == 1 
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
%                         Pb_on(2,1) <= 0
                    Pb_on(2,1) <= OPTIONS.Pb_Max
%                         Pb_on(2,1) <= upper_bound_ESM2N
                    Pb_on(2,1) >= OPTIONS.Pb_Min
                elseif Delta_PL_island2 <=0 && delta_g(2,index_time) == 0 
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
                    Pb_on(2,1) <= 0
                    Pb_on(2,1) >= OPTIONS.Pb_Min
                elseif Delta_PL_island2 <=0 && delta_g(2,index_time) == 1 
                    % the scaler factor is related with load power in each
                    % island part and the adjusting factor.
%                         Pb_on(2,1) <= 0
                    Pb_on(2,1) <= OPTIONS.Pb_Max
%                         Pb_on(2,1) <= upper_bound_ESM2N
                    Pb_on(2,1) >= OPTIONS.Pb_Min
                end
            end

            % charging and discharging
            E_on(1,1) <=  OPTIONS.E_Max(1)
            E_on(2,1) <=  OPTIONS.E_Max(2)
            E_on(1,1) >= 0
            E_on(2,1) >= 0

            % ESM output power and the capacity constraints
            if index_time == 1
                OPTIONS.E_Max(1) - Pb_on(1,1) == E_on(1,1)
                OPTIONS.E_Max(2) - Pb_on(2,1) == E_on(2,1)
            else
                E(1,index_time-1) - Pb_on(1,1) == E_on(1,1)
                E(2,index_time-1) - Pb_on(2,1) == E_on(2,1)
            end

            % system power balance & load shedding
            if operation_mode <= 3
                 OPTIONS.P_L_TIME_off(1,index_time) - sum(load_shedding_on(1,1)) + Ppr_on(1,1) == sum( Pg_on(1:OPTIONS.N_g,1) ) + sum(Pb_on(1:OPTIONS.N_e,1))
                 load_shedding_on(2,1) == 0

                % load shedding range
                load_shedding_on(1) == zeros(1,1)
                load_shedding_on(2) == zeros(1,1)

            elseif operation_mode <= 7
                for t_index = index_time
                    redundent_sw_s(1:2,1).'*OPTIONS.Coupled_load(:,index_time) + OPTIONS.island1_load(1, t_index) - load_shedding_on(1,1) + Ppr_on(1,1) == (Pg_on(1,1)) + Pb_on(1,1) 
                    redundent_sw_s(3:4,1).'*OPTIONS.Coupled_load(:,index_time) + OPTIONS.island2_load(1, t_index) - load_shedding_on(2,1) == (Pg_on(2,1)) + Pb_on(2,1)
                end

                % load shedding range
                load_shedding_on(1) <= OPTIONS.island1_non_load(1,index_time)
                load_shedding_on(2) <= OPTIONS.island2_non_load(1,index_time)
            end

            temp_dual_Sp : redundent_sw_s(1:2,1) == redundent_sw(1:2,index_time)
            temp_dual_Ss : redundent_sw_s(3:4,1) == redundent_sw(3:4,index_time)

            % voyage planning
            switch operation_mode
                case 0
                    Ppr_on(1,1) == OPTIONS.P_pr_avg;
                    Pb_on(1:2,1) == 0;
                case 1
                    (Ppr_on(1,1)./2.2e-3).^(1/3) >= Distance_slot_obj - reduced_distance_on;
                    Pb_on(1:2,1) == 0;
                case 2
                    Ppr_on(1,1) == OPTIONS.P_pr_avg;
                case 3
                    (Ppr_on(1,1)./2.2e-3).^(1/3) >= Distance_slot_obj - reduced_distance_on; 
                case 4
                    Ppr_on(1,1) == OPTIONS.P_pr_avg;
                    Pb(1:2,1) == 0;
                case 5
                    sum((Ppr_on(1,1)./2.2e-3).^(1/3)) >= Distance_slot_obj - reduced_distance_on;
                    Pb_on(1:2,1) == 0;
                case 6
                    Ppr_on(1,1) == OPTIONS.P_pr_avg;
                case 7
                    sum((Ppr_on(1,1)./2.2e-3).^(1/3)) >= Distance_slot_obj - reduced_distance_on; 
            end

    cvx_end

    if strcmp(cvx_status, 'Infeasible')
        disp('stop');
    end

    Pg(1,index_time) = Pg_on(1,1);
    Pg(2,index_time) = Pg_on(2,1);
    Pb(1,index_time) = roundn(Pb_on(1,1), -4);
    Pb(2,index_time) = roundn(Pb_on(2,1), -4);
    Ppr(1,index_time) = Ppr_on(1,1);
    E(1,index_time) = E_on(1,1);
    E(2,index_time) = E_on(2,1);
    load_shedding(1:2,index_time) = load_shedding_on(1:2,1);

    Rest_Distance =  Rest_Distance - Distance_slot_obj + reduced_distance_on;
    Rest_velocity_avg = Rest_Distance/(OPTIONS.N_t-index_time);
    Rest_ppr_avg = (Rest_velocity_avg).^3*2.2e-3;

    sub_optval_on = sub_optval_on + cvx_optval;

    dual_Sp(1:2,index_time) = temp_dual_Sp;
    dual_Ss(1:2,index_time) = temp_dual_Ss;
    dual_delta_g1(1,index_time) = temp_dual_delta_g1;
    dual_delta_g2(1,index_time) = temp_dual_delta_g2;
end
%     sub_delta_g = ones(2,OPTIONS.N_t);
reduced_distance = Rest_Distance;
sub_delta_g(1:2,1:OPTIONS.N_t) = delta_g(1:2,1:OPTIONS.N_t);
sub_optval = sub_optval_on;

%% startup cost: first startup cost and other startup cost
startup = ( (sub_delta_g(1:OPTIONS.N_g,1:OPTIONS.N_t) - [zeros(OPTIONS.N_g,1) sub_delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t-1)] ) >=1);
% shutdown = ( (sub_delta_g(1:OPTIONS.N_g,1:OPTIONS.N_t) - [ sub_delta_g(1:OPTIONS.N_g, 2:OPTIONS.N_t) zeros(OPTIONS.N_g,1) ] ) >=1);
startup_cost = sum(OPTIONS.C_ss(1:OPTIONS.N_g,1).'* startup ,2);
% shutdown_cost = OPTIONS.weight(1)*sum(OPTIONS.C_ss(1:OPTIONS.N_g,1).'* shutdown ,2);
one_item_cost = sum(OPTIONS.G(1:OPTIONS.N_g,3).' * sub_delta_g(1:OPTIONS.N_g,1:OPTIONS.N_t) ,2) ...
              + sum(OPTIONS.E(1:OPTIONS.N_e,2) * sub_load_shedding(1:OPTIONS.N_e,1:OPTIONS.N_t) ,2) ;

objval_upperbound = sub_optval + startup_cost + one_item_cost;

sub_dual.switch(1:2,:) = dual_Sp;
sub_dual.switch(3:4,:) = dual_Ss;
sub_dual.delta_g(1,:) = dual_delta_g1;
sub_dual.delta_g(2,:) = dual_delta_g2;

if objval_upperbound == 5.5164e+04
    disp(objval_upperbound);
end

end