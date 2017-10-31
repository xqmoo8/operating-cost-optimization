function [optimal_cost, time_val_problem ] = cost_optimization_for_test_benders( time_slot, voya_distance, accelerate_flag_input, online_offline_input, operation_mode_input )
global OPTIONS online_offline accelerate_flag decomposition_flag
%% 
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
    accelerate_flag = 3;
else
    accelerate_flag = accelerate_flag_input;
end

if ~exist('online_offline_input', 'var')
    online_offline = 1;
else
    online_offline = online_offline_input;
end

%%  operation_mode:
% 0 (normal wo PPA ESMC) 
% 1 (normal w PPA wo ESMC) 
% 2 (normal wo PPA w ESMC) 
% 3 (normal w PPA w ESMC)
% 4 (Fault wo PPA ESMC)
% 5 (Fault w PPA wo ESMC) 
% 6 (Fault wo PPA w ESMC) 
% 7 (Fault w PPA w ESMC)
if ~exist('operation_mode_input', 'var')
    operation_mode = 3;
else
    operation_mode = operation_mode_input;
end

%% initial the parameters
initial_parameters();
load_demand(operation_mode);

% startup variable used in master problem
Pg = zeros(2,OPTIONS.N_t);
Pb = zeros(2,OPTIONS.N_t);
delta = zeros(2, OPTIONS.N_t);
% delta = ones(2, OPTIONS.N_t);
redundent_sw = [ones(1, OPTIONS.N_t); zeros(1, OPTIONS.N_t); zeros(1, OPTIONS.N_t); ones(1, OPTIONS.N_t)];
objval_upperbound = 0;
objval_lowerbound = 0;
best_upperbound = inf;
best_lowerbound = 0;
benders_cut_lowerbound = 0;
error = 0;
dual_gap = 0;

decomposition_flag = 1;

%%
% % test code for benders decompostion
% delta = ones(2, OPTIONS.N_t);
% [objval_upperbound, sub_optval, dual_delta, dual_switch, Pg, Pb, Ppr ] = optimization_subproblem( operation_mode, delta, redundent_sw, Pg, Pb );
% total_sub(1).sub_optval = sub_optval;
% total_sub(1).delta = delta;
% total_sub(1).redundent_sw = redundent_sw;
% total_dual(1).delta = dual_delta;
% total_dual(1).switch = dual_switch;
% [objval_lowerbound(1), delta, redundent_sw ] = optimization_masterproblem( operation_mode, total_sub, total_dual );

% % benchmark for testing benders optimization
% delta = ones(2, OPTIONS.N_t);
% [result_benchmark(1), Pg, Pb, Ppr ] = total_cost_optimization( operation_mode, delta, redundent_sw, Pg, Pb );
% for startup_index = 1: 12+30
%     delta = starup_slot(startup_index);
%     [result_benchmark(startup_index+1), Pg, Pb, Ppr ] = total_cost_optimization( operation_mode, delta, redundent_sw, Pg, Pb );
% end

%% start the iteration of the algorithm
for benders_index = 1: 800
    tic;
    % sub_problem optimization based on benders decomposition
    [objval_upperbound(benders_index), sub_optval, dual_delta, dual_switch, Pg, Pb, Ppr ] = optimization_subproblem( operation_mode, delta, redundent_sw );
    % storage and passing of objective value and dual variable 
    total_sub(benders_index).sub_optval = sub_optval;
    total_sub(benders_index).delta = delta;
    total_sub(benders_index).redundent_sw = redundent_sw;
    total_dual(benders_index).delta = dual_delta;
    total_dual(benders_index).switch = dual_switch;
    
    time_val_sub(benders_index) = toc;
    tic;
    % master_problem optimization based on benders decomposition
    [objval_lowerbound(benders_index), delta, redundent_sw, benders_cut ] = optimization_masterproblem( operation_mode, total_sub, total_dual, benders_cut_lowerbound, Ppr );
    if benders_cut > benders_cut_lowerbound
        benders_cut_lowerbound = benders_cut;
    end

    time_val_master(benders_index) = toc;
    time_val_problem(1:3,benders_index) = [time_val_sub(benders_index); time_val_master(benders_index); time_val_sub(benders_index)+time_val_master(benders_index) ];
    
    if benders_index >= 2
        % performance comparison
        error(benders_index) = objval_upperbound(benders_index) - objval_lowerbound(benders_index-1);
        dual_gap(benders_index) = 100*error(benders_index)/objval_lowerbound(benders_index);
        disp('upperbound, lowerbound, error, dual_gap');
        disp([objval_upperbound(benders_index) objval_lowerbound(benders_index) error(benders_index) dual_gap(benders_index)]);
        
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
        if error(benders_index) <= 1e-3
            break;
        end
    end

end

optimal_cost = [objval_upperbound(2:end); objval_lowerbound(1:end-1); error(2:end); dual_gap(2:end); time_val_sub(2:end); time_val_master(2:end)];
% result = [best_upperbound(2:end); best_lowerbound(1:end-1); error(2:end); dual_gap(2:end);];
save('optimal_cost.mat','optimal_cost');

plot_result(optimal_cost);

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
P_vs_base = ones(1,OPTIONS.zone).'*2.4/OPTIONS.zone;
P_no_base = 1.2;

if operation_mode <= 3
    P_vs = P_vs_base * OPTIONS.P_vs_t_Scale;
    P_no = P_no_base * OPTIONS.P_no_t_Scale;

    P_total_time = sum(P_vs) + P_no;
    P_total_average = sum(P_total_time(1:OPTIONS.N_t))/OPTIONS.N_t;
    
    OPTIONS.P_L_TIME_off = P_total_time(1:OPTIONS.N_t);
    OPTIONS.P_L_TIME_off_avg = P_total_average;
elseif operation_mode <= 7
    % the default fault is at 5-6 bus and 7-8bus. Thus,the scale factors of
    % vital and semi-vital loads in different island parts are 3/6 and 1/6;
    % The scale factors of non-vital loads are (5+9-6)/12 and 1-(5+9-6)/12
    % (0.75 and 0.25)
    P_vs_island1 = sum(P_vs_base(1:3)) * OPTIONS.P_vs_t_Scale;
    P_vs_island2 = sum(P_vs_base(6)) * OPTIONS.P_vs_t_Scale;
    P_no_island1 = 0.75*P_no_base * OPTIONS.P_no_t_Scale;
    P_no_island2 = 0.25*P_no_base * OPTIONS.P_no_t_Scale;
    P_coupled_load = (P_vs_base(4:5)) * OPTIONS.P_vs_t_Scale;
    
    OPTIONS.Coupled_load = P_coupled_load;
    OPTIONS.island1_load = P_vs_island1 + P_no_island1;
    OPTIONS.island1_load_average = OPTIONS.island1_load(1:OPTIONS.N_t)/OPTIONS.N_t;
    OPTIONS.island1_non_load = P_no_island1;
    OPTIONS.island2_load = P_vs_island2 + P_no_island2;
    OPTIONS.island2_load_average = OPTIONS.island2_load(1:OPTIONS.N_t)/OPTIONS.N_t;
    OPTIONS.island2_non_load = P_no_island2;
end

% % the load demand without random feature
% P_L_Scale_on = interp1(1:13,OPTIONS.P_L_Scale_t(1:13),1:0.5:13,'spline') + 0.1*rand(1, 25); 
% % the load demand with random feature
% % P_L_Scale_off = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline');
% P_L_Scale_off = interp1(1:13,OPTIONS.P_L_Scale_t(1:13),1:0.5:13,'spline');
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
OPTIONS.Pg_Min(2) = 0.5;

OPTIONS.Xi_E = 1;
OPTIONS.Ppr_Max = (OPTIONS.velocity(1)).^3*2.2e-3;
OPTIONS.Pb_Max(1) = 1;
OPTIONS.Pb_Min(1) = -1;
OPTIONS.E_Max = [2 2];

% generator function parameters
OPTIONS.G(1,1:3) = [10 30 220];
OPTIONS.G(2,1:3) = [14.5 12 170];
OPTIONS.E(1:2,1:2) = [1 5; 1 5];
OPTIONS.C_ss = [90; 45];
OPTIONS.R_G = 1;
OPTIONS.error = 1e-3;
OPTIONS.Penalty = max([2*OPTIONS.G(1,1)*OPTIONS.Pg_Max(1) + OPTIONS.G(1,2) 2*OPTIONS.G(2,1)*OPTIONS.Pg_Max(2) + OPTIONS.G(2,2)]);
A_max = max(OPTIONS.G(1:2,1));
B_max = max(OPTIONS.G(1:2,2));
P_G_max = max(OPTIONS.Pg_Max);
P_E_max = max(OPTIONS.Pb_Max);

OPTIONS.Penalty_L = max(2*A_max*P_G_max + B_max, 2*OPTIONS.Xi_E*OPTIONS.E(1,1)*OPTIONS.Pb_Max);
OPTIONS.Penalty_D = 3*(2.2e-3)^(1/3)*OPTIONS.Penalty_L/(OPTIONS.Ppr_Max)^(1/3-1);

OPTIONS.alpha = 2/3;

% load_information = [OPTIONS.P_L_TIME_off; OPTIONS.P_L_TIME_on; OPTIONS.Coupled_load];
% save('load_information');
end

%% heuristic algorithm
function [delta] = starup_slot(startup_index)
global OPTIONS
delta = zeros(2, OPTIONS.N_t);
% one slot shutdown for two generators
if startup_index <= 12
    delta(startup_index) = 1;
% two slots shutdown for two generators
elseif startup_index <= 42
    index = startup_index - 12;
    delta(1,ceil(index/5)) = 0;
    if mod(index-1,5)+1 < ceil(index/5)
        delta(2,mod(index-1,5)+1) = 0;
    elseif mod(index-1,5)+1 > ceil(index/5)
        delta(2,mod(index-1,5)+2) = 0;
    elseif mod(index-1,5)+1 == ceil(index/5)
        delta(2,mod(index-1,5)+2) = 0;
    end
% % three slots shutdown for two generators
% elseif startup_index <= 122
%     index = startup_index - 42;
%     delta(1,ceil(index/5)) = 0;
%     if mod(index,5)<ceil(index/5)
%         delta(2,mod(index,5)) = 0;
%     elseif mod(index,5)>=ceil(index/5)
%         delta(2,mod(index,5)+1) = 0;
%     end
end
% disp(delta);
% start_index = OPTIONS.N_t - size(sequence,2);
end

%% the subproblem which is used to calculate the optimal power of generators and ESMs
function [ cvx_optval, Pg, Pb, Ppr ] = total_cost_optimization( operation_mode, delta, redundent_sw,  Pg_m, Pb )  
global OPTIONS

startup = (delta(1:2,2:end) - delta(1:2,1:end-1)>=1);

cvx_begin
% cvx_solver SeDuMi
cvx_begin quiet
    variable Ppr(1,OPTIONS.N_t) nonnegative
    variable Pb(2,OPTIONS.N_t)
    variable E(2,OPTIONS.N_t) nonnegative
    variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    minimize( sum(OPTIONS.G(1,1) * delta(1,1:OPTIONS.N_t) * power(Pg(1,1:OPTIONS.N_t).',2) ,2) ...
            + sum(OPTIONS.G(2,1) * delta(2,1:OPTIONS.N_t) * power(Pg(2,1:OPTIONS.N_t).',2) ,2) ...
            + sum(OPTIONS.G(1,2) * delta(1,1:OPTIONS.N_t) * power(Pg(1,1:OPTIONS.N_t).',1) ,2) ...
            + sum(OPTIONS.G(2,2) * delta(2,1:OPTIONS.N_t) * power(Pg(2,1:OPTIONS.N_t).',1) ,2) ...
            + sum(OPTIONS.G(1,3) * delta(1,1:OPTIONS.N_t) ,2) ...
            + sum(OPTIONS.G(2,3) * delta(2,1:OPTIONS.N_t) ,2) ...
            + OPTIONS.Xi * sum(OPTIONS.E(1:2,1).'* power(Pb(1:OPTIONS.N_e,1:OPTIONS.N_t),2) ,2) ...
            + OPTIONS.Xi * sum(OPTIONS.E(1:2,2).'* ones(OPTIONS.N_g,OPTIONS.N_t) ,2) ...
            + sum(OPTIONS.C_ss(1:2,1).'* startup ,2) )
 
    subject to
        % the range constraints of all the variables
        Pg(1,1:OPTIONS.N_t) <= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1)
        Pg(2,1:OPTIONS.N_t) <= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Max(2)
        Pg(1,1:OPTIONS.N_t) >= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)
        Pg(2,1:OPTIONS.N_t) >= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(2)

        % ramping rate power of generators
        Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) <= OPTIONS.R_G
        Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) >= -OPTIONS.R_G
        Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) <= OPTIONS.R_G
        Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) >= -OPTIONS.R_G

        % propulsion power limitation
        Ppr(1,1:OPTIONS.N_t) <= OPTIONS.Ppr_Max * ones(1,OPTIONS.N_t)
        
        % ESM limitation
        Pb(1,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
        Pb(2,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
        Pb(1,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
        Pb(2,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
        
        % charging and discharging
        E(1,1:OPTIONS.N_t) <=  OPTIONS.E_Max(1) * ones(1,OPTIONS.N_t)
        E(2,1:OPTIONS.N_t) <=  OPTIONS.E_Max(2) * ones(1,OPTIONS.N_t)
        E(1,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)
        E(2,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)
        
        % ESM output power and the capacity constraints
        OPTIONS.E_Max(1) - Pb(1,1) == E(1,1)
        OPTIONS.E_Max(2) - Pb(2,1) == E(2,1)
        for t_index = 1:OPTIONS.N_t-1
            E(1,t_index) - Pb(1,t_index+1) == E(1,t_index+1)
            E(2,t_index) - Pb(2,t_index+1) == E(2,t_index+1)
        end
        
        % system power balance
        if operation_mode <= 3 
            for t_index = 1:OPTIONS.N_t
                 OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == sum( Pg(1:OPTIONS.N_g,t_index) ) + sum(Pb(1:OPTIONS.N_e,t_index))
            end
        elseif operation_mode <= 7
            for t_index = 1:OPTIONS.N_t
                redundent_sw(1,:)*OPTIONS.Coupled_load(:,t_index) +  (1-OPTIONS.alpha - 2/6) * OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == (Pg(1,t_index)) + Pb(1,t_index) 
                redundent_sw(2,:)*OPTIONS.Coupled_load(:,t_index) + OPTIONS.alpha * OPTIONS.P_L_TIME_off(1,t_index)  == (Pg(2,t_index)) + Pb(2,t_index)
            end
        end
        
        % voyage planning            
        if operation_mode ==0
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
            Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==1
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance
        elseif operation_mode ==2
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
        elseif operation_mode ==3
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance 
        elseif operation_mode ==4
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
            Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==5
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance
            Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==6
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
        else
            sum((Ppr(5:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - 4*( OPTIONS.P_pr_avg ./2.2e-3).^(1/3)
        end
cvx_end

disp('optimal');
disp(cvx_optval);

end

%% the subproblem which is used to calculate the optimal power of generators and ESMs
function [objval_upperbound, sub_optval, dual_delta, dual_switch, Pg, Pb, Ppr ] = optimization_subproblem( operation_mode, delta, redundent_sw )  
global OPTIONS online_offline

lowerbound(1,:) = delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1);
lowerbound(2,:) = delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(2);

%% offline algorithm
if online_offline ==0
    for index = 1:1
%         cvx_begin
        % cvx_solver SeDuMi
        cvx_begin quiet
            variable Ppr(1,OPTIONS.N_t) nonnegative
            variable Pb(2,OPTIONS.N_t)
            variable E(2,OPTIONS.N_t) nonnegative
            variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
            variable redundent_sw_s(4,OPTIONS.N_t)
            variable delta_s(OPTIONS.N_g,OPTIONS.N_t)
            variable load_shedding(2, OPTIONS.N_t) nonnegative
            variable reduced_distance nonnegative
            variable Pc(OPTIONS.N_g, OPTIONS.N_t) nonnegative
            variable Pd(OPTIONS.N_g, OPTIONS.N_t) nonnegative
            dual variable dual_delta_g1
            dual variable dual_delta_g2
            dual variable dual_Sp
            dual variable dual_Ss
            minimize( sum(OPTIONS.G(1,1) * power(Pg(1,1:OPTIONS.N_t).',2) ,1) ...
                    + sum(OPTIONS.G(2,1) * power(Pg(2,1:OPTIONS.N_t).',2) ,1) ...
                    + sum(OPTIONS.G(1,2) * power(Pg(1,1:OPTIONS.N_t).',1) ,1) ...
                    + sum(OPTIONS.G(2,2) * power(Pg(2,1:OPTIONS.N_t).',1) ,1) ...
                    + OPTIONS.Xi_E * sum(OPTIONS.E(1:2,1).'* power(Pb(1:OPTIONS.N_e,1:OPTIONS.N_t),2) ,2) ...
                    + OPTIONS.Xi_E * sum(OPTIONS.E(1:2,2).'* ones(OPTIONS.N_g,OPTIONS.N_t) ,2)...
                    + sum(OPTIONS.Penalty_L * sum(load_shedding(1:2,1:OPTIONS.N_t),2) ,1) ...
                    + OPTIONS.Penalty_D * reduced_distance ...
                    + sum(OPTIONS.Penalty * Pc(1,1:OPTIONS.N_t).' ,1) ...
                    + sum(OPTIONS.Penalty * Pc(2,1:OPTIONS.N_t).' ,1) ...
                    + sum(OPTIONS.Penalty * Pd(1,1:OPTIONS.N_t).' ,1) ...
                    + sum(OPTIONS.Penalty * Pd(2,1:OPTIONS.N_t).' ,1) )

            subject to
                % the range constraints of all the variables
                Pg(1,1:OPTIONS.N_t) <= delta_s(1,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1) + Pc(1,1:OPTIONS.N_t)
                Pg(1,1:OPTIONS.N_t) >= delta_s(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1) - Pd(1,1:OPTIONS.N_t)
                Pg(2,1:OPTIONS.N_t) <= delta_s(2,1:OPTIONS.N_t) * OPTIONS.Pg_Max(2) + Pc(2,1:OPTIONS.N_t)
                Pg(2,1:OPTIONS.N_t) >= delta_s(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(2) - Pd(2,1:OPTIONS.N_t)

                dual_delta_g1 : delta_s(1,1:OPTIONS.N_t) == delta(1,1:OPTIONS.N_t)
                dual_delta_g2 : delta_s(2,1:OPTIONS.N_t) == delta(2,1:OPTIONS.N_t)

                % ramping rate power of generators
                Pg(1,1) <= OPTIONS.Pg_Max(1)
                Pg(2,1) <= OPTIONS.Pg_Max(2)
                Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) <= OPTIONS.R_G
                Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) >= -OPTIONS.R_G
                Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) <= OPTIONS.R_G
                Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) >= -OPTIONS.R_G
                Pg(1,OPTIONS.N_t) <= OPTIONS.R_G
                Pg(2,OPTIONS.N_t) <= OPTIONS.R_G

                % propulsion power limitation
                Ppr(1,1:OPTIONS.N_t) <= OPTIONS.Ppr_Max * ones(1,OPTIONS.N_t)

                % ESM limitation
                Pb(1,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
                Pb(2,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
                Pb(1,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
                Pb(2,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)

                % charging and discharging
                E(1,1:OPTIONS.N_t) <=  OPTIONS.E_Max(1) * ones(1,OPTIONS.N_t)
                E(2,1:OPTIONS.N_t) <=  OPTIONS.E_Max(2) * ones(1,OPTIONS.N_t)
                E(1,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)
                E(2,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)

                % ESM output power and the capacity constraints
                OPTIONS.E_Max(1) - Pb(1,1) == E(1,1)
                OPTIONS.E_Max(2) - Pb(2,1) == E(2,1)
                for t_index = 1:OPTIONS.N_t-2
                    E(1,t_index) - Pb(1,t_index+1) == E(1,t_index+1)
                    E(2,t_index) - Pb(2,t_index+1) == E(2,t_index+1)
                end

                % system power balance
                if operation_mode <= 3
                    for t_index = 1:OPTIONS.N_t
                         OPTIONS.P_L_TIME_off(1,t_index) - sum(load_shedding(1:2,t_index)) + Ppr(1,t_index) == sum( Pg(1:OPTIONS.N_g,t_index) ) + sum(Pb(1:OPTIONS.N_e,t_index))
                    end
                elseif operation_mode <= 7
                    for t_index = 1:OPTIONS.N_t
                        redundent_sw_s(1:2, t_index).'*OPTIONS.Coupled_load(:, t_index) + OPTIONS.island1_load(1, t_index) - load_shedding(1,t_index) + Ppr(1, t_index) == (Pg(1, t_index)) + Pb(1, t_index) 
                        redundent_sw_s(3:4, t_index).'*OPTIONS.Coupled_load(:, t_index) + OPTIONS.island2_load(1, t_index) - load_shedding(2,t_index) == (Pg(2, t_index)) + Pb(2, t_index)
                    end
                end

                dual_Sp : redundent_sw_s(1:2, :) == redundent_sw(1:2, :)
                dual_Ss : redundent_sw_s(3:4, :) == redundent_sw(3:4, :)

                % voyage planning
                switch operation_mode
                    case 0
                        Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
                        Pb(1:2,1:OPTIONS.N_t) == 0;
%                         break;
                    case 1
                        sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance;
                        Pb(1:2,1:OPTIONS.N_t) == 0;
                    case 2
                        Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
                    case 3
                        sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance; 
                    case 4
                        Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
                        Pb(1:2,1:OPTIONS.N_t) == 0;
                    case 5
                        sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance;
                        Pb(1:2,1:OPTIONS.N_t) == 0;
                    case 6
                        Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
                    case 7
                        sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - reduced_distance; 
                end
                
        cvx_end
        sub_optval = cvx_optval;
    end

%% online algorithm
elseif online_offline == 1
    % inital the rest parameters
    Rest_Distance = OPTIONS.Distance;
    Rest_ppr_avg = OPTIONS.P_pr_avg;
    sub_optval_on = 0;
    
    % update the rest parameters: rest distance, rest average ppr, and the
    % power deviation of load at next time slot. 
    for index_time = 1:OPTIONS.N_t
        if operation_mode <= 3
            Delta_PL =  OPTIONS.P_L_TIME_off(index_time) - OPTIONS.P_L_TIME_off_avg ;
            Distance_slot_obj = ((Rest_ppr_avg - (1 - OPTIONS.alpha)*Delta_PL)/2.2e-3)^(1/3);
        elseif operation_mode <= 7
            Delta_PL_island1 = OPTIONS.island1_load(index_time) - OPTIONS.island1_load_average;
            Delta_PL_island2 = OPTIONS.island2_load(index_time) - OPTIONS.island2_load_average;
            Distance_slot_obj = ((Rest_ppr_avg - (1 - OPTIONS.alpha)*Delta_PL_island1)/2.2e-3)^(1/3);
        end
        
        % cvx_begin
        % cvx_solver SeDuMi
        cvx_begin quiet
            variable Pg_on(OPTIONS.N_g) nonnegative
            variable Pb_on(2,1)
            variable E_on(2,1) nonnegative
            variable Ppr_on(1) nonnegative
            variable delta_s(OPTIONS.N_g)
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
                    + OPTIONS.Xi_E * sum(OPTIONS.E(1:2,2).'* ones(OPTIONS.N_g,1) ,2)...
                    + sum(OPTIONS.Penalty_L * sum(load_shedding_on(1:2,1),2) ,1) ...
                    + OPTIONS.Penalty_D * reduced_distance_on ...
                    + sum(OPTIONS.Penalty * Pc(1).' ,1) ...
                    + sum(OPTIONS.Penalty * Pc(2).' ,1) ...
                    + sum(OPTIONS.Penalty * Pd(1).' ,1) ...
                    + sum(OPTIONS.Penalty * Pd(2).' ,1))

            subject to
                % the range constraints of all the variables
                Pg_on(1,1) <= delta_s(1,1) * OPTIONS.Pg_Max(1) + Pc(1)
                Pg_on(1,1) >= delta_s(1,1) * OPTIONS.Pg_Min(1) - Pd(1)
                Pg_on(2,1) <= delta_s(2,1) * OPTIONS.Pg_Max(2) + Pc(2)
                Pg_on(2,1) >= delta_s(2,1) * OPTIONS.Pg_Min(2) - Pd(2)

                temp_dual_delta_g1 : delta_s(1,1) == delta(1,1)
                temp_dual_delta_g2 : delta_s(2,1) == delta(2,1)

                % ramping rate power of generators
                if index_time == 1
                    Pg_on(1,1) <= OPTIONS.R_G + Pc(1)
                    Pg_on(2,1) <= OPTIONS.R_G + Pc(2)
                else
                    Pg_on(1,1) - Pg(1,index_time-1) <= OPTIONS.R_G + Pc(1)
                    Pg_on(1,1) - Pg(1,index_time-1) >= -OPTIONS.R_G - Pd(1)
                    Pg_on(2,1) - Pg(2,index_time-1) <= OPTIONS.R_G + Pc(2)
                    Pg_on(2,1) - Pg(2,index_time-1) >= -OPTIONS.R_G - Pd(2)
                end

                % propulsion power limitation
                Ppr_on(1,1) <= OPTIONS.Ppr_Max

                % ESM limitation constraints (38)
                if operation_mode <= 3
                    if Delta_PL >= 0
                        % the scaler factor is related with load power in each
                        % island part and the adjusting factor.
                        Pb_on(1,1) + Pb_on(2,1) <= 2 * OPTIONS.Pb_Max
                        Pb_on(1,1) + Pb_on(2,1) <= OPTIONS.alpha * (OPTIONS.P_L_TIME_off(index_time) - OPTIONS.P_L_TIME_off_avg)
                        Pb_on(1,1) + Pb_on(2,1) >= 0
                    elseif Delta_PL < 0
                        Pb_on(1,1) + Pb_on(2,1) <= 0
                        Pb_on(1,1) + Pb_on(2,1) >= OPTIONS.alpha * (OPTIONS.P_L_TIME_off(index_time) - OPTIONS.P_L_TIME_off_avg)
                        Pb_on(1,1) + Pb_on(2,1) >= 2 * OPTIONS.Pb_Min
                    end
                elseif operation_mode <= 7
                    if Delta_PL_island1 >=0
                        % the scaler factor is related with load power in each
                        % island part and the adjusting factor.
                        Pb_on(1,1) <= OPTIONS.Pb_Max
                        Pb_on(1,1) <= OPTIONS.alpha*(OPTIONS.island1_load - OPTIONS.island1_load_average)
                        Pb_on(1,1) >= 0
                    elseif Delta_PL_island1 <0
                        % the scaler factor is related with load power in each
                        % island part and the adjusting factor.
                        Pb_on(1,1) <= 0
                        Pb_on(1,1) >= OPTIONS.alpha*(OPTIONS.island1_load - OPTIONS.island1_load_average)
                        Pb_on(1,1) >= OPTIONS.Pb_Min 
                    end
                    
                    if Delta_PL_island2 >=0
                        % the scaler factor is related with load power in each
                        % island part and the adjusting factor.
                        Pb_on(2,1) <= OPTIONS.Pb_Max
                        Pb_on(2,1) <= OPTIONS.alpha*(OPTIONS.island2_load - OPTIONS.island2_load_average)
                        Pb_on(2,1) >= 0
                    elseif Delta_PL_island2 <0
                        % the scaler factor is related with load power in each
                        % island part and the adjusting factor.
                        Pb_on(2,1) <= 0
                        Pb_on(2,1) >= OPTIONS.alpha*(OPTIONS.island2_load - OPTIONS.island2_load_average)
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

                % system power balance
                if operation_mode <= 3
                     OPTIONS.P_L_TIME_off(1,index_time) - sum(load_shedding_on(1,1)) + Ppr_on(1,1) == sum( Pg_on(1:OPTIONS.N_g,1) ) + sum(Pb_on(1:OPTIONS.N_e,1))
                     load_shedding_on(2,1) == 0
                elseif operation_mode <= 7
                    for t_index = 1
                        redundent_sw_s(1,:)*OPTIONS.Coupled_load(:,index_time) + OPTIONS.island1_load(1, t_index) - load_shedding_on(1,1) + Ppr_on(1,1) == (Pg_on(1,1)) + Pb_on(1,1) 
                        redundent_sw_s(2,:)*OPTIONS.Coupled_load(:,index_time) + OPTIONS.island2_load(1, t_index) - load_shedding_on(2,1) == (Pg_on(2,1)) + Pb_on(2,1)
                    end
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
        
        Pg(1,index_time) = Pg_on(1,1);
        Pg(2,index_time) = Pg_on(2,1);
        Pb(1,index_time) = Pb_on(1,1);
        Pb(2,index_time) = Pb_on(2,1);
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
%     delta_s = ones(2,OPTIONS.N_t);
    delta_s(1:2,1:OPTIONS.N_t) = delta(1:2,1:OPTIONS.N_t);
    sub_optval = sub_optval_on;
end

%% startup cost: first startup cost and other startup cost
startup = (delta_s(1:2,2:OPTIONS.N_t) - delta_s(1:2,1:OPTIONS.N_t-1)>=1);
startup_cost = sum(OPTIONS.C_ss(1:2,1).'* startup ,2) + sum(OPTIONS.C_ss(1:2,1));
one_item_cost = sum(OPTIONS.G(1,3) * delta_s(1,1:OPTIONS.N_t) ,2) + sum(OPTIONS.G(2,3) * delta_s(2,1:OPTIONS.N_t) ,2);

objval_upperbound = sub_optval + startup_cost + one_item_cost;

dual_switch(1:2,:) = dual_Sp;
dual_switch(3:4,:) = dual_Ss;
dual_delta(1,:) = dual_delta_g1;
dual_delta(2,:) = dual_delta_g2;

% disp('upperbound');
% disp(objval_upperbound);

end

%% the master problem which is used to determine the redundent switches and ??
function [cvx_optval, master_delta, master_redundent_switch, benders_cut ] = optimization_masterproblem( operation_mode, total_sub, total_dual, benders_cut_lowerbound, Ppr )
global OPTIONS accelerate_flag

if ~exist('benders_cut', 'var')
    benders_cut = 0;
end

if ~exist('Pd', 'var')
    Pd = 0;
end

% %% dual variable and upperbound
% lambda_B(1,:) = -2*OPTIONS.G(1,1)*Pg(1,:) - OPTIONS.G(1,2)*Pg(1,:);
% lambda_B(2,:) = -2*OPTIONS.G(2,1)*Pg(2,:) - OPTIONS.G(2,2)*Pg(2,:);
% % lambda_Pg(1,1:OPTIONS.N_t) = - 2*OPTIONS.G(1,1)*Pg(1,:) - lambda_B(1,:);
% % lambda_Pg(2,1:OPTIONS.N_t) = - 2*OPTIONS.G(1,1)*Pg(2,:) - lambda_B(2,:);
% lambda_delta(1,:) = -OPTIONS.G(1,1)*power(Pg(1,:),2) - OPTIONS.G(1,2)*Pg(1,:);
% lambda_delta(2,:) = -OPTIONS.G(2,1)*power(Pg(2,:),2) - OPTIONS.G(2,2)*Pg(2,:);
% for index_t = 1:OPTIONS.N_t
%     lambda_Px(1, index_t) = lambda_B(1, index_t).'*OPTIONS.Coupled_load(1, index_t);
%     lambda_Px(2, index_t) = lambda_B(1, index_t).'*OPTIONS.Coupled_load(2, index_t);
%     lambda_Sy(1, index_t) = lambda_B(2, index_t).'*OPTIONS.Coupled_load(1, index_t);
%     lambda_Sy(2, index_t) = lambda_B(2, index_t).'*OPTIONS.Coupled_load(2, index_t);
% end

% cvx_begin
cvx_solver Mosek
cvx_begin quiet
    variable master_delta(OPTIONS.N_g, OPTIONS.N_t) binary
    variable master_redundent_switch(4, OPTIONS.N_t) binary
    variable startup(OPTIONS.N_g, OPTIONS.N_t-1) binary
    variable switch_change(4, OPTIONS.N_t-1) binary
    variable benders_cut 
%     variable E(2,OPTIONS.N_t) nonnegative
%     variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    minimize( sum(OPTIONS.G(1,3) * master_delta(1,1:OPTIONS.N_t) ,2) ...
            + sum(OPTIONS.G(2,3) * master_delta(2,1:OPTIONS.N_t) ,2) ...
            + sum(OPTIONS.C_ss(1:2,1).'* startup ,2) + sum(OPTIONS.C_ss(1:2,1)) + benders_cut )

    subject to
        % startup detect
        startup(1,1:OPTIONS.N_t-1) >= (master_delta(1,2:OPTIONS.N_t) - master_delta(1,1:OPTIONS.N_t-1))
        startup(2,1:OPTIONS.N_t-1) >= (master_delta(2,2:OPTIONS.N_t) - master_delta(2,1:OPTIONS.N_t-1))
        
        %  generator switch time
        for t_index = 1:1:OPTIONS.N_t - OPTIONS.Tmin_g+1
            master_delta(1,t_index:t_index+OPTIONS.Tmin_g-1) >= OPTIONS.Tmin_g * startup(t_index)
            master_delta(2,t_index:t_index+OPTIONS.Tmin_g-1) >= OPTIONS.Tmin_g * startup(t_index)
        end
                 
        %  Redundent_switch
        master_redundent_switch(1,:) + master_redundent_switch(2,:) == ones(1, OPTIONS.N_t)
        master_redundent_switch(3,:) + master_redundent_switch(4,:) == ones(1, OPTIONS.N_t)
        
%         % switch change detect
%         switch_change(1,1:OPTIONS.N_t-1) >= (master_redundent_switch(1,2:OPTIONS.N_t) - master_redundent_switch(1,1:OPTIONS.N_t-1))
%         switch_change(2,1:OPTIONS.N_t-1) >= (master_redundent_switch(2,2:OPTIONS.N_t) - master_redundent_switch(2,1:OPTIONS.N_t-1))
%         switch_change(3,1:OPTIONS.N_t-1) >= (master_redundent_switch(3,2:OPTIONS.N_t) - master_redundent_switch(3,1:OPTIONS.N_t-1))
%         switch_change(4,1:OPTIONS.N_t-1) >= (master_redundent_switch(4,2:OPTIONS.N_t) - master_redundent_switch(4,1:OPTIONS.N_t-1))
%         
%         %  Redundent switch time
%         for t_index = 1:1:OPTIONS.N_t - OPTIONS.Tmin_sw+1
%             master_delta(1,t_index:t_index+OPTIONS.Tmin_sw-1) >= OPTIONS.Tmin_sw * startup(t_index)
%             master_delta(2,t_index:t_index+OPTIONS.Tmin_sw-1) >= OPTIONS.Tmin_sw * startup(t_index)
%         end
%         
         % benders cuts
%         for index_benders = max([size(total_sub,2)-2 1] ):size(total_sub,2)
        for index_benders = 1:size(total_sub,2)
            benders_cut >=  total_sub(index_benders).sub_optval ...
                            + total_dual(index_benders).delta(1,:)*(master_delta(1,:).' - total_sub(index_benders).delta(1,:).')...
                            + total_dual(index_benders).delta(2,:)*(master_delta(2,:).' - total_sub(index_benders).delta(2,:).')... 
                            + total_dual(index_benders).switch(1,:)*(master_redundent_switch(1,:).' - total_sub(index_benders).redundent_sw(1,:).')...
                            + total_dual(index_benders).switch(2,:)*(master_redundent_switch(2,:).' - total_sub(index_benders).redundent_sw(2,:).')... 
                            + total_dual(index_benders).switch(3,:)*(master_redundent_switch(3,:).' - total_sub(index_benders).redundent_sw(3,:).')...
                            + total_dual(index_benders).switch(4,:)*(master_redundent_switch(4,:).' - total_sub(index_benders).redundent_sw(4,:).');
%                             + total_dual(index_benders).delta(1,:)*master_delta(1,:).' - total_dual(index_benders).delta(1,:)*total_sub(index_benders).delta(1,:).'...
%                             + total_dual(index_benders).delta(2,:)*master_delta(2,:).' - total_dual(index_benders).delta(2,:)*total_sub(index_benders).delta(2,:).' ...
%                             + total_dual(index_benders).switch(1,:)*master_redundent_switch(1,:).' - total_dual(index_benders).switch(1,:)*total_sub(index_benders).redundent_sw(1,:).'...
%                             + total_dual(index_benders).switch(2,:)*master_redundent_switch(2,:).' - total_dual(index_benders).switch(2,:)*total_sub(index_benders).redundent_sw(2,:).';
        end
 
        switch accelerate_flag
            case 1 
                % only power range
                % speedup constraints: power range
                master_delta(1,1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + master_delta(2,1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) +  OPTIONS.Pb_Max * ones(1,OPTIONS.N_t) >= OPTIONS.P_L_TIME_off + Ppr
                master_delta(1,1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + master_delta(2,1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) +  OPTIONS.Pb_Min * ones(1,OPTIONS.N_t) <= OPTIONS.P_L_TIME_off + Ppr
            case 2 
                % only lower bound
                % speedup constraint: lower bound
                benders_cut >= benders_cut_lowerbound
            case 3 
                % have two constraints
                % speedup constraints: power range
                master_delta(1,1:OPTIONS.N_t)*OPTIONS.Pg_Max(1) + master_delta(2,1:OPTIONS.N_t)*OPTIONS.Pg_Max(2) +  OPTIONS.Pb_Max * ones(1,OPTIONS.N_t) >= OPTIONS.P_L_TIME_off + Ppr
                master_delta(1,1:OPTIONS.N_t)*OPTIONS.Pg_Min(1) + master_delta(2,1:OPTIONS.N_t)*OPTIONS.Pg_Min(2) +  OPTIONS.Pb_Min * ones(1,OPTIONS.N_t) <= OPTIONS.P_L_TIME_off + Ppr

                % speedup constraint: lower bound
                benders_cut >= benders_cut_lowerbound
            otherwise % non-constriant
        end        
cvx_end

redundent_sw = master_redundent_switch;

% disp('optimal');
% disp(cvx_optval);

end
