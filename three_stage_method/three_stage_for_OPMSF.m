function three_stage_for_OPMSF()
global OPTIONS
%% the adjusted parameters 
dbstop if error
if ~exist('time_slot', 'var')
    OPTIONS.N_t = 10;
else
    OPTIONS.N_t = time_slot;
end

if ~exist('voya_distance', 'var')
    OPTIONS.Distance = 150;
else
    OPTIONS.Distance = voya_distance;
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
    No_test = 3;
else
    No_test = No_test_in;
end

% initial the parameters
initial_parameters(parameter_test);
% initial the load demand
load_demand(operation_mode);

%% First stage: ESM power calculation at determined load demand
load_demand_first_stage = OPTIONS.P_L_TIME_off(1, 1:OPTIONS.N_t) + OPTIONS.Ppr_avg;
[Pb_first_stage] = dynamic_programming_first_stage(load_demand_first_stage);

%% Second stage: propulsion power deviation calculation
load_demand_second_stage = OPTIONS.P_L_TIME(1, 1:OPTIONS.N_t) - Pb_first_stage;
[Delta_Ppr_second_stage] = dynamic_programming_second_stage(load_demand_second_stage);

%% Third stage: meet the travel distance constraint by adjust propulsion power deviation 
load_demand_third_stage = load_demand_second_stage + Delta_Ppr_second_stage;
[Pg, Delta_Ppr_third_stage] = PSO(load_demand_third_stage);

total_P.Pg = Pg;
total_P.Pb = Pb_first_stage;
total_P.Ppr = OPTIONS.Ppr_avg + Delta_Ppr_second_stage + Delta_Ppr_third_stage;

cost_for_comparison = save_optimal_cost_information(No_test, near_opt_optimal, optimal_cost, total_P, delta_g, operation_mode );
end

%% Save optimal operation status
function [cost_for_comparison] = save_optimal_cost_information(No_test, near_opt_optimal, optimal_cost, total_P, delta_g, operation_mode )
global OPTIONS accelerate_flag
operating_cost = OPTIONS.G(1:OPTIONS.N_g, 1).' * power(total_P.Pg(1:OPTIONS.N_g, 1:OPTIONS.N_t), 2) ...
               + OPTIONS.G(1:OPTIONS.N_g, 2).' * total_P.Pg(1:OPTIONS.N_g, 1:OPTIONS.N_t) ...
               + OPTIONS.G(1:OPTIONS.N_g, 3).' * delta_g ...
               + OPTIONS.Xi_E * OPTIONS.E(1:OPTIONS.N_e, 1).' * power(total_P.Pb(1:OPTIONS.N_e, 1:OPTIONS.N_t), 2) ...
               + OPTIONS.Xi_E * OPTIONS.E(1:OPTIONS.N_e, 2).' * ones(2, OPTIONS.N_t);
% Penalty_cost_LS =  OPTIONS.Penalty_L * sum(total_P.load_shedding(1:2, 1:OPTIONS.N_t),1);
% Penalty_cost_D = 1.02 * OPTIONS.Penalty_D * reduced_distance;
startup_g =  delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t) - [zeros(OPTIONS.N_g,1) delta_g(1:OPTIONS.N_g, 1:OPTIONS.N_t-1)];
startup_cost = OPTIONS.C_ss(1:OPTIONS.N_g,1).'* (round(startup_g) > 0);
operating_cost = operating_cost + startup_cost;

% total cost obtained from the master problem
data.cost(1, 1) = optimal_cost(1, end);
% operating cost calculated based on the results such as Pg, Ppr and Pe
data.cost(2, 1:OPTIONS.N_t) = operating_cost;
% total operating cost of generators and ESMs
data.cost(3, 1) = sum(operating_cost, 2);
% penalty cost of load shedding
% data.cost(3, 2) = sum(Penalty_cost_LS, 2);
% penalty cost of reduced distance
% data.cost(3, 3) = Penalty_cost_D;
% total cost of operating cost, two terms of penalty cost
% data.cost(3, 4) = sum(operating_cost) + sum(Penalty_cost_LS) + Penalty_cost_D;

% data for parameter passing
cost_for_comparison = data.cost(3, 1);

% total_switching = sum( abs(redundent_sw(1,1) - OPTIONS.initial_redundent_sw(1,1)) ) + sum(abs(redundent_sw(3,1) - OPTIONS.initial_redundent_sw(3,1))) ...
%         + sum( abs(redundent_sw(1,2:OPTIONS.N_t) - redundent_sw(1,1:OPTIONS.N_t-1)) ) + sum(abs(redundent_sw(3,2:OPTIONS.N_t) - redundent_sw(3,1:OPTIONS.N_t-1)));

data.power(1:OPTIONS.N_g, 1:OPTIONS.N_t) =  total_P.Pg;
data.power(OPTIONS.N_g+1 : OPTIONS.N_g+OPTIONS.N_e, 1:OPTIONS.N_t) =  total_P.Pb;
data.power(5, 1:OPTIONS.N_t) =  total_P.Ppr;
% data.power(6:7, 1:OPTIONS.N_t) =  total_P.load_shedding;
data.power(7, 1:OPTIONS.N_t) =  sum(total_P.Pg + total_P.Pb, 1) - total_P.Ppr;
% data.power(8:8+OPTIONS.N_g-1, 1:OPTIONS.N_t) =  total_P.Pc;
% data.power(10,:) = OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t);

data.status(1:OPTIONS.N_g, 1:OPTIONS.N_t) = delta_g;
% data.status(OPTIONS.N_g+1:OPTIONS.N_g+2, 1:OPTIONS.N_t) = [redundent_sw(1,1:OPTIONS.N_t); redundent_sw(3,1:OPTIONS.N_t)];
% data.status(OPTIONS.N_g+3, 1) = total_switching;

% data.distance(1, 1:2) = [OPTIONS.Distance reduced_distance];

filename = ['3_Stage_mode_',num2str(operation_mode),'_D.',num2str(OPTIONS.Distance),...
            '_T.',num2str(OPTIONS.N_t),'_No.',num2str(No_test),'.mat'];
save(filename,'data');

end

%% Load data generation
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
    
    P_vs = P_vs_base * OPTIONS.P_vs_t_Scale;
    P_no = P_no_base * OPTIONS.P_no_t_Scale;
    P_total_vs = sum(P_vs);
    
    P_total_time = P_total_vs + P_no;
    OPTIONS.P_L_TIME_off = P_total_time(1:OPTIONS.N_t);
end

end

%% Inital the parameters for optimization
function [] = initial_parameters(parameter_test)
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
OPTIONS.Ppr_avg = (OPTIONS.velocity_avg).^3*2.2e-3;
OPTIONS.Pg_Max(1) = 8;
OPTIONS.Pg_Min(1) = 1;
OPTIONS.Pg_Max(2) = 4;
OPTIONS.Pg_Min(2) = 0.5;

OPTIONS.Xi_E = 1;
OPTIONS.Ppr_Max = (OPTIONS.velocity(1)).^3*2.2e-3;
OPTIONS.Pb_Max(1) = 1;
OPTIONS.Pb_Min(1) = -1;
OPTIONS.E_Max = [2 2];
OPTIONS.E_Min = [0.2 0.2];

% generator function parameters
switch parameter_test
    case 1 % 
        OPTIONS.G(1,1:3) = [10 19 170];
        OPTIONS.G(2,1:3) = [10 19 170];
        OPTIONS.G(3,1:3) = [28.5 28 190];
        
    case 2 % similar minimum SFC (160), large SFC deviation (250 and 150) between lower bound Pmin and efficient point Peff 
        OPTIONS.G(1,1:3) = [12 60 350];
        OPTIONS.G(2,1:3) = [12 60 350];
        OPTIONS.G(3,1:3) = [30 12 170];
        
    % Kanellos: OPMS in IEEE Transactions on Sustainable Energy, Electric power systems research, and Inventions     
    case 3 
        OPTIONS.G(1,1:3) = [5.4 61.5 390];
        OPTIONS.G(2,1:3) = [5.4 63 400];
        OPTIONS.G(3,1:3) = [13.1 12 430];
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
end