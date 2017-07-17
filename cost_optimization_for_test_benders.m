function [optimal_cost_comparison ] = cost_optimization_benders_2G_2E_cvx( voya_distance )
global OPTIONS
 
if ~exist('voya_distance', 'var')
    OPTIONS.Distance = 70;
else
    OPTIONS.Distance = voya_distance;
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
operation_mode = 3;

initial_parameters();

lambda_delta = zeros(OPTIONS.N_g,OPTIONS.N_t);
lambda_Pg = zeros(OPTIONS.N_g,OPTIONS.N_t);
% startup variable used in master problem
delta_g = ones(OPTIONS.N_g,OPTIONS.N_t);
Pg = zeros(2,OPTIONS.N_t);
Pb = zeros(2,OPTIONS.N_t);
Redundent_switch(1,1:2) = [1 0];
objval_upperbound = 0;
objval_lowerbound = 0;

% startup variable used in subproblem
delta = ones(2, OPTIONS.N_t);
[result(1), Pg, Pb, Ppr ] = total_cost_optimization( operation_mode, delta, Redundent_switch, Pg, Pb );
[upperbound(1), sub_optval, Pg, Pb, Ppr ] = optimization_subproblem( operation_mode, delta, Redundent_switch, Pg, Pb );
[objval_lowerbound(1), delta, Redundent_switch, Pg, b_benders ] = optimization_masterproblem( operation_mode, delta, Redundent_switch, Pg, Pb, Ppr, sub_optval);


for startup_index = 1: 12+30
    delta = starup_slot(startup_index);
    [result(startup_index+1), Pg, Pb, Ppr ] = total_cost_optimization( operation_mode, delta, Redundent_switch, Pg, Pb );
    [upperbound(startup_index+1), cvx_optval, Pg, Pb, Ppr ] = optimization_subproblem( operation_mode, delta, Redundent_switch, Pg, Pb );
    [objval_lowerbound(startup_index+1), delta, Redundent_switch, Pg, b_benders ] = optimization_masterproblem( operation_mode, delta, Redundent_switch, Pg, Pb, Ppr, Pd, cvx_optval, b_benders );
end

save('result.mat','result');
save('upperbound.mat','upperbound');


% b_benders = 0;
% for index_sm = 1:20
%     for operation_mode = 7:7
%         [objval_upperbound(index_sm), cvx_optval, Pg, Pb, Ppr, Pd] = optimization_subproblem( operation_mode, delta, Redundent_switch, Pg, Pb );
%         if isnan(Pg(1,1))
%             break;
%         end
%         [objval_lowerbound(index_sm), delta, Redundent_switch, Pg, b_benders ] = optimization_masterproblem( operation_mode, delta, Redundent_switch, Pg, Pb, Ppr, Pd, cvx_optval, b_benders );
%         if isnan(Pg(1,1))
%             break;
%         end
%         ul_error = objval_upperbound(index_sm) - objval_lowerbound(index_sm);
%         if  objval_upperbound(index_sm) - objval_lowerbound(index_sm) < OPTIONS.error
% %             break;
%         end
%     end
% end


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
OPTIONS.N_t = 6;
% minmum time between startup and shuntdown
OPTIONS.Tmin = 1;
% time of faults happen
OPTIONS.Fault_time = 4;

OPTIONS.velocity_avg = OPTIONS.Distance/OPTIONS.N_t;
OPTIONS.P_pr_avg = (OPTIONS.velocity_avg).^3*2.2e-3;
OPTIONS.Pg_Max(1) = 8;
OPTIONS.Pg_Min(1) = 1;
OPTIONS.Pg_Max(2) = 4;
OPTIONS.Pg_Min(2) = 0.5;

OPTIONS.Xi = 1;
OPTIONS.Ppr_Max = (OPTIONS.velocity(1)).^3*2.2e-3;
OPTIONS.Pb_Max(1) = 1;
OPTIONS.Pb_Min(1) = -1;
OPTIONS.E_Max = [3 3];

OPTIONS.P_L = [2.7 0.9]; % P_Generater
OPTIONS.P_L_Scale = [0.5 0.6 0.8 0.82 0.7 0.6 0.4 0.35 0.3 0.33 0.4 0.5 0.4]; 

% the load demand without random feature
P_L_Scale_on = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline') + 0.1*rand(1, 25); 
% the load demand with random feature
% P_L_Scale_off = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline');
P_L_Scale_off = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline');

OPTIONS.P_L_TIME_off = sum(OPTIONS.P_L.'* P_L_Scale_off(:,1:OPTIONS.N_t), 1);
OPTIONS.P_L_TIME_on= sum(OPTIONS.P_L.'* P_L_Scale_on(:,1:OPTIONS.N_t), 1);

OPTIONS.Coupled_load = zeros(2, OPTIONS.N_t);
OPTIONS.Coupled_load(1, :) = 1 * OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t)./6;
OPTIONS.Coupled_load(2, :) = 1 * OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t)/6;

% generator function parameters
OPTIONS.G(1,1:3) = [10 30 220];
OPTIONS.G(2,1:3) = [14.5 12 170];
OPTIONS.E(1:2,1:3) = [1 90 0; 1 90 0];
OPTIONS.alpha = 1/6;
OPTIONS.C_ss = [90; 45];
OPTIONS.R_G = 1;
OPTIONS.error = 1e-3;

% load_information = [OPTIONS.P_L_TIME_off; OPTIONS.P_L_TIME_on; OPTIONS.Coupled_load];
% save('load_information');

end

%% heuristic algorithm
function [delta] = starup_slot(startup_index)
global OPTIONS
delta = ones(2, OPTIONS.N_t);
% one slot shutdown for two generators
if startup_index <= 12
    delta(startup_index) = 0;
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
function [ cvx_optval, Pg, Pb, Ppr ] = total_cost_optimization( operation_mode, delta, Redundent_switch,  Pg_m, Pb )  
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
                Redundent_switch*OPTIONS.Coupled_load(:,t_index) +  (1-OPTIONS.alpha - 2/6) * OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == (Pg(1,t_index)) + Pb(1,t_index) 
                ~Redundent_switch*OPTIONS.Coupled_load(:,t_index) + OPTIONS.alpha * OPTIONS.P_L_TIME_off(1,t_index)  == (Pg(2,t_index)) + Pb(2,t_index)
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
function [objval_upperbound, cvx_optval, Pg, Pb, Ppr ] = optimization_subproblem( operation_mode, delta_s, Redundent_switch,  Pg_m, Pb )  
global OPTIONS

lowerbound(1,:) = delta_s(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1);
lowerbound(2,:) = delta_s(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(2);

cvx_begin
% cvx_solver SeDuMi
cvx_begin quiet
    variable Ppr(1,OPTIONS.N_t) nonnegative
    variable Pb(2,OPTIONS.N_t)
    variable E(2,OPTIONS.N_t) nonnegative
    variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    minimize( sum(OPTIONS.G(1,1) * delta_s(1,1:OPTIONS.N_t) * power(Pg(1,1:OPTIONS.N_t).',2) ,2) ...
            + sum(OPTIONS.G(2,1) * delta_s(2,1:OPTIONS.N_t) * power(Pg(2,1:OPTIONS.N_t).',2) ,2) ...
            + sum(OPTIONS.G(1,2) * delta_s(1,1:OPTIONS.N_t) * power(Pg(1,1:OPTIONS.N_t).',1) ,2) ...
            + sum(OPTIONS.G(2,2) * delta_s(2,1:OPTIONS.N_t) * power(Pg(2,1:OPTIONS.N_t).',1) ,2) ...
            + OPTIONS.Xi * sum(OPTIONS.E(1:2,1).'* power(Pb(1:OPTIONS.N_e,1:OPTIONS.N_t),2) ,2) ...
            + OPTIONS.Xi * sum(OPTIONS.E(1:2,2).'* ones(OPTIONS.N_g,OPTIONS.N_t) ,2)  )
 
    subject to
        % the range constraints of all the variables
        Pg(1,1:OPTIONS.N_t) <= delta_s(1,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1) 
        Pg(2,1:OPTIONS.N_t) <= delta_s(2,1:OPTIONS.N_t) * OPTIONS.Pg_Max(2)
        Pg(1,1:OPTIONS.N_t) >= delta_s(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)
        Pg(2,1:OPTIONS.N_t) >= delta_s(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(2)

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
                Redundent_switch*OPTIONS.Coupled_load(:,t_index) +  (1-OPTIONS.alpha - 2/6) * OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == (Pg(1,t_index)) + Pb(1,t_index) 
                ~Redundent_switch*OPTIONS.Coupled_load(:,t_index) + OPTIONS.alpha * OPTIONS.P_L_TIME_off(1,t_index)  == (Pg(2,t_index)) + Pb(2,t_index)
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


startup = (delta_s(1:2,2:OPTIONS.N_t) - delta_s(1:2,1:OPTIONS.N_t-1)>=1);
startup_cost = sum(OPTIONS.C_ss(1:2,1).'* startup ,2);
one_item_cost = sum(OPTIONS.G(1,3) * delta_s(1,1:OPTIONS.N_t) ,2) + sum(OPTIONS.G(2,3) * delta_s(2,1:OPTIONS.N_t) ,2);

objval_upperbound = cvx_optval + startup_cost + one_item_cost;

disp('upperbound');
disp(objval_upperbound);

end

%% the master problem which is used to determine the redundent switches and ??
function [objval_lowerbound, delta_g, Redundent_switch, Pg_m, benders_cut ] = optimization_masterproblem( operation_mode, delta, Redundent_switch, Pg, Pb, Ppr, sub_optval, benders_cut, Pd )
global OPTIONS
 
if ~exist('benders_cut', 'var')
    benders_cut = 0;
end

if ~exist('Pd', 'var')
    Pd = 0;
end

%% dual variable and upperbound
lambda_B(1,:) = -2*OPTIONS.G(1,1)*Pg(1,:) - OPTIONS.G(1,2)*Pg(1,:);
lambda_B(2,:) = -2*OPTIONS.G(2,1)*Pg(2,:) - OPTIONS.G(2,2)*Pg(2,:);
% lambda_Pg(1,1:OPTIONS.N_t) = - 2*OPTIONS.G(1,1)*Pg(1,:) - lambda_B(1,:);
% lambda_Pg(2,1:OPTIONS.N_t) = - 2*OPTIONS.G(1,1)*Pg(2,:) - lambda_B(2,:);
lambda_delta(1,:) = -OPTIONS.G(1,1)*power(Pg(1,:),2) - OPTIONS.G(1,2)*Pg(1,:);
lambda_delta(2,:) = -OPTIONS.G(2,1)*power(Pg(2,:),2) - OPTIONS.G(2,2)*Pg(2,:);
for index_t = 1:OPTIONS.N_t
    lambda_Px(1, index_t) = lambda_B(1, index_t).'*OPTIONS.Coupled_load(1, index_t);
    lambda_Px(2, index_t) = lambda_B(1, index_t).'*OPTIONS.Coupled_load(2, index_t);
    lambda_Sy(1, index_t) = lambda_B(2, index_t).'*OPTIONS.Coupled_load(1, index_t);
    lambda_Sy(2, index_t) = lambda_B(2, index_t).'*OPTIONS.Coupled_load(2, index_t);
end

cvx_begin
cvx_solver Mosek
% cvx_begin quiet
    variable delta_master(OPTIONS.N_g, OPTIONS.N_t) binary
    variable startup(OPTIONS.N_g, OPTIONS.N_t-1) binary
    variable benders_cut 
%     variable E(2,OPTIONS.N_t) nonnegative
%     variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    minimize( sum(OPTIONS.G(1,3) * delta_master(1,1:OPTIONS.N_t) ,2) ...
            + sum(OPTIONS.G(2,3) * delta_master(2,1:OPTIONS.N_t) ,2) ...
            + sum(OPTIONS.C_ss(1:2,1).'* startup ,2) + benders_cut )
 
    subject to
        % startup detect
        startup(1,1:OPTIONS.N_t-1) >= (delta_master(1,2:OPTIONS.N_t) - delta_master(1,1:OPTIONS.N_t-1))
        startup(2,1:OPTIONS.N_t-1) >= (delta_master(2,2:OPTIONS.N_t) - delta_master(2,1:OPTIONS.N_t-1))

%         % the range constraints of all the variables
%         Pg(1,1:OPTIONS.N_t) <= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1)
%         Pg(2,1:OPTIONS.N_t) <= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Max(2)
%         Pg(1,1:OPTIONS.N_t) >= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)
%         Pg(2,1:OPTIONS.N_t) >= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(2)

%         %  Redundent_switch
%         Redundent_switch(1,1) + Redundent_switch(1,3) == 1
%         Redundent_switch(1,2) + Redundent_switch(1,4) == 1
        
         % benders cuts
        benders_cut >= sub_optval + lambda_delta(1,:)*delta_master(1,:).' - lambda_delta(1,:)*delta(1,:).' + lambda_delta(2,:)*delta_master(2,:).' - lambda_delta(2,:)*delta(2,:).'
%                        + sum(Redundent_switch_m(1,1)*lambda_Px(1,:) + Redundent_switch_m(1,2)*lambda_Px(1,:) - Redundent_switch(1,1)*lambda_Px(1,:) - Redundent_switch(1,2)*lambda_Px(1,:)...
%                        + Redundent_switch_m(1,3)*lambda_Sy(1,:) + Redundent_switch_m(1,4)*lambda_Sy(1,:) - ~Redundent_switch(1,1)*lambda_Sy(1,:) - ~Redundent_switch(1,2)*lambda_Sy(1,:));
                   
cvx_end

disp('optimal');
disp(cvx_optval);

end
