% clear;
function [optimal_cost_comparison ] = cost_optimization_offline_benders_4G_2E( voya_distance )
global OPTIONS Parameter
 
if ~exist('voya_distance', 'var')
    OPTIONS.Distance = 130;
else
    OPTIONS.Distance = voya_distance;
end

OPTIONS.velocity = [25 0];
OPTIONS.N_e = 2;
OPTIONS.N_g = 4; 
OPTIONS.N_t = 8;

OPTIONS.velocity_avg = OPTIONS.Distance/OPTIONS.N_t;
OPTIONS.P_pr_avg = (OPTIONS.velocity_avg).^3*2.2e-3;
OPTIONS.Pg_Max(1) = 8;
OPTIONS.Pg_Min(1) = 1;
OPTIONS.Pg_Max(2) = 4;
OPTIONS.Pg_Min(2) = 0.5;

OPTIONS.Ppr_Max = (OPTIONS.velocity(1)).^3*2.2e-3;
OPTIONS.Pb_Max(1) = 1;
OPTIONS.Pb_Min(1) = -1;
OPTIONS.E_Max = [3 3];

OPTIONS.P_L = [2.7 0.9]; % P_Generater
OPTIONS.P_L_Scale = [0.5 0.6 0.8 0.82 0.7 0.6 0.4 0.35 0.3 0.33 0.4 0.5 0.4]; 

% the load demand without random feature
P_L_Scale_on = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline'); 
% the load demand with random feature
P_L_Scale_off = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline') + 0.1*rand(1, 25);

OPTIONS.P_L_TIME_off = sum(OPTIONS.P_L.'* P_L_Scale_off(:,1:OPTIONS.N_t), 1);
OPTIONS.P_L_TIME_on= sum(OPTIONS.P_L.'* P_L_Scale_on(:,1:OPTIONS.N_t), 1);

OPTIONS.Coupled_load = zeros(2, OPTIONS.N_t);
OPTIONS.Coupled_load(1, :) = 1 * OPTIONS.P_L_TIME_on(1,1:OPTIONS.N_t)./6;
OPTIONS.Coupled_load(2, :) = 1 * OPTIONS.P_L_TIME_on(1,1:OPTIONS.N_t)/6;
Redundent_switch(1,1:2) = 1;

% generator function parameters
Parameter.G(1,1:3) = [10 30 220];
Parameter.G(2,1:3) = [14.5 12 170];
Parameter.E(1,1:3) = [1 90 0];
Parameter.alpha = 1/6;
Parameter.C_ss = 90;
Parameter.R_G = 1;
Parameter.error = 1e-3;

load_information = [OPTIONS.P_L_TIME_off; OPTIONS.P_L_TIME_on; OPTIONS.Coupled_load];
save('load_information');

lambda_delta = zeros(4,OPTIONS.N_t);
lambda_Pg = zeros(4,OPTIONS.N_t);
delta_g = ones(4,OPTIONS.N_t);
% Pb = OPTIONS.Pb_Max/2*ones(1,OPTIONS.N_t);
delta = ones(OPTIONS.N_g, OPTIONS.N_t);

%%  operation_mode:
% 0 (normal wo PPA ESMC) 
% 1 (normal w PPA wo ESMC) 
% 2 (normal wo PPA w ESMC) 
% 3 (normal w PPA w ESMC)

% 4 (Fault wo PPA ESMC) 
% 5 (Fault w PPA wo ESMC) 
% 6 (Fault wo PPA w ESMC) 
% 7 (Fault w PPA w ESMC)

% for operation_mode = 5:5
%     [objval, Pg, Pb, Ppr ] = cost_optimization_problem( delta, Redundent_switch, operation_mode );
%     optimal_cost_comparison(1,operation_mode+1) = objval;
% end
% 
% optimal_cost_comparison(2,:) = (optimal_cost_comparison(1,:) -optimal_cost_comparison(1))*100/optimal_cost_comparison(1);
% b_benders = zeros(1,index_sm);
b_benders = 0;

for index_sm = 1:20
    for operation_mode = 7:7
        [objval_upperbound, Pg, Pb, Ppr, cvx_optval ] = cost_optimization_subproblem( delta, Redundent_switch, operation_mode );
        [objval_lowerbound, Pg, delta, Redundent_switch, b_benders ] = cost_optimization_masterproblem( Pg, Pb, Ppr, delta, Redundent_switch, operation_mode, cvx_optval, b_benders );
        objval_upperbound - objval_lowerbound
        if  objval_upperbound - objval_lowerbound < Parameter.error
%             break;
        end
    end
end
end


function [objval_upperbound, Pg, Pb, Ppr, cvx_optval ] = cost_optimization_subproblem( delta, Redundent_switch, Pb, operation_mode )  
global OPTIONS Parameter
%% subproblem
% cvx_begin
cvx_begin quiet 
    variable Ppr(1,OPTIONS.N_t) nonnegative
%     variable Pb(2,OPTIONS.N_t)
%     variable E(2,OPTIONS.N_t) nonnegative
%     variable delta(N_g,N_t) binary
    variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    minimize( sum(  sum(Parameter.G(1,1)* power(Pg(1:2,1:OPTIONS.N_t),2) ,1) ...
                    + sum(Parameter.G(2,1)* power(Pg(3:OPTIONS.N_g,1:OPTIONS.N_t),2) ,1) ) )
%                     + sum(Parameter.E(1,1)* power(Pb(1:2,1:OPTIONS.N_t),2), 1) ...
%                     + sum(Parameter.E(1,2)*Pb(1:2,1:OPTIONS.N_t),1) ...
%                     + sum(Parameter.E(1,3)*ones(2,OPTIONS.N_t),1) ) )
    subject to
        % the range constraints of all the variables
        Pg(1,1:OPTIONS.N_t) <= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1)
        Pg(2,1:OPTIONS.N_t) <= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1)
        Pg(3,1:OPTIONS.N_t) <= delta(3,1:OPTIONS.N_t) * OPTIONS.Pg_Max(2)
        Pg(4,1:OPTIONS.N_t) <= delta(4,1:OPTIONS.N_t) * OPTIONS.Pg_Max(2)
%         Pg(1,1:OPTIONS.N_t) >= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)
%         Pg(2,1:OPTIONS.N_t) >= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)

        Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) >= -Parameter.R_G
        Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) >= -Parameter.R_G
        Pg(3,2:OPTIONS.N_t) -Pg(3,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(3,2:OPTIONS.N_t) -Pg(3,1:OPTIONS.N_t-1) >= -Parameter.R_G
        Pg(4,2:OPTIONS.N_t) -Pg(4,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(4,2:OPTIONS.N_t) -Pg(4,1:OPTIONS.N_t-1) >= -Parameter.R_G

        Ppr(1,1:OPTIONS.N_t) <= OPTIONS.Ppr_Max * ones(1,OPTIONS.N_t)
%         Pb(1,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
%         Pb(1,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
%         Pb(2,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
%         Pb(2,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
      
        
%         E(1,1:OPTIONS.N_t) <=  OPTIONS.E_Max(1) * ones(1,OPTIONS.N_t)
%         E(1,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)
%         E(2,1:OPTIONS.N_t) <=  OPTIONS.E_Max(2) * ones(1,OPTIONS.N_t)
%         E(2,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)
        
%         % ESM output power and the capacity constraints
%         OPTIONS.E_Max(1) - Pb(1,1) == E(1,1)
%         OPTIONS.E_Max(2) - Pb(2,1) == E(2,1)
%         for t_index = 1:OPTIONS.N_t-1
%             E(1,t_index) - Pb(1,t_index+1) == E(1,t_index+1)
%             E(2,t_index) - Pb(2,t_index+1) == E(2,t_index+1)
%         end
        
        % system power balance
        if operation_mode <= 3 
            for t_index = 1:OPTIONS.N_t
                 OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == sum( Pg(1:4,t_index) ) + sum(Pb(1:2,t_index))
            end
        elseif operation_mode <= 7
            for t_index = 1:4
                 OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == sum( Pg(1:4,t_index) ) + sum(Pb(1:2,t_index))
            end
            for t_index = 5:OPTIONS.N_t
                Redundent_switch*OPTIONS.Coupled_load(:,t_index) +  (1-Parameter.alpha - 2/6) * OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == sum(Pg(1:2,t_index)) + Pb(1,t_index)
                ~Redundent_switch*OPTIONS.Coupled_load(:,t_index) + Parameter.alpha * OPTIONS.P_L_TIME_off(1,t_index)  == sum(Pg(3:4,t_index)) + Pb(2,t_index)
            end
        end

        % voyage planning            
        if operation_mode ==0                    
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
%             Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==1
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance                    
%             Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==2                    
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
        elseif operation_mode ==3
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance 
        elseif operation_mode ==4
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
%             Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==5
            sum((Ppr(5:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - 4*( OPTIONS.P_pr_avg ./2.2e-3).^(1/3)
            Ppr(1,1:4) == OPTIONS.P_pr_avg 
%             Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==6
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance 
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
%             Pb(1:2,1:4) == 0
        else
            sum((Ppr(5:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - 4*( OPTIONS.P_pr_avg ./2.2e-3).^(1/3)
%             Pb(1:2,1:4) == 0
            Ppr(1,1:4) == OPTIONS.P_pr_avg
        end
        
cvx_end


y = size(find(delta(:,2:OPTIONS.N_t) - delta(:,1:OPTIONS.N_t-1)==1),2);
objval_upperbound= cvx_optval + y*Parameter.C_ss(1)  ...
%                     + sum(sum(Parameter.G(1,2)*Pg(1:2,1:OPTIONS.N_t) 
                    + Parameter.G(1,3)*delta(1:2,1:OPTIONS.N_t)  ...
%                     + Parameter.G(2,2)*Pg(3:OPTIONS.N_g,1:OPTIONS.N_t) 
                    + Parameter.G(2,3)*delta(3:OPTIONS.N_g,1:OPTIONS.N_t),1),2);
% %% FIGURE PLOT
% figure
% hold on
% bar([ Pg(1,1:OPTIONS.N_t); Pg(2,1:OPTIONS.N_t); Pg(3,1:OPTIONS.N_t); Pg(4,1:OPTIONS.N_t); Pb(1,1:OPTIONS.N_t); Pb(2,1:OPTIONS.N_t)].','stacked');
% plot(Ppr(1:OPTIONS.N_t),'linewidth',1.5);
% plot(OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t),'linewidth',1.5);
% 
% xlim([0 OPTIONS.N_t+1]);
% legend('P_{G1}','P_{G2}','P_{pr}','P_{l}','Orientation','horizontal');
% ylabel('Active Power (MW)');
% xlabel('Time (hours)');
% 
% legend('P_{g_1}','P_{g_2}','P_{b_1}','P_{b_2}','P_{PR}','P_{L}');
% hold off

end

function [objval_lowerbound, Pg, delta, Redundent_switch, b_benders ] = cost_optimization_masterproblem( Pg, Pb, Ppr, delta, Redundent_switch, operation_mode, cvx_optval, b_benders )
global OPTIONS Parameter

%% dual variable and upperbound
lambda_B = -2*Parameter.E(1,1)*Pb - Parameter.E(1,2);
lambda_Pg(1,1:OPTIONS.N_t) = -2*Parameter.G(1,1)*Pg(1,1:OPTIONS.N_t) - lambda_B(1,:) ;
lambda_Pg(2,1:OPTIONS.N_t) = -2*Parameter.G(1,1)*Pg(2,1:OPTIONS.N_t) - lambda_B(1,:) ;
lambda_Pg(3,1:OPTIONS.N_t) = -2*Parameter.G(2,1)*Pg(3,1:OPTIONS.N_t) - lambda_B(2,:) ;
lambda_Pg(4,1:OPTIONS.N_t) = -2*Parameter.G(2,1)*Pg(4,1:OPTIONS.N_t) - lambda_B(2,:) ;
lambda_delta = zeros(OPTIONS.N_g,OPTIONS.N_t);
for index_t = 1:OPTIONS.N_t
    lambda_Px(1, index_t) = lambda_B(1, index_t).'*OPTIONS.Coupled_load(1, index_t);
    lambda_Sy(1, index_t) = lambda_B(2, index_t).'*OPTIONS.Coupled_load(2, index_t);
end
% 
% cvx_optval_2 =  sum(Parameter.G(1,1)* power(Pg(2,1:OPTIONS.N_t),2)  + Parameter.G(1,2)*Pg(2,1:OPTIONS.N_t) + Parameter.G(1,3)*ones(1,OPTIONS.N_t) ...
%                     + Parameter.E(1,1)* power(Pb(2,1:OPTIONS.N_t),2) );
% y_2 = size(find(delta(2,2:OPTIONS.N_t) - delta(1,1:OPTIONS.N_t-1)==1),2);
% objval_upperbound_2= cvx_optval_2 + y_2*Parameter.C_ss(1) + sum( Parameter.E(1,2)*Pb(2,1:OPTIONS.N_t) + Parameter.E(1,3)*ones(1,OPTIONS.N_t)) ;

%% build master problem
%% delta(4,4*Nt); y(4,4*(Nt-1)); Pg(4,4*Nt); mu(1); redundant_switch(1,4) 
% % total number = 3*N_g*Nt -N_g +5
% delta and Pg_min
A_delta_bound = zeros(2*OPTIONS.N_g*OPTIONS.N_t, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5);
for index_t= 0:OPTIONS.N_t-1
    for index_g = 0: OPTIONS.N_g/2-1
            A_delta_bound(8*index_t+2*index_g+1, index_t+index_g*OPTIONS.N_t+1) = -8;
            A_delta_bound(8*index_t+2*index_g+1, index_t+index_g*OPTIONS.N_t+8*OPTIONS.N_t-3) = 1;
            A_delta_bound(8*index_t+2*index_g+2, index_t+index_g*OPTIONS.N_t+1) = 0.5;
            A_delta_bound(8*index_t+2*index_g+2, index_t+index_g*OPTIONS.N_t+8*OPTIONS.N_t-3) = -1;
%             A_delta_bound(8*index_t+2*index_g+2, index_t+index_g*OPTIONS.N_t+1) = 0;
%             A_delta_bound(8*index_t+2*index_g+2, index_t+index_g*OPTIONS.N_t+8*OPTIONS.N_t-3) = 0;
    end    
    for index_g = OPTIONS.N_g/2: OPTIONS.N_g-1
            A_delta_bound(8*index_t+2*index_g+1, index_t+index_g*OPTIONS.N_t+1) = -4;
            A_delta_bound(8*index_t+2*index_g+1, index_t+index_g*OPTIONS.N_t+8*OPTIONS.N_t-3) = 1;
            A_delta_bound(8*index_t+2*index_g+2, index_t+index_g*OPTIONS.N_t+1) = 0.5;
            A_delta_bound(8*index_t+2*index_g+2, index_t+index_g*OPTIONS.N_t+8*OPTIONS.N_t-3) = -1;
%             A_delta_bound(8*index_t+2*index_g+2, index_t+index_g*OPTIONS.N_t+1) = 0;
%             A_delta_bound(8*index_t+2*index_g+2, index_t+index_g*OPTIONS.N_t+8*OPTIONS.N_t-3) = 0;
    end
end
% b_delta_bound = (Pg(2,1:OPTIONS.N_t).')./OPTIONS.Pg_Min(2);
b_delta_bound = zeros(2*OPTIONS.N_g*OPTIONS.N_t, 1);

% delta and y
A_delta = zeros((OPTIONS.N_t-1)*OPTIONS.N_g, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5);
for index_t = 0:OPTIONS.N_t-2
    for index_g = 0: OPTIONS.N_g-1
        A_delta(4*index_t+index_g+1, index_t+index_g*(OPTIONS.N_t)+1) = -1;
        A_delta(4*index_t+index_g+1, index_t+index_g*(OPTIONS.N_t)+2) = 1;
        A_delta(4*index_t+index_g+1, index_t+index_g*(OPTIONS.N_t-1)+4*OPTIONS.N_t+1) = -1;
    end
end
b_delta = zeros((OPTIONS.N_t-1)*OPTIONS.N_g,1);

% benders cuts
A_benders = zeros(1, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5);
for index_g = 0:OPTIONS.N_g-1
    A_benders(1, (index_g*OPTIONS.N_t+1):(index_g*OPTIONS.N_t+OPTIONS.N_t)) = lambda_delta(index_g+1,:); 
end
for index_g = 0:OPTIONS.N_g-1
    A_benders(1,index_g*OPTIONS.N_t+8*OPTIONS.N_t-3:index_g*OPTIONS.N_t+9*OPTIONS.N_t-4) = lambda_Pg(index_g+1,:);
end
A_benders(1, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 1) = -1; 

% do not have the value of last loop
% if index_sm == 1
benders_cut = 0;
for index_g = 1:OPTIONS.N_g
    benders_cut = benders_cut+ lambda_Pg(index_g,:)*Pg(index_g,:)'+ lambda_delta(1,:)*delta(1,:)';
end
tmp_b_benders =  - cvx_optval  + benders_cut;
if b_benders == 0
    b_benders = tmp_b_benders;
else
    b_benders = min([tmp_b_benders b_benders]);
end
% else
%     b_benders =  - cvx_optval + sum(lambda_delta*delta') + sum(lambda_Pg*Pg');
% end


A = [A_delta; A_benders; A_delta_bound];
b = [b_delta; b_benders; b_delta_bound];

% balance equality constraint
Aeq_B = zeros(OPTIONS.N_t, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5);
beq_B = zeros(OPTIONS.N_t, 1);
% for index_t = 0: 0
for index_t = 0: OPTIONS.N_t-1
    Aeq_B(2*index_t+1, index_t+0*OPTIONS.N_t+8*OPTIONS.N_t-4+1) = 1;
    Aeq_B(2*index_t+1, index_t+1*OPTIONS.N_t+8*OPTIONS.N_t-4+1) = 1;
    Aeq_B(2*index_t+1, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 2) = -OPTIONS.Coupled_load(1,index_t+1);
    Aeq_B(2*index_t+1, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 3) = -OPTIONS.Coupled_load(2,index_t+1);
    Aeq_B(2*index_t+1, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 4) = -OPTIONS.Coupled_load(1,index_t+1);
    Aeq_B(2*index_t+1, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5) = -OPTIONS.Coupled_load(2,index_t+1);
    Aeq_B(2*index_t+2, index_t+2*OPTIONS.N_t+8*OPTIONS.N_t-4+1) = 1;
    Aeq_B(2*index_t+2, index_t+3*OPTIONS.N_t+8*OPTIONS.N_t-4+1) = 1;
    Aeq_B(2*index_t+2, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 2) = -OPTIONS.Coupled_load(1,index_t+1);
    Aeq_B(2*index_t+2, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 3) = -OPTIONS.Coupled_load(2,index_t+1);
    Aeq_B(2*index_t+2, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 4) = -OPTIONS.Coupled_load(1,index_t+1);
    Aeq_B(2*index_t+2, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5) = -OPTIONS.Coupled_load(2,index_t+1);

    beq_B(2*index_t+1, 1) =  (1-Parameter.alpha - 2/6) * OPTIONS.P_L_TIME_off(1,index_t+1) + Ppr(index_t+1) -Pb(1,index_t+1);
    beq_B(2*index_t+2, 1) =  Parameter.alpha * OPTIONS.P_L_TIME_off(1,index_t+1) - Pb(2,index_t+1);
end

% redundancy switch constraint
Aeq_SW = zeros(OPTIONS.N_t, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5);
for index_t = 0: OPTIONS.N_t-1
    Aeq_SW(2*index_t+1, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 2) = 1;
    Aeq_SW(2*index_t+1, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 4) = 1;
    Aeq_SW(2*index_t+2, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 3) = 1;
    Aeq_SW(2*index_t+2, 3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5) = 1;
end
beq_SW = ones(2*OPTIONS.N_t,1);

Aeq = [ Aeq_B; Aeq_SW];
beq = [ beq_B; beq_SW];
% Aeq = [ Aeq_SW];
% beq = [ beq_SW];

% upper and lower bound
lb = [zeros(1,4*OPTIONS.N_t) zeros(1,4*(OPTIONS.N_t-1)) zeros(1,2*OPTIONS.N_t)                   zeros(1,2*OPTIONS.N_t)                    0    zeros(1,4)];
ub = [ones(1,4*OPTIONS.N_t)  ones(1,4*(OPTIONS.N_t-1))  OPTIONS.Pg_Max(1)*ones(1,2*OPTIONS.N_t)  OPTIONS.Pg_Max(2)*ones(1,2*OPTIONS.N_t)   inf  ones(1,4)];

f = [Parameter.G(1,3)*ones(1,2*OPTIONS.N_t) Parameter.G(2,3)*ones(1,2*OPTIONS.N_t) Parameter.C_ss*ones(1,4*(OPTIONS.N_t-1)) ...
    Parameter.G(1,2)*ones(1,2*OPTIONS.N_t) Parameter.G(2,2)*ones(1,2*OPTIONS.N_t) 1 zeros(1,4)];

intcon = [ 1: 8*(OPTIONS.N_t)-4  3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 2:3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 5];
% intcon = [];
options = optimoptions('intlinprog','Display','iter');
[x,objval_lowerbound] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, options);
% x = x.';
% objval_lowerbound = f(1,OPTIONS.N_t+1:2*OPTIONS.N_t-1)*x(OPTIONS.N_t+1:2*OPTIONS.N_t-1,1) + f(1,2*OPTIONS.N_t:3*OPTIONS.N_t-1)*x(2*OPTIONS.N_t:3*OPTIONS.N_t-1,1) + x(end);

delta =  round([x(1:1*OPTIONS.N_t,1).'; x(1*OPTIONS.N_t+1:2*OPTIONS.N_t,1).'; x(2*OPTIONS.N_t+1:3*OPTIONS.N_t,1).'; x(3*OPTIONS.N_t+1:4*OPTIONS.N_t,1).';]);
Pg(1,:) =  x(0*OPTIONS.N_t+8*OPTIONS.N_t-4+1:1*OPTIONS.N_t+8*OPTIONS.N_t-4,1).';
Pg(2,:) =  x(1*OPTIONS.N_t+8*OPTIONS.N_t-4+1:2*OPTIONS.N_t+8*OPTIONS.N_t-4,1).';
Pg(3,:) =  x(2*OPTIONS.N_t+8*OPTIONS.N_t-4+1:3*OPTIONS.N_t+8*OPTIONS.N_t-4,1).';
Pg(4,:) =  x(3*OPTIONS.N_t+8*OPTIONS.N_t-4+1:4*OPTIONS.N_t+8*OPTIONS.N_t-4,1).';
Redundent_switch = x(3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 2:3*OPTIONS.N_g*OPTIONS.N_t - OPTIONS.N_g + 3,1).';

% %% FIGURE PLOT
% figure
% hold on
% bar([ Pg(1,1:OPTIONS.N_t);  Pg(2,1:OPTIONS.N_t);  Pg(3,1:OPTIONS.N_t);  Pg(4,1:OPTIONS.N_t); Pb(1,1:OPTIONS.N_t); Pb(2,1:OPTIONS.N_t)].','stacked');
% plot(Ppr(1:OPTIONS.N_t),'linewidth',1.5);
% plot(OPTIONS.P_L_TIME_on(1,1:OPTIONS.N_t),'linewidth',1.5);
% 
% xlim([0 OPTIONS.N_t+1]);
% % ylim([0 10]);
% % plot([1 12], [P_prop P_prop], 'linewidth', 2, 'Color',[0.0,0.6,0.9]);
% % plot(1:12, P_load, 'linewidth', 2, 'Color',[0.67,0,1.0]);
% legend('P_{G1}','P_{G2}','P_{G3}','P_{G4}','P_{E1}','P_{E2}','P_{pr}','P_{l}','Orientation','horizontal');
% ylabel('Active Power (MW)');
% xlabel('Time (hours)');
% 
% legend('P_{g_1}','P_{b_1}','P_{PR}','P_{L}');
% hold off

end



function [cvx_optval, Pg, Pb, Ppr ] = cost_optimization( delta, Redundent_switch, operation_mode )  
global OPTIONS Parameter
%% subproblem
cvx_begin
% cvx_begin quiet 
    variable Ppr(1,OPTIONS.N_t) nonnegative
    variable Pb(2,OPTIONS.N_t)
    variable E(2,OPTIONS.N_t) nonnegative
%     variable delta(N_g,N_t) binary
    variable Pg(OPTIONS.N_g, OPTIONS.N_t) nonnegative
    minimize( sum(  sum(Parameter.G(1,1)* power(Pg(1:2,1:OPTIONS.N_t),2)  + Parameter.G(1,2)*Pg(1:2,1:OPTIONS.N_t)  ...
                        + Parameter.G(1,3)*ones(2,OPTIONS.N_t) ,1) ...
                  + sum(Parameter.G(2,1)* power(Pg(3:OPTIONS.N_g,1:OPTIONS.N_t),2)  + Parameter.G(2,2)*Pg(3:OPTIONS.N_g,1:OPTIONS.N_t)  ...
                        + Parameter.G(2,3)*ones(2,OPTIONS.N_t) ,1) ...
                        + sum(Parameter.E(1,1)* power(Pb(1:2,1:OPTIONS.N_t),2), 1), 2 ) )
    subject to
        % the range constraints of all the variables
        Pg(1,1:OPTIONS.N_t) <= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1)
        Pg(2,1:OPTIONS.N_t) <= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1)
        Pg(3,1:OPTIONS.N_t) <= delta(3,1:OPTIONS.N_t) * OPTIONS.Pg_Max(2)
        Pg(4,1:OPTIONS.N_t) <= delta(4,1:OPTIONS.N_t) * OPTIONS.Pg_Max(2)
%         Pg(1,1:OPTIONS.N_t) >= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)
%         Pg(2,1:OPTIONS.N_t) >= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)

        Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) >= -Parameter.R_G
        Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) >= -Parameter.R_G
        Pg(3,2:OPTIONS.N_t) -Pg(3,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(3,2:OPTIONS.N_t) -Pg(3,1:OPTIONS.N_t-1) >= -Parameter.R_G
        Pg(4,2:OPTIONS.N_t) -Pg(4,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(4,2:OPTIONS.N_t) -Pg(4,1:OPTIONS.N_t-1) >= -Parameter.R_G

        Ppr(1,1:OPTIONS.N_t) <= OPTIONS.Ppr_Max * ones(1,OPTIONS.N_t)
        Pb(1,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
        Pb(1,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
        Pb(2,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
        Pb(2,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
      
        
        E(1,1:OPTIONS.N_t) <=  OPTIONS.E_Max(1) * ones(1,OPTIONS.N_t)
        E(1,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)
        E(2,1:OPTIONS.N_t) <=  OPTIONS.E_Max(2) * ones(1,OPTIONS.N_t)
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
                 OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == sum( Pg(1:4,t_index) ) + sum(Pb(1:2,t_index))
            end
        elseif operation_mode <= 7
            for t_index = 1:4
                 OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == sum( Pg(1:4,t_index) ) + sum(Pb(1:2,t_index))
            end
            for t_index = 5:OPTIONS.N_t
                Redundent_switch*OPTIONS.Coupled_load(:,t_index) +  Parameter.alpha * OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == sum(Pg(1:2,t_index)) + Pb(1,t_index)
                ~Redundent_switch*OPTIONS.Coupled_load(:,t_index) + (1-Parameter.alpha - 2/6) * OPTIONS.P_L_TIME_off(1,t_index)  == sum(Pg(3:4,t_index)) + Pb(2,t_index)
            end
        end

        % voyage planning            
        if operation_mode ==0                    
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
            Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==1
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance                    
            Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==2                    
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
        elseif operation_mode ==3
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance 
        elseif operation_mode ==4
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
            Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==5
            sum((Ppr(5:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - 4*( OPTIONS.P_pr_avg ./2.2e-3).^(1/3)
            Ppr(1,1:4) == OPTIONS.P_pr_avg 
            Pb(1:2,1:OPTIONS.N_t) == 0
        elseif operation_mode ==6
            sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance 
            Ppr(1,1:OPTIONS.N_t) == OPTIONS.P_pr_avg;
            Pb(1:2,1:4) == 0
        else
            sum((Ppr(5:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance - 4*( OPTIONS.P_pr_avg ./2.2e-3).^(1/3)
            Pb(1:2,1:4) == 0
            Ppr(1,1:4) == OPTIONS.P_pr_avg
        end
        
cvx_end


y = size(find(delta(1:2,2:OPTIONS.N_t) - delta(1:2,1:OPTIONS.N_t-1)==1),2);
cvx_optval = cvx_optval + y*Parameter.C_ss(1) + sum(sum( Parameter.E(1,2)*Pb(1:2,1:OPTIONS.N_t) + Parameter.E(1,3)*ones(2,OPTIONS.N_t),1),2) ;
  + Parameter.G(1,2)*Pg(1:2,1:OPTIONS.N_t) + Parameter.G(2,2)*Pg(3:OPTIONS.N_g,1:OPTIONS.N_t)

% %% FIGURE PLOT
% figure
% hold on
% bar([ Pg(1,1:OPTIONS.N_t); Pg(2,1:OPTIONS.N_t); Pg(3,1:OPTIONS.N_t); Pg(4,1:OPTIONS.N_t); Pb(1,1:OPTIONS.N_t); Pb(2,1:OPTIONS.N_t)].','stacked');
% plot(Ppr(1:OPTIONS.N_t),'linewidth',1.5);
% plot(OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t),'linewidth',1.5);
% 
% xlim([0 OPTIONS.N_t+1]);
% legend('P_{G1}','P_{G2}','P_{pr}','P_{l}','Orientation','horizontal');
% ylabel('Active Power (MW)');
% xlabel('Time (hours)');
% 
% legend('P_{g_1}','P_{g_2}','P_{b_1}','P_{b_2}','P_{PR}','P_{L}');
% hold off

end