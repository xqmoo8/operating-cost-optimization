function [Pg, Pb, delta, Ppr, lambda_d, lambda_b1, lambda_b2] = cost_optimization_off()
global OPTIONS Parameter
% m = 20; n = 10; p = 4;
% A = randn(m,n); b = randn(m,1);
% C = randn(p,n); d = randn(p,1); e = rand;
% cvx_begin
%     variable x(n)
%     minimize( norm( A * x - b, 2 ) )
%     subject to
%         C * x == d
%         norm( x, Inf ) <= e
% cvx_end

OPTIONS.Distance = 100;
OPTIONS.velocity = [17 0];
OPTIONS.Operation_Time = 12;

OPTIONS.velocity_avg = OPTIONS.Distance/OPTIONS.Operation_Time;
OPTIONS.P_pr_avg = (OPTIONS.velocity_avg).^3*2.2e-3;
OPTIONS.Delta_P_pr = 2;

OPTIONS.P_L = [2.7 0.9]; % P_Generater
OPTIONS.P_L_Scale = [0.5 0.6 0.8 0.8 0.7 0.6 0.4 0.35 0.3 0.33 0.4 0.5]; 
OPTIONS.P_L_TIME = sum(OPTIONS.P_L.'* OPTIONS.P_L_Scale, 1);

OPTIONS.P_prop = 2.2e-3*(OPTIONS.Distance/12).^3;
OPTIONS.P_demand = OPTIONS.P_prop + OPTIONS.P_L_TIME;

OPTIONS.N_e = 2;
OPTIONS.N_g = 2; 
OPTIONS.N_t = OPTIONS.Operation_Time;


% generator function parameters
Parameter.G(1,1:3) = [13.5 10 300];
Parameter.G(2,1:3) = [6 30 250];
Parameter.E(1,1:3) = [150 2 0];
Parameter.alpha = 0.6;
Parameter.C_ss = 100;

N_e = OPTIONS.N_e;
N_g = OPTIONS.N_g;
N_t = OPTIONS.N_t;
C_ss = Parameter.C_ss;

Pg_Max(1:N_t) = 8;
Ppr_Max(1:N_t) = 12;
Pb_Max(1:N_t) = 1;
Pb_Min(1:N_t) = -1;
E_Max(1:N_t) = 2;
Pg_constant(1:N_t) = 1;

Pb(1,1:N_t) = 0;
E(2,1:N_t) = 2;

delta(1:2,1:N_t) = 1;
switch_PS(1:2) = 1;
error_primal_dual = 10;

% if error_primal_dual <= 1e-3
    
%% SUBPROBLEM OPTIMIZATION
cvx_begin
    variable Ppr(1,N_t) nonnegative
    variable Pb(2,N_t)
    variable E(2,N_t) nonnegative
%     variable delta(N_g,N_t) binary
    variable Pg(2,N_t) nonnegative
    minimize( sum(  Parameter.G(1,1)* power(Pg(1,1:N_t),2)  + Parameter.G(1,2)*Pg(1,1:N_t) + Parameter.G(1,3)*Pg_constant(1:N_t) ...
                    + Parameter.G(2,1)* power(Pg(2,1:N_t),2)  + Parameter.G(2,2)*Pg(2,1:N_t) + Parameter.G(2,3)*Pg_constant(1:N_t) ...
                    + Parameter.E(1,1)* power(Pb(1,1:N_t),2) + Parameter.E(1,1)* power(Pb(2,1:N_t),2) ) )
%                     + Parameter.E(1,1)* power(Pb(1,1:N_t),1) + Parameter.E(1,2)* power(Pb(1,1:N_t),2) + Parameter.E(1,1)* power(Pb(2,1:N_t),1) + Parameter.E(1,2)* power(Pb(2,1:N_t),2) ) )
    subject to         
        % the range constraints of all the variables
        Pg(1,1:N_t) <= delta(1,N_t).*Pg_Max(1:N_t)
        Pg(2,1:N_t) <= delta(2,N_t).*Pg_Max(1:N_t)
        Ppr(1:N_t) <= Ppr_Max(1:N_t)
        Pb(1,1:N_t) <= Pb_Max(1:N_t)
        Pb(2,1:N_t) <= Pb_Max(1:N_t)
        Pb(1,1:N_t) >= Pb_Min(1:N_t)
        Pb(2,1:N_t) >= Pb_Min(1:N_t)
        E(1,1:N_t) <= E_Max(1:N_t)
        E(2,1:N_t) <= E_Max(1:N_t)
%         Ppr(1:N_t) <= Ppr_Max(1:N_t)

        % system power balance
        for t_index = 1:12
            Parameter.alpha * OPTIONS.P_L_TIME(1,t_index) + Ppr(t_index) == Pg(1,t_index) + Pb(1,t_index)
            (1-Parameter.alpha) * OPTIONS.P_L_TIME(1,t_index) == Pg(2,t_index) + Pb(2,t_index)
        end
        
        % ESM output power and the capacity constraints        
        2 - Pb(1,1) == E(1,1)
        2 - Pb(2,1) == E(2,1)
        for t_index = 1:11
            E(1,t_index) - Pb(1,t_index+1) == E(1,t_index+1)
            E(2,t_index) - Pb(2,t_index+1) == E(2,t_index+1)
        end
        sum((Ppr(1:N_t)/2.2e-3).^(1/3)) >= OPTIONS.Distance
cvx_end

objval = cvx_optval;

% obtain the largrange multiplier

cvx_begin
    variable lambda_d nonnegative
    variable lambda_b1(1, N_t) nonnegative
    variable lambda_b2(1, N_t) nonnegative
    maximize(    sum(  Parameter.G(1,1)* power(Pg(1,1:N_t),2)  + Parameter.G(1,2)*Pg(1,1:N_t) + Parameter.G(1,3)*Pg_constant(1:N_t) ...
            + Parameter.G(2,1)* power(Pg(2,1:N_t),2)  + Parameter.G(2,2)*Pg(2,1:N_t) + Parameter.G(2,3)*Pg_constant(1:N_t) ...
            + Parameter.E(1,1)* power(Pb(1,1:N_t),2) + Parameter.E(1,1)* power(Pb(2,1:N_t),2) ) ...
            + lambda_d*( OPTIONS.Distance - sum((Ppr(1:N_t)/2.2e-3).^(1/3)) ) + sum( lambda_b1(1:N_t) .* ( Pg(1:N_t) + Pb(1:N_t) - Parameter.alpha * OPTIONS.P_L_TIME(1:N_t) - Ppr(1:N_t) ) ...
            + lambda_b2(1:N_t) .* ( Pg(2,1:N_t) + Pb(2,1:N_t) - (1 - Parameter.alpha) * OPTIONS.P_L_TIME(1:N_t) ) ) )
    subject to  
            lambda_d <= 10 
            lambda_b1(1:N_t) <= 5
            lambda_b2(1:N_t) <= 5
            
cvx_end


number_of_startup_shutdown = sum(abs(delta(1:N_g,2:N_t) - delta(1:N_g,1:N_t-1))); 
% obtain the largrange function
UB = objval + sum( C_ss * number_of_startup_shutdown) + sum( Parameter.E(1,2)* power(Pb(1,1:N_t),1) + Parameter.E(1,2)* power(Pb(2,1:N_t),1) ) ...
            + lambda_d*( OPTIONS.Distance - sum((Ppr(1:N_t)/2.2e-3).^(1/3)) ) + sum( lambda_b1(1:N_t) .* ( Pg(1:N_t) + Pb(1:N_t) - Parameter.alpha * OPTIONS.P_L_TIME(1:N_t) - Ppr(1:N_t) ) ...
            + lambda_b2(1:N_t) .* ( Pg(2,1:N_t) + Pb(2,1:N_t) - (1-Parameter.alpha) * OPTIONS.P_L_TIME(1:N_t) ) );

% LB =  master_problem(Pg, Pb, sub_delta, Ppr, lambda_d, lambda_b1, lambda_b2 );

% error_primal_dual = UB-LB;

%% FIGURE PLOT
figure
plot(Ppr,'linewidth',1.5);
hold on
plot(OPTIONS.P_L_TIME(1,:),'linewidth',1.5);
hold on
% plot(Ppr(1,1:N_t)+P_L_TIME(1,1:N_t),'k','linewidth',2);
hold on
plot(Pb(1,1:N_t),'linewidth',1.5);
hold on
plot(Pb(2,1:N_t),'linewidth',1.5);
hold on
plot(Pg(1,1:N_t),'linewidth',1.5);
hold on
plot(Pg(2,1:N_t),'linewidth',1.5);
hold on
% plot(Pb(1,1:N_t)+Pg(1,1:N_t),'r');
% ylim([0 5]);

legend('P_{PR}','P_{L}','P_{b_1}','P_{b_2}','P_{g_1}','P_{g_2}');
% legend('P_{PR}','P_{L}','P_{B1}','P_{B2}','P_G');
% end
end

%% NORMAL OPERATION COST FUNCTION
function total_minicost = minicost_in_normal(P_demand)
    PG(1,:) = linspace(1,12,40);

    total_minicost = 0;
    for index = 1:1:12
    %     cost(index,:) = 320 + 10 *PG(1,:) + 22 *PG(1,:).^2 + 340 + 20*(P_demand(index)-PG(1,:)) + 32*(P_demand(index)-PG(1,:)).^2;
        cost(index,:) = 320 + 10 *PG(1,:) + 18 *PG(1,:).^2 + 340 + 20*(P_demand(index)-PG(1,:)) + 32*(P_demand(index)-PG(1,:)).^2;
    %     figure
    %     plot(PG1,cost);
        mini_PG(1,index) = PG(1,find(cost(index,:) == min( cost(index,:))) );
        mini_PG(2,index) = P_demand(index) - mini_PG(1,index);
        total_minicost = total_minicost + min(cost(index,:));
    end
end