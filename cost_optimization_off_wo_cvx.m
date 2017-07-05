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

OPTIONS.N_e = 1;
OPTIONS.N_g = 1; 
OPTIONS.N_t = OPTIONS.Operation_Time;


% generator function parameters
Parameter.G(1,1:3) = [13.5 10 300];
Parameter.G(2,1:3) = [6 30 250];
Parameter.E(1,1:3) = [2 2 0];
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

Pb(1,1:N_t) = 1;
E(2,1:N_t) = 2;

delta(1:2,1:N_t) = 1;
switch_PS(1:2) = 1;
error_primal_dual = 10;

step_size = 0.0001;

lambda_d = 1.2;
lambda_b(1:N_t) = 0.1;

Pg(1:N_t) = OPTIONS.P_pr_avg+4;
Ppr(1:N_t) = OPTIONS.P_pr_avg;
Pb(1:N_t) = 0;

for iteration_index = 1:1000
    %% variable update, partial of largrange function
    % P_G 
    for g_index = 1:N_g
        for t_index = 1:N_t
            Pg_t = Pg(g_index, t_index);
            lambda_bt = lambda_b( t_index );
%             d_Pg_t = 2*Parameter.G(1,1)*Pg_t + Parameter.G(1,2) + lambda_bt ;
            d_Pg_t = 2*Parameter.G(1,1)*Pg_t + Parameter.G(1,2) -lambda_d*( (Pg(t_index) + Pb(t_index) - Parameter.alpha * OPTIONS.P_L_TIME(t_index)) /2.2e-3)^(-2/3) ;
            Pg(g_index, t_index) = Pg(g_index, t_index) - step_size*d_Pg_t;
        end
    end

    % P_E
    for e_index = 1:N_e
        for t_index = 1:N_t
            Pb_t = Pb(e_index, t_index);
            lambda_bt = lambda_b( t_index );
%             d_Pb_t = 2*Parameter.E(1,1)*Pb_t + lambda_bt ;
            d_Pb_t = 2*Parameter.E(1,1)*Pb_t -lambda_d*( (Pg(t_index) + Pb(t_index) - Parameter.alpha * OPTIONS.P_L_TIME(t_index)) /2.2e-3)^(-2/3) 
            Pb(e_index, t_index) = Pb(e_index, t_index) - step_size*d_Pb_t;
        end
    end
    
    Ppr(1:N_t) = Pg(1:N_t) + Pb(1:N_t) - Parameter.alpha * OPTIONS.P_L_TIME(1:N_t);

%     % P_pr
%     for t_index = 1:N_t
%         Ppr_t = Ppr(t_index);
%         lambda_bt = lambda_b( t_index );
%         d_Ppr_t = -lambda_d*(Ppr_t/2.2e-3)^(-2/3) - lambda_bt ;  
%         Ppr(t_index) = Ppr(t_index) - step_size*d_Ppr_t;
%     end

    % E_b
    for t_index = 1:N_t
        Ppr_t = E_b(t_index);
        lambda_bt = lambda_b( t_index );
        d_Ppr_t = -lambda_d*(Ppr_t/2.2e-3)^(-2/3) - lambda_bt ;  
        Ppr(t_index) = Ppr(t_index) - step_size*d_Ppr_t;
    end

    primal_obj(iteration_index) = sum(  Parameter.G(1,1)* power(Pg(1,1:N_t),2)  + Parameter.G(1,2)*Pg(1,1:N_t) + Parameter.G(1,3)*Pg_constant(1:N_t) ...
                + Parameter.E(1,1)* power(Pb(1,1:N_t),2) + lambda_d*( OPTIONS.Distance - sum((Ppr(1:N_t)/2.2e-3).^(1/3)) ) ...
                + sum( lambda_b(1:N_t) .* ( Pg(1:N_t) + Pb(1:N_t) - Parameter.alpha * OPTIONS.P_L_TIME(1:N_t) - Ppr(1:N_t) ) ) );

    % lambda_d
    lambda_d = lambda_d + step_size * ( OPTIONS.Distance - sum((Ppr(1:N_t)/2.2e-3).^(1/3)) ) ;

    % lambda_b
    for t_index = 1:N_t
    %     lambda_b_t = lambda_b( t_index );
        d_lamda_b_t = ( Pg(1,t_index) + Pb(1,t_index) - (1 - Parameter.alpha) * OPTIONS.P_L_TIME(t_index) - Ppr(t_index)) ;
        lambda_b( t_index ) = lambda_b( t_index ) + step_size*d_lamda_b_t;
    end

    dual_obj(iteration_index) = sum(  Parameter.G(1,1)* power(Pg(1,1:N_t),2)  + Parameter.G(1,2)*Pg(1,1:N_t) + Parameter.G(1,3)*Pg_constant(1:N_t) ...
                + Parameter.E(1,1)* power(Pb(1,1:N_t),2) + lambda_d*( OPTIONS.Distance - sum((Ppr(1:N_t)/2.2e-3).^(1/3)) ) ...
                + sum( lambda_b(1:N_t) .* ( Pg(1:N_t) + Pb(1:N_t) - Parameter.alpha * OPTIONS.P_L_TIME(1:N_t) - Ppr(1:N_t) ) ) );

    % if error_primal_dual <= 1e-3
    if abs(primal_obj(iteration_index) - dual_obj(iteration_index))< 0.0001
        break;
    end
    
end    
    

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
% plot(Pb(2,1:N_t),'linewidth',1.5);
% hold on
plot(Pg(1,1:N_t),'linewidth',1.5);
hold on
% plot(Pg(2,1:N_t),'linewidth',1.5);
% hold on
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