function [objval, Pg, Pb,  Ppr, lambda_d, lambda_b1, lambda_b2] = cost_optimization_offline(algorithm_mode)
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

if ~exist('algorithm_mode', 'var')
% defaut: offline algorithm mode
    algorithm_mode =1;
end

if algorithm_mode ==1
    % the load demand without random feature
    P_L_Scale = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline'); 
elseif algorithm_mode ==0
    % the load demand with random feature
    P_L_Scale = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline') + 0.1*rand(1, 25);
end

OPTIONS.P_L_TIME = sum(OPTIONS.P_L.'* P_L_Scale, 1);

N_e = OPTIONS.N_e;
N_g = OPTIONS.N_g;
N_t = OPTIONS.N_t;
C_ss = Parameter.C_ss;

Pg_Max(1) = OPTIONS.Pg_Max(1);
Pg_Min(1) = OPTIONS.Pg_Min(1);
Pg_Max(2) = OPTIONS.Pg_Max(2);
Pg_Min(2) = OPTIONS.Pg_Min(2);

Ppr_Max(1) = OPTIONS.Ppr_Max(1);
Pb_Max(1) = OPTIONS.Pb_Max(1);
Pb_Min(1) = OPTIONS.Pb_Min(1);
E_Max(1) = OPTIONS.E_Max(1);
Pg_constant(1) = OPTIONS.Pg_constant(1);
R_G = Parameter.R_G;

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
    minimize( sum(  Parameter.G(1,1)* power(Pg(1,1:N_t),2)  + Parameter.G(1,2)*Pg(1,1:N_t) + Parameter.G(1,3)*Pg_constant*ones(1,N_t) ...
                    + Parameter.E(1,1)* power(Pb(1,1:N_t),2)  ) )
%                     + Parameter.G(2,1)* power(Pg(2,1:N_t),2)  + Parameter.G(2,2)*Pg(2,1:N_t) + Parameter.G(2,3)*Pg_constant(1:N_t) ...
    subject to         
        % the range constraints of all the variables
        Pg(1,1:N_t) <= delta(1:N_t)*Pg_Max(1)
        Pg(2,1:N_t) <= delta(1:N_t)*Pg_Max(2)
%         Pg(1,1:N_t) >= delta(1,N_t).*Pg_Min(1:N_t)
%         Pg(2,1:N_t) >= delta(2,N_t).*Pg_Min(2:N_t)
        
        for t_index = 1:OPTIONS.N_t-1
            abs( Pg(1,t_index) - Pg(1,t_index+1) ) * delta(1,t_index) * delta(1,t_index+1) <= R_G
            abs( Pg(2,t_index) - Pg(2,t_index+1) ) * delta(2,t_index) * delta(2,t_index+1) <= R_G
        end
        
        Ppr(1:N_t) <= Ppr_Max(1)*ones(1,N_t)
        Pb(1,1:N_t) <= Pb_Max(1)*ones(1,N_t)
        Pb(2,1:N_t) <= Pb_Max(1)*ones(1,N_t)
        Pb(1,1:N_t) >= Pb_Min(1)*ones(1,N_t)
        Pb(2,1:N_t) >= Pb_Min(1)*ones(1,N_t)
        
        E(1,1:N_t) <= E_Max(1)*ones(1,N_t)
        E(2,1:N_t) <= E_Max(1)*ones(1,N_t)
%         Ppr(1:N_t) <= Ppr_Max(1:N_t)

        % system power balance
        for t_index = 1:OPTIONS.N_t
            Parameter.alpha * OPTIONS.P_L_TIME(1,t_index) + Ppr(t_index) == Pg(1,t_index) + Pb(1,t_index)
            (1-Parameter.alpha) * OPTIONS.P_L_TIME(1,t_index) == Pg(2,t_index) + Pb(2,t_index)
        end
        
        % ESM output power and the capacity constraints
        2 - Pb(1,1) == E(1,1)
        2 - Pb(2,1) == E(2,1)
        for t_index = 1:OPTIONS.N_t-1
            E(1,t_index) - Pb(1,t_index+1) == E(1,t_index+1)
            E(2,t_index) - Pb(2,t_index+1) == E(2,t_index+1)
        end
        sum((Ppr(1:N_t)/2.2e-3).^(1/3)) >= OPTIONS.Distance
cvx_end

objval = cvx_optval;


% %% FIGURE PLOT
% figure
% hold on
% bar([ Pg(1,1:N_t); Pg(2,1:N_t); Pb(1,1:N_t); Pb(2,1:N_t)].','stacked');
% plot(Ppr(1:OPTIONS.N_t),'linewidth',1.5);
% plot(OPTIONS.P_L_TIME(1,1:OPTIONS.N_t),'linewidth',1.5);
% 
% xlim([0 25]);
% % ylim([0 10]);
% % plot([1 12], [P_prop P_prop], 'linewidth', 2, 'Color',[0.0,0.6,0.9]);
% % plot(1:12, P_load, 'linewidth', 2, 'Color',[0.67,0,1.0]);
% legend('P_{G1}','P_{G2}','P_{pr}','P_{l}','Orientation','horizontal');
% ylabel('Active Power (MW)');
% xlabel('Time (hours)');
% 
% legend('P_{g_1}','P_{g_2}','P_{b_1}','P_{b_2}','P_{PR}','P_{L}');
% hold off

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