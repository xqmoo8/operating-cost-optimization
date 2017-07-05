function cost_optimization_on()
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
OPTIONS.P_L_Scale = [0.5 0.6 0.8 0.8 0.7 0.6 0.4 0.35 0.3 0.33 0.4 0.5]; % P_Generater
% OPTIONS.P_V_SV = 0.8;
P_L_TIME = sum(OPTIONS.P_L.'* OPTIONS.P_L_Scale, 1);

P_prop = 2.2e-3*(OPTIONS.Distance/12).^3;
% Population_one.velocity(timeindex) = (Population_one.chrom(timeindex)/2.2e-3).^(1/3);
% P_load = sum(OPTIONS.P_L.'* OPTIONS.P_L_Scale, 1);
P_demand = P_prop + P_L_TIME;

N_e = 2; N_g = 2; N_t = OPTIONS.Operation_Time;

% generator function parameters
Parameter.G(1,1:3) = [13.5 10 300];
Parameter.G(2,1:3) = [6 30 250];
Parameter.E(1,1:3) = [150 2 0];
alpha = 0.6;

Pg_Max(1:N_t) = 8;
Ppr_Max(1:N_t) = 12;
Pb_Max(1:N_t) = 2;
Pb_Min(1:N_t) = -2;
E_Max(1:N_t) = 2;
Pg_constant(1:N_t) = 1;

Pb(1,1:N_t) = 0;
E(2,1:N_t) = 2;
% Ppr(1:N_t) = OPTIONS.P_pr_avg;


%% ISLAND PART #2# SOLVER
cvx_begin
    variable Ppr(1,N_t) nonnegative
    variable Pb(2,N_t)
    variable E(2,N_t) nonnegative
%     variable delta(N_g,N_t) binary
    variable Pg(2,N_t) nonnegative
    minimize( sum(  Parameter.G(1,1)* power(Pg(1,1:N_t),2)  + Parameter.G(1,2)*Pg(1,1:N_t) + Parameter.G(1,3)*Pg_constant(1:N_t) ...
                    + Parameter.G(2,1)* power(Pg(2,1:N_t),2)  + Parameter.G(2,2)*Pg(2,1:N_t) + Parameter.G(2,3)*Pg_constant(1:N_t) ...
                    + Parameter.E(1,1)* power(Pb(1,1:N_t),1) + Parameter.E(1,2)* power(Pb(1,1:N_t),2) + Parameter.E(1,1)* power(Pb(2,1:N_t),1) + Parameter.E(1,2)* power(Pb(2,1:N_t),2) ) )
%                     + Parameter.E(1,1)* power(Pb(1,1:N_t),1) + Parameter.E(1,1)* power(Pb(2,1:N_t),1) ) )
    subject to   
        % the range constraints of all the variables
        Pg(1,1:N_t) <= Pg_Max(1:N_t)
        Pg(2,1:N_t) <= Pg_Max(1:N_t)
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
            alpha * P_L_TIME(1,t_index) + Ppr(t_index) == Pg(1,t_index) + Pb(1,t_index)
            (1-alpha) * P_L_TIME(1,t_index) == Pg(2,t_index) + Pb(2,t_index)
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

%% FIGURE PLOT
figure
plot(Ppr,'linewidth',1.5);
hold on
plot(P_L_TIME(1,:),'linewidth',1.5);
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