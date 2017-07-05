% benders_decomposition_for_OPMS
clear
global OPTIONS Parameter

OPTIONS.Distance = 300;
OPTIONS.velocity = [17 0];
OPTIONS.N_e = 2;
OPTIONS.N_g = 2; 
OPTIONS.N_t = 24;

OPTIONS.velocity_avg = OPTIONS.Distance/OPTIONS.N_t;
OPTIONS.P_pr_avg = (OPTIONS.velocity_avg).^3*2.2e-3;
OPTIONS.Delta_P_pr = 2;

OPTIONS.P_L = [2.7 0.9]; % P_Generater
OPTIONS.P_L_Scale = [0.5 0.6 0.8 0.82 0.7 0.6 0.4 0.35 0.3 0.33 0.4 0.5 0.4]; 

% the load demand without random feature
P_L_Scale_off = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline'); 
% the load demand with random feature
P_L_Scale_on = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline') + 0.1*rand(1, 25);

OPTIONS.P_L_TIME_off = sum(OPTIONS.P_L.'* P_L_Scale_off(:,1:OPTIONS.N_t), 1);
OPTIONS.P_L_TIME_on  = sum(OPTIONS.P_L.'* P_L_Scale_on(:,1:OPTIONS.N_t), 1);
OPTIONS.Delta_Load = OPTIONS.P_L_TIME_on - mean( OPTIONS.P_L_TIME_off );

OPTIONS.P_prop = 2.2e-3*(OPTIONS.Distance/12).^3;
% OPTIONS.P_demand = OPTIONS.P_prop + OPTIONS.P_L_TIME;NS.error = 2e-3;

OPTIONS.Pg_Max(1) = 8;
OPTIONS.Pg_Min(1) = 4;
OPTIONS.Pg_Max(2) = 4;
OPTIONS.Pg_Min(2) = 2;

OPTIONS.Ppr_Max = 12;
OPTIONS.Pb_Max(1) = 1;
OPTIONS.Pb_Min(1) = -1;
OPTIONS.E_Max(1) = 2;
OPTIONS.Pg_constant(1:OPTIONS.N_t) = 1;

% generator function parameters
Parameter.G(1,1:3) = [13.5 10 300];
Parameter.G(2,1:3) = [6 30 250];
Parameter.E(1,1:3) = [150 10 0];
Parameter.alpha = 0.6;
Parameter.C_ss = 100;
Parameter.R_G = 1;
Parameter.error = 1e-3;

% % variables in subproblem
% obj_up = inf;
% PG(2:OPTIONS.N_t) = 0;
% PPR(1:OPTIONS.N_t) = 0;
% PB(1:OPTIONS.N_t) = 0;
% % variables in master problem
% obj_down = -inf;
% SW(1:2) = 1;
% Delta_g(2:OPTIONS.N_t) = 2;

objval = 0;

for index_loop = 1:OPTIONS.N_t
%     [Pg, Pb, delta, Ppr, lambda_d, lambda_b1, lambda_b2] = cost_optimization_off_w_cvx();
%     objval = objval + cost_optimization_online(index_loop);
    cost_optimization_offline();
    

    if obj_up - obj_down < Parameter.error
        break;
    end
end


%% FIGURE PLOT
figure
hold on
bar([ Pg(1,1:OPTIONS.N_t); Pg(2,1:OPTIONS.N_t); Pb(1,1:OPTIONS.N_t); Pb(2,1:OPTIONS.N_t)].','stacked');
plot(Ppr(1:OPTIONS.N_t),'linewidth',1.5);
plot(OPTIONS.P_L_TIME(1,1:OPTIONS.N_t),'linewidth',1.5);

xlim([0 25]);
% ylim([0 10]);
% plot([1 12], [P_prop P_prop], 'linewidth', 2, 'Color',[0.0,0.6,0.9]);
% plot(1:12, P_load, 'linewidth', 2, 'Color',[0.67,0,1.0]);
legend('P_{G1}','P_{G2}','P_{pr}','P_{l}','Orientation','horizontal');
ylabel('Active Power (MW)');
xlabel('Time (hours)');

legend('P_{g_1}','P_{g_2}','P_{b_1}','P_{b_2}','P_{PR}','P_{L}');
hold off