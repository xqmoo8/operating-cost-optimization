clear;

% global OPTIONS
objval = 0;

OPTIONS.Distance = 250;
OPTIONS.velocity = [17 0];
OPTIONS.N_e = 2;
OPTIONS.N_g = 2; 
OPTIONS.N_t = 24;

OPTIONS.Pg_Max(1) = 8;
OPTIONS.Pg_Min(1) = 1;
OPTIONS.Pg_Max(2) = 4;
OPTIONS.Pg_Min(2) = 0.5;

OPTIONS.Ppr_Max = 12;
OPTIONS.Pb_Max(1) = 1;
OPTIONS.Pb_Min(1) = -1;
OPTIONS.E_Max(1) = 3;

OPTIONS.P_L = [2.7 0.9]; % P_Generater
OPTIONS.P_L_Scale = [0.5 0.6 0.8 0.82 0.7 0.6 0.4 0.35 0.3 0.33 0.4 0.5 0.4]; 

% the load demand without random feature
P_L_Scale_off = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline'); 
% the load demand with random feature
P_L_Scale_on = interp1(1:13,OPTIONS.P_L_Scale,1:0.5:13,'spline') + 0.1*rand(1, 25);

OPTIONS.P_L_TIME_off = sum(OPTIONS.P_L.'* P_L_Scale_off(:,1:OPTIONS.N_t), 1);
OPTIONS.P_L_TIME_on  = sum(OPTIONS.P_L.'* P_L_Scale_on(:,1:OPTIONS.N_t), 1);

Coupled_load(1, :) = 0.05*sum(OPTIONS.P_L_Scale.'*P_L_Scale_off(:,1:OPTIONS.N_t), 1);
Coupled_load(2, :) = 0.05*sum(OPTIONS.P_L_Scale.'*P_L_Scale_off(:,1:OPTIONS.N_t), 1);
Redundent_switch(1,1:2) = 1;

% generator function parameters
Parameter.G(1,1:3) = [13.5 10 300];
Parameter.G(2,1:3) = [6 30 250];
Parameter.E(1,1:3) = [10 10 0];
Parameter.alpha = 0.6;
Parameter.C_ss = 100;
Parameter.R_G = 1;
Parameter.error = 1e-3;

lambda_delta = zeros(1,OPTIONS.N_t);
lambda_Pb = zeros(1,OPTIONS.N_t);
delta_g = ones(1,OPTIONS.N_t);
% Pb = OPTIONS.Pb_Max/2*ones(1,OPTIONS.N_t);

delta(1,1:OPTIONS.N_t ) = 1;
delta(2,1:OPTIONS.N_t ) = 1;

for index_sm = 1:20
%% subproblem
cvx_begin
% cvx_begin quiet 
    variable Ppr(1,OPTIONS.N_t) nonnegative
    variable Pb(2,OPTIONS.N_t)
    variable E(2,OPTIONS.N_t) nonnegative
%     variable delta(N_g,N_t) binary
    variable Pg(2,OPTIONS.N_t) nonnegative
    minimize( sum(  sum(Parameter.G(1,1)* power(Pg(1:2,1:OPTIONS.N_t),2)  + Parameter.G(1,2)*Pg(1:2,1:OPTIONS.N_t) + Parameter.G(1,3)*ones(2,OPTIONS.N_t) ...
                    + Parameter.E(1,1)* power(Pb(1:2,1:OPTIONS.N_t),2) ,2) ,1 ) )
    subject to
        % the range constraints of all the variables
        Pg(1,1:OPTIONS.N_t) <= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1)
        Pg(2,1:OPTIONS.N_t) <= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Max(1)
%         Pg(1,1:OPTIONS.N_t) >= delta(1,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)
%         Pg(2,1:OPTIONS.N_t) >= delta(2,1:OPTIONS.N_t) * OPTIONS.Pg_Min(1)

        Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(1,2:OPTIONS.N_t) -Pg(1,1:OPTIONS.N_t-1) >= -Parameter.R_G
        Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) <= Parameter.R_G
        Pg(2,2:OPTIONS.N_t) -Pg(2,1:OPTIONS.N_t-1) >= -Parameter.R_G

        Ppr(1,1:OPTIONS.N_t) <= OPTIONS.Ppr_Max * ones(1,OPTIONS.N_t)
        Pb(1,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
        Pb(1,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
        Pb(2,1:OPTIONS.N_t) <= OPTIONS.Pb_Max * ones(1,OPTIONS.N_t)
        Pb(2,1:OPTIONS.N_t) >= OPTIONS.Pb_Min * ones(1,OPTIONS.N_t)
      
        
        E(1,1:OPTIONS.N_t) <=  OPTIONS.E_Max * ones(1,OPTIONS.N_t)
        E(1,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)
        E(2,1:OPTIONS.N_t) <=  OPTIONS.E_Max * ones(1,OPTIONS.N_t)
        E(2,1:OPTIONS.N_t) >= zeros(1,OPTIONS.N_t)
       
        % system power balance
        for t_index = 1:OPTIONS.N_t
            Redundent_switch*Coupled_load(:,t_index) +  Parameter.alpha * OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == Pg(1,t_index) + Pb(1,t_index)
%              Parameter.alpha * OPTIONS.P_L_TIME_off(1,t_index) + Ppr(1,t_index) == Pg(1,t_index) + Pb(1,t_index)
            ~Redundent_switch*Coupled_load(:,t_index) + (1-Parameter.alpha) * OPTIONS.P_L_TIME_off(1,t_index)  == Pg(2,t_index) + Pb(2,t_index)
        end
        
        % ESM output power and the capacity constraints
        2 - Pb(1,1) == E(1,1)
        2 - Pb(2,1) == E(2,1)
        for t_index = 1:OPTIONS.N_t-1
            E(1,t_index) - Pb(1,t_index+1) == E(1,t_index+1)
            E(2,t_index) - Pb(2,t_index+1) == E(2,t_index+1)
        end
        
        sum((Ppr(1:OPTIONS.N_t)./2.2e-3).^(1/3)) >= OPTIONS.Distance
cvx_end

b_delta_bound = (Pg(2,1:OPTIONS.N_t).')./OPTIONS.Pg_Min(2);


%% FIGURE PLOT
figure
hold on
bar([ Pg(1,1:OPTIONS.N_t); Pg(2,1:OPTIONS.N_t); Pb(1,1:OPTIONS.N_t); Pb(2,1:OPTIONS.N_t)].','stacked');
plot(Ppr(1:OPTIONS.N_t),'linewidth',1.5);
plot(OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t),'linewidth',1.5);

xlim([0 OPTIONS.N_t+1]);
% ylim([0 10]);
% plot([1 12], [P_prop P_prop], 'linewidth', 2, 'Color',[0.0,0.6,0.9]);
% plot(1:12, P_load, 'linewidth', 2, 'Color',[0.67,0,1.0]);
legend('P_{G1}','P_{G2}','P_{pr}','P_{l}','Orientation','horizontal');
ylabel('Active Power (MW)');
xlabel('Time (hours)');

legend('P_{g_1}','P_{g_2}','P_{b_1}','P_{b_2}','P_{PR}','P_{L}');
hold off


%% dual variable and upperbound
lambda_Pb_1 = -Parameter.G(1,1)*Pg(1,1:OPTIONS.N_t) - Parameter.G(1,2) ;
lambda_delta = zeros(1,OPTIONS.N_t);
lambda_Pb_2 = -Parameter.G(1,1)*Pg(2,1:OPTIONS.N_t) - Parameter.G(1,2) ;
lambda_delta = zeros(1,OPTIONS.N_t);

y = size(find(delta(1:2,2:OPTIONS.N_t) - delta(1:2,1:OPTIONS.N_t-1)==1),2);
objval_upperbound= cvx_optval + y*Parameter.C_ss(1) + sum(sum( Parameter.E(1,2)*Pb(1:2,1:OPTIONS.N_t) + Parameter.E(1,3)*ones(2,OPTIONS.N_t),1),2) ;

cvx_optval_2 =  sum(Parameter.G(1,1)* power(Pg(2,1:OPTIONS.N_t),2)  + Parameter.G(1,2)*Pg(2,1:OPTIONS.N_t) + Parameter.G(1,3)*ones(1,OPTIONS.N_t) ...
                    + Parameter.E(1,1)* power(Pb(2,1:OPTIONS.N_t),2) );
y_2 = size(find(delta(2,2:OPTIONS.N_t) - delta(1,1:OPTIONS.N_t-1)==1),2);
objval_upperbound_2= cvx_optval_2 + y_2*Parameter.C_ss(1) + sum( Parameter.E(1,2)*Pb(2,1:OPTIONS.N_t) + Parameter.E(1,3)*ones(1,OPTIONS.N_t)) ;

%% build master problem
%% delta(1,Nt); y(1,Nt-1); Pb(1,Nt); E(1,Nt); mu(1); redundant_switch(1,2) total number = 4*Nt + 2
% delta and Pg_min
A_delta_bound = zeros(OPTIONS.N_t, 4*OPTIONS.N_t)
for index_x = 1:OPTIONS.N_t
    A_delta_bound(index_x,index_x) = 1;
end
b_delta_bound = (Pg(2,1:OPTIONS.N_t).')./OPTIONS.Pg_Min(2);

% delta and y
A_delta = zeros(OPTIONS.N_t-1, 4*OPTIONS.N_t)
for index_x = 1:OPTIONS.N_t-1
    if index_x <= OPTIONS.N_t
        A_delta(index_x,index_x) = -1;
        A_delta(index_x,index_x+1) = 1;
    end
    A_delta(index_x,OPTIONS.N_t+index_x) = -1;
end
b_y = zeros(OPTIONS.N_t-1,1);

% benders cuts
A_benders = zeros(1, 4*OPTIONS.N_t);
for index_x = 1:OPTIONS.N_t
    A_benders(1, index_x) = - lambda_delta(1,index_x); 
end
for index_x = 2*OPTIONS.N_t:3*OPTIONS.N_t-1
    A_benders(1, index_x) = - lambda_Pb_2(1,index_x-2*OPTIONS.N_t+1);
end
A_benders(1, 4*OPTIONS.N_t) = -1; 

% do not have the value of last loop
if index_sm == 1
    b_benders =  - cvx_optval_2 + sum(lambda_delta*delta(2,:)') + sum(lambda_Pb_2*Pb(2,:)');
else
    b_benders =  - cvx_optval_2 + sum(lambda_delta*delta(2,:)') + sum(lambda_Pb_2*Pb(2,:)');
end

% ESM inequality constraint
A_E_start = [zeros(1,2*OPTIONS.N_t-1) -1 zeros(1,OPTIONS.N_t-1) -1 zeros(1,OPTIONS.N_t)];
b_E_start = -OPTIONS.E_Max;
% A_E_end = [zeros(1,3*OPTIONS.N_t-1) -1 zeros(1,OPTIONS.N_t-2) 1 0];
% b_E_end = OPTIONS.E_Max;

% A = [A_delta; A_E_start; A_delta_bound];
% b = [b_y; b_E_start; b_delta_bound];

A = [A_delta; A_benders; A_E_start; A_delta_bound];
b = [b_y; b_benders; b_E_start; b_delta_bound];

% ESM equality constraint
Aeq_E = zeros(OPTIONS.N_t-1, 4*OPTIONS.N_t);
for index_x = 1: OPTIONS.N_t-1
    Aeq_E(index_x, 2*OPTIONS.N_t+index_x) = - 1; 
    Aeq_E(index_x, 3*OPTIONS.N_t+index_x-1) =   1; 
    Aeq_E(index_x, 3*OPTIONS.N_t+index_x) = - 1; 
end
beq_E = zeros(OPTIONS.N_t-1,1);

% balance equality constraint
Aeq_B = zeros(OPTIONS.N_t,4*OPTIONS.N_t);
for index_x = 1: OPTIONS.N_t
    Aeq_B(index_x, index_x) = Pg(2,index_x);
    Aeq_B(index_x, 2*OPTIONS.N_t+index_x-1) = 1;
end
beq_B = (Parameter.alpha * OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t)).';

Aeq = [Aeq_E; Aeq_B];
beq = [beq_E; beq_B];

% upper and lower bound
lb = [zeros(1,OPTIONS.N_t) zeros(1,OPTIONS.N_t-1) -OPTIONS.Pb_Max*ones(1,OPTIONS.N_t) zeros(1,OPTIONS.N_t)                0  ];
ub = [ones(1,OPTIONS.N_t)  ones(1,OPTIONS.N_t-1)  OPTIONS.Pb_Max*ones(1,OPTIONS.N_t)  OPTIONS.E_Max*ones(1,OPTIONS.N_t)   inf ];

f = [zeros(1,OPTIONS.N_t) ones(1,OPTIONS.N_t-1) Parameter.E(2)*ones(1,OPTIONS.N_t)  zeros(1,OPTIONS.N_t) 1];

intcon = [ 1: 2*OPTIONS.N_t-1 ];
options = optimoptions('intlinprog','Display','iter');
x = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub, options);

objval_lowerbound = f(1,OPTIONS.N_t+1:2*OPTIONS.N_t-1)*x(OPTIONS.N_t+1:2*OPTIONS.N_t-1,1) + f(1,2*OPTIONS.N_t:3*OPTIONS.N_t-1)*x(2*OPTIONS.N_t:3*OPTIONS.N_t-1,1) + x(end);

delta =  x(1:OPTIONS.N_t,1).';

%% FIGURE PLOT
figure
hold on
bar([ Pg(1,1:OPTIONS.N_t); Pb(1,1:OPTIONS.N_t)].','stacked');
plot(Ppr(1:OPTIONS.N_t),'linewidth',1.5);
plot(OPTIONS.P_L_TIME_on(1,1:OPTIONS.N_t),'linewidth',1.5);

xlim([0 OPTIONS.N_t+1]);
% ylim([0 10]);
% plot([1 12], [P_prop P_prop], 'linewidth', 2, 'Color',[0.0,0.6,0.9]);
% plot(1:12, P_load, 'linewidth', 2, 'Color',[0.67,0,1.0]);
legend('P_{G1}','P_{G2}','P_{pr}','P_{l}','Orientation','horizontal');
ylabel('Active Power (MW)');
xlabel('Time (hours)');

legend('P_{g_1}','P_{b_1}','P_{PR}','P_{L}');
hold off

if abs(objval_upperbound-objval_lowerbound) < 1e-2
    break;
end
objval_upperbound-objval_lowerbound
end
         