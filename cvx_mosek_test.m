cvx_solver mosek
cvx_begin
    variable Pb(2,5)
    variable Pg(2,5) nonnegative
    minimize( sum(  power(Pg(1,1:5),2) +  power(Pb(1,1:5),1)  ,2))
    subject to
        Pg(1,1:5) <= ones(1,5) * 5
        Pb(1,1:5) <= ones(1,5) * 1
        Pg(1,1:5) >= ones(1,5) * 0
        Pb(1,1:5) >= ones(1,5) * -1
        
        for index = 1:5
            Pg(1,index) + Pb(2,1:5) == [1 3 4 2 5.2]
        end
            
cvx_end

for index_t = 1:OPTIONS.N_t-1
    startup_test(1,index_t) = delta(1,index_t+1) - delta(1,index_t) ;
    startup_test(2,index_t) = delta(2,index_t+1) - delta(2,index_t) ;
end

for index_t = 1:OPTIONS.N_t-OPTIONS.Tmin-1 % only for OPTIONS.Tmin=2
    up_time(1,index_t) = delta(1,index_t+1) + delta(1,index_t+2) 
    up_time(2,index_t) = delta(2,index_t+1) + delta(2,index_t+2) 
end

OPTIONS.Tmin*startup_test(1,:)
OPTIONS.Tmin*startup_test(2,:)

y = size(find(delta(:,2:OPTIONS.N_t) - delta(:,1:OPTIONS.N_t-1)==1),2);
sum( sum(  Parameter.G(1,3)*sum(delta(1:OPTIONS.N_g,1:OPTIONS.N_t),2)  ...
               + sum(Parameter.E(1,2)*Pb(1:2,1:OPTIONS.N_t),2) ,1) + cvx_optval_sub ) + Parameter.C_ss*y 
           
sum( sum(  Parameter.G(1,3)*sum(delta(1:OPTIONS.N_g,1:OPTIONS.N_t),2)   ...
               + sum(Parameter.E(1,2)*Pb(1:2,1:OPTIONS.N_t),2) ,1)+ cvx_optval + sum(Parameter.C_ss*y,2))

benders_cut = cvx_optval_sub + lambda_Pe(1,:)*Pb_m(1,:).' + lambda_Pe(2,:)*Pb_m(2,:).' - lambda_Pe(1,:)*Pb(1,:).' - lambda_Pe(2,:)*Pb(2,:).' ...
           + lambda_delta(1,:)*delta_g(1,:).' + lambda_delta(2,:)*delta_g(2,:).' - lambda_delta(1,:)*delta(1,:).' - lambda_delta(2,:)*delta(2,:).'...
           + sum(Redundent_switch_m(1,1)*lambda_Px(1,:) + Redundent_switch_m(1,2)*lambda_Px(1,:) - Redundent_switch(1,1)*lambda_Px(1,:) - Redundent_switch(1,2)*lambda_Px(1,:)...
           + Redundent_switch_m(1,3)*lambda_Sy(1,:) + Redundent_switch_m(1,4)*lambda_Sy(1,:) - ~Redundent_switch(1,1)*lambda_Sy(1,:) - ~Redundent_switch(1,2)*lambda_Sy(1,:))

%% FIGURE PLOT
figure
hold on
bar([ Pg(1,1:OPTIONS.N_t); Pg(2,1:OPTIONS.N_t); Pb(1,1:OPTIONS.N_t); Pb(2,1:OPTIONS.N_t)].','stacked');
plot(Ppr(1:OPTIONS.N_t),'linewidth',1.5);
plot(OPTIONS.P_L_TIME_off(1,1:OPTIONS.N_t),'linewidth',1.5);

xlim([0 OPTIONS.N_t+1]);
legend('P_{G1}','P_{G2}','P_{pr}','P_{l}','Orientation','horizontal');
ylabel('Active Power (MW)');
xlabel('Time (hours)');

legend('P_{g_1}','P_{g_2}','P_{b_1}','P_{b_2}','P_{PR}','P_{L}');
hold off