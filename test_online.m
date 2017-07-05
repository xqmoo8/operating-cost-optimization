function [operationg_cost_total ] = operating_cost_online( )
%% global variable
global OPTIONS Parameter

load('load_information');
objval = 0;

% the weight factor for 
H = linspace(2e-6,5e-4,10);
Xi = linspace(0.0,0.8,8);

OPTIONS.Delta_Load = OPTIONS.P_L_TIME_on - mean( OPTIONS.P_L_TIME_off );

lambda_delta = zeros(1,OPTIONS.N_t);
lambda_Pb = zeros(1,OPTIONS.N_t);
delta_g = ones(1,OPTIONS.N_t);
% Pb = OPTIONS.Pb_Max/2*ones(1,OPTIONS.N_t);

delta(1,1:OPTIONS.N_t ) = 1;
delta(2,1:OPTIONS.N_t ) = 1;

% operation_mode =3;
% Redundent_switch(1,1:2) = [0 1];
% cost = cost_optimization( delta, Redundent_switch, operation_mode );

%%  operation_mode:
% 0 (normal wo PPA ESMC) 
% 1 (normal w PPA wo ESMC) 
% 2 (normal wo PPA w ESMC) 
% 3 (normal w PPA w ESMC)

% 4 (Fault wo PPA ESMC) 
% 5 (Fault w PPA wo ESMC) 
% 6 (Fault wo PPA w ESMC)
% 7 (Fault w PPA w ESMC)

for index_H = 1:size(H,2)
    for index_Xi = 1:1
%     for index_Xi = 1:size(Xi,2)
        [operationg_cost, Pg_total, Pb_total, Ppr_total ] = cost_optimization(delta, [1 0], 3, H(index_H), Xi(index_Xi));
        operationg_cost_total(index_H, index_Xi) = operationg_cost; 
    end
end

end


function [operationg_cost, Pg_total, Pb_total, Ppr_total, X_E ] = cost_optimization( delta, Redundent_switch, operation_mode, H, Xi)
global OPTIONS Parameter
%%
    Pg_total = zeros(2,OPTIONS.N_t);
    Pb_total = zeros(2,OPTIONS.N_t);
    E_total = zeros(2,OPTIONS.N_t);
    Ppr_total = zeros(1,OPTIONS.N_t);
    X_E = zeros(2,OPTIONS.N_t);
    operationg_cost = 0;
    
%     X_E(1,1:4) = OPTIONS.Pb_Max;
%     X_E(:,1) = OPTIONS.E_Max;

    rest_velocity_avg = OPTIONS.Distance/OPTIONS.N_t;
    rest_propulsion_power_avg = 2.2e-3*rest_velocity_avg^3;
    rest_distance = OPTIONS.Distance;
    
    for t_index = 1:OPTIONS.N_t
        if t_index>1
          if isnan(Pg_total(1,t_index-1) )
              Pg_total =inf;
              Pb_total =inf;
              Ppr_total =inf;
              cvx_optval =inf;
              return
          else
              Pg_bound(1,1) = Parameter.R_G + Pg_total(1,t_index-1);
              Pg_bound(1,2) =   -Parameter.R_G + Pg_total(1,t_index-1);
              Pg_bound(2,1) = Parameter.R_G + Pg_total(2,t_index-1);
              Pg_bound(2,2) =   -Parameter.R_G + Pg_total(2,t_index-1);
          end
        end

        if isempty( delta(1,t_index) ) || delta(1,t_index)==inf || isnan( delta(1,t_index) )
          return
        end
      
%     cvx_begin
        cvx_begin quiet
            variable Ppr(1) nonnegative
            variable Pb(2)
            variable Pg(2) nonnegative
            minimize( H*sum(  Parameter.G(1:2,1).'* power(Pg(1:2,1),2)  + Parameter.G(1:2,2).'*Pg(1:2,1) + Parameter.G(1:2,3).'*ones(2,1)...
                            + sum( Parameter.E(1,1)* power(Pb(1:2,1),2) + Parameter.E(1,2)* power(Pb(1:2,1),1) ,1) ) + sum(Pb(1:2,1).'*X_E(:,t_index)) )
%             minimize( H*sum(  Parameter.G(1:2,1).'* power(Pg(1:2,1),2)  + Parameter.G(1:2,2).'*Pg(1:2,1) + Parameter.G(1:2,3).'*ones(2,1)...
%                             + sum( Parameter.E(1,1)* power(Pb(1:2,1),2) ,1) ) + sum(Pb.'*X_E(:,t_index)) )
%             minimize( H*sum(  Parameter.G(1:2,1).'* power(Pg(1:2,1),2)  + Parameter.G(1:2,2).'*Pg(1:2,1) + Parameter.G(1:2,3).'*ones(2,1)...
%                              ) + sum(Pb.'*X_E(:,t_index)) )
            subject to
                % the range constraints of all the variables
                Pg(1,1) <= delta(1,t_index) * OPTIONS.Pg_Max(1)
                Pg(2,1) <= delta(2,t_index) * OPTIONS.Pg_Max(2)
        %         Pg(1,1) >= delta(1,1) * OPTIONS.Pg_Min(1)
        %         Pg(2,1) >= delta(2,1) * OPTIONS.Pg_Min(1)
                if t_index ==1
                elseif abs ( delta(1,t_index) -  delta(1,t_index-1) ) == 1
                    Pg(1,1) >= 0
                    Pg(2,1) >= 0
                else 
                    Pg(1,1) <= Pg_bound(1,1)
                    Pg(1,1) >= Pg_bound(1,2)
                    Pg(2,1) <= Pg_bound(2,1)
                    Pg(2,1) >= Pg_bound(2,2)
                end

                Ppr(1,1) <= OPTIONS.Ppr_Max
                Pb(1,1) <= OPTIONS.Pb_Max
                Pb(1,1) >= OPTIONS.Pb_Min
                Pb(2,1) <= OPTIONS.Pb_Max
                Pb(2,1) >= OPTIONS.Pb_Min

                if t_index == 1
                    OPTIONS.E_Max(1) + Pb(1,1) <=  OPTIONS.E_Max * ones(1,OPTIONS.N_t)
                    OPTIONS.E_Max(1) + Pb(1,1) >= zeros(1,OPTIONS.N_t)
                    OPTIONS.E_Max(1) + Pb(2,1) <=  OPTIONS.E_Max * ones(1,OPTIONS.N_t)
                    OPTIONS.E_Max(1) + Pb(2,1) >= zeros(1,OPTIONS.N_t)
                else
                    E_total(1,t_index) + Pb(1,1) <=  OPTIONS.E_Max * ones(1,OPTIONS.N_t)
                    E_total(1,t_index) + Pb(1,1) >= zeros(1,OPTIONS.N_t)
                    E_total(2,t_index) + Pb(2,1) <=  OPTIONS.E_Max * ones(1,OPTIONS.N_t)
                    E_total(2,t_index) + Pb(2,1) >= zeros(1,OPTIONS.N_t)
                end

                % system power balance
                if operation_mode <= 3
                    if operation_mode ==0
                        Ppr(1,1) == OPTIONS.P_pr_avg;
                        Pb(1:2,1) == 0
                    elseif operation_mode ==1
                        ((Ppr(1)/2.2e-3).^(1/3)) >= ((rest_propulsion_power_avg - Xi*OPTIONS.Delta_Load(t_index))/2.2e-3).^(1/3)          
                        Pb(1:2,1) == 0
                    elseif operation_mode ==2                    
                        Ppr(1,1) == OPTIONS.P_pr_avg;
                    else
                        ((Ppr(1)/2.2e-3).^(1/3)) >= ((rest_propulsion_power_avg- Xi*OPTIONS.Delta_Load(t_index))/2.2e-3).^(1/3)
                    end

                    sum(OPTIONS.Coupled_load(:,t_index)) +  OPTIONS.P_L_TIME_on(1,t_index) + Ppr(1) == sum(Pg(1:2)) - sum(Pb(1:2))

                elseif operation_mode <= 7
                    if operation_mode ==4                    
                        Ppr(1,1) == OPTIONS.P_pr_avg;
                        Pb(1:2,1) == 0
                    elseif operation_mode ==5
                        ((Ppr(1)/2.2e-3).^(1/3)) >= ((rest_propulsion_power_avg- Xi*OPTIONS.Delta_Load(t_index))/2.2e-3).^(1/3)                    
                        Pb(1:2,1) == 0
                    elseif operation_mode ==6
                        Ppr(1,1) == OPTIONS.P_pr_avg;
                    else
                        ((Ppr(1)/2.2e-3).^(1/3)) >= ((rest_propulsion_power_avg- Xi*OPTIONS.Delta_Load(t_index))/2.2e-3).^(1/3)
                    end

                    if t_index <= 4
                        sum(OPTIONS.Coupled_load(:,t_index)) +  OPTIONS.P_L_TIME_on(1,t_index) + Ppr(1) == sum(Pg(1:2)) - sum(Pb(1:2))
                    else
                        Redundent_switch*OPTIONS.Coupled_load(:,t_index) +  Parameter.alpha * OPTIONS.P_L_TIME_on(1,t_index) + Ppr(1) == Pg(1) - Pb(1)
            %              Parameter.alpha * OPTIONS.P_L_TIME_on(1,t_index) + Ppr(1,t_index) == Pg(1,t_index) + Pb(1,t_index)
                        ~Redundent_switch*OPTIONS.Coupled_load(:,t_index) + (1-Parameter.alpha) * OPTIONS.P_L_TIME_on(1,t_index)  == Pg(2) - Pb(2)
                    end
                end
        cvx_end

        Pg_total(1:2,t_index) = Pg;
        Pb_total(1:2,t_index) = Pb;
        X_E(:,t_index+1) = X_E(:,t_index) + Pb;
        if t_index == 1
            E_total(:,t_index+1) = OPTIONS.E_Max(1) + Pb;
        else
            E_total(:,t_index+1) = E_total(:,t_index) + Pb;
        end
%         E_total(1:2,t_index) = E;
        Ppr_total(1,t_index) = Ppr;
        rest_distance = rest_distance - (Ppr_total(1,t_index)/2.2e-3).^(1/3);
        rest_velocity_avg = rest_distance/(OPTIONS.N_t - t_index);
        rest_propulsion_power_avg = 2.2e-3*rest_velocity_avg^3;
        
%         operationg_cost = operationg_cost + cvx_optval ;

    end
    
%     operationg_cost = sum(  Parameter.G(1:2,1).'* power(Pg_total(1:2,1:OPTIONS.N_t),2)  + Parameter.G(1:2,2).'*Pg_total(1:2,1:OPTIONS.N_t) + Parameter.G(1:2,3).'*ones(2,OPTIONS.N_t) ...
%                         + sum(Parameter.E(1,1)* power(Pb_total(1:2,1:OPTIONS.N_t),2) + Parameter.E(1,2)* power(Pb_total(1:2,1:OPTIONS.N_t),1) ,1) ,2 );
%     operationg_cost = sum(  Parameter.G(1:2,1).'* power(Pg_total(1:2,1:OPTIONS.N_t),2)  + Parameter.G(1:2,2).'*Pg_total(1:2,1:OPTIONS.N_t) + Parameter.G(1:2,3).'*ones(2,OPTIONS.N_t) ...
%                         + sum(Parameter.E(1,1)* power(Pb_total(1:2,1:OPTIONS.N_t),2)  ,1) ,2 );
%     operationg_cost = sum(  Parameter.G(1:2,1).'* power(Pg_total(1:2,1:OPTIONS.N_t),2)  + Parameter.G(1:2,2).'*Pg_total(1:2,1:OPTIONS.N_t) + Parameter.G(1:2,3).'*ones(2,OPTIONS.N_t) ,2 );

    y = size(find(delta(1:2,2:OPTIONS.N_t) - delta(1:2,1:OPTIONS.N_t-1)==1),2);
    operationg_cost = operationg_cost + y*Parameter.C_ss(1);
    
end