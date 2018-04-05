dbstop if error
operation_mode = 7;
index_time = 12;
Rest_ppr_avg = [3.85611815044505];
delta = [1,1,1,1,1,1,1,1,1,1,1,1;1,1,1,1,1,1,1,1,1,1,1,1];
E = [2,1,0,0,0,0,0,0,0,0,0;2,1.991,1.955,1.904,1.864,1.855,1.875,1.901,1.941,1.970,1.989];
Pb = [0,1,1,0,0,0,0,0,0,0,0;0,0.009,0.037,0.051,0.04,0.009,-0.019,-0.026,-0.04,-0.029,-0.019];
Pg = [2,4,5.24,6.34,6.41,6.41,6.35,6.28,6.14,5.98,5.80; 0.550,0.660,0.880,0.990,0.900,0.660,0.440,0.380,0.270,0.360,0.440];

Ppr = [4.36,4.186,3.861,3.744,3.959,4.399,4.753,4.791,4.883,4.595,4.31];
load_shedding = [0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0];
Rest_Distance = [12.0570913247533];
Rest_ppr_avg = [4.23537198095950];
redundent_sw = [1,1,1,1,1,1,1,1,1,1,1,1;0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0;1,1,1,1,1,1,1,1,1,1,1,1];


% OPTIONS.R_G = OPTIONS.Pg_Max*0.75;

% switch index_time
%     case 12
%         varphi_sub = 1;
%     otherwise
%         varphi_sub = OPTIONS.varphi;
% end
% 
% Delta_PL_island1 = OPTIONS.island1_load(index_time) - OPTIONS.island1_load_average;
% Delta_PL_island2 = OPTIONS.island2_load(index_time) - OPTIONS.island2_load_average;
% Distance_slot_obj = ((Rest_ppr_avg - (1 - varphi_sub)*Delta_PL_island1)/2.2e-3)^(1/3);
% 
% % ESM output power and the capacity constraints
% if index_time == 1
%     rest_pmax_ESM1 = OPTIONS.E_Max(1) - OPTIONS.E_Min(1) ;
%     rest_pmax_ESM2 = OPTIONS.E_Max(2) - OPTIONS.E_Min(2) ;
%     rest_pmin_ESM1 = OPTIONS.E_Max(1) - OPTIONS.E_Max(1) ;
%     rest_pmin_ESM2 = OPTIONS.E_Max(2) - OPTIONS.E_Max(2) ;
% else
%     rest_pmax_ESM1 = E(1,index_time-1) - OPTIONS.E_Min(1) ;
%     rest_pmax_ESM2 = E(2,index_time-1) - OPTIONS.E_Min(2) ;
%     rest_pmin_ESM1 = E(1,index_time-1) - OPTIONS.E_Max(1) ;
%     rest_pmin_ESM2 = E(2,index_time-1) - OPTIONS.E_Max(2) ;
% end
% upper_bound_ESM1P = roundn( min( varphi_sub * Delta_PL_island1, rest_pmax_ESM1 ), -4);
% upper_bound_ESM1N = roundn( max( varphi_sub * Delta_PL_island1, rest_pmin_ESM1 ), -4);
% upper_bound_ESM2P = roundn( min( OPTIONS.varphi * Delta_PL_island2, rest_pmax_ESM2 ), -4);
% upper_bound_ESM2N = roundn( max( OPTIONS.varphi * Delta_PL_island2, rest_pmin_ESM2 ), -4);

        switch index_time
            case OPTIONS.N_t
                varphi_sub = 1;
            otherwise
                varphi_sub = OPTIONS.varphi;
        end
        
        if operation_mode <= 3
            Delta_PL =  OPTIONS.P_L_TIME_off(index_time) - OPTIONS.P_L_TIME_off_avg ;
            Distance_slot_obj = ((Rest_ppr_avg - (1 - varphi_sub)*Delta_PL)/2.2e-3)^(1/3);
            
            % ESM output power and the capacity constraints
            if index_time == 1
                rest_pmax_ESM = OPTIONS.E_Max(1) + OPTIONS.E_Max(2) - OPTIONS.E_Min(1) - OPTIONS.E_Min(2);
                rest_pmin_ESM = OPTIONS.E_Max(1) + OPTIONS.E_Max(2)  -OPTIONS.E_Max(1) - OPTIONS.E_Max(2);
            else
                rest_pmax_ESM = E(1,index_time-1) + E(2,index_time-1) - OPTIONS.E_Min(1) - OPTIONS.E_Min(2) ;
                rest_pmin_ESM = E(1,index_time-1) + E(2,index_time-1) - OPTIONS.E_Max(1) - OPTIONS.E_Max(2) ;
            end
            upper_bound_ESM_P = roundn( min( varphi_sub * Delta_PL, rest_pmax_ESM ), -4);
            upper_bound_ESM_N = roundn( max( varphi_sub * Delta_PL, rest_pmin_ESM ), -4);
        elseif operation_mode <= 7
            Delta_PL_island1 = OPTIONS.island1_load(index_time) - OPTIONS.island1_load_average;
            Delta_PL_island2 = OPTIONS.island2_load(index_time) - OPTIONS.island2_load_average;
            Distance_slot_obj = ((Rest_ppr_avg - (1 - varphi_sub)*Delta_PL_island1)/2.2e-3)^(1/3);
            
            % ESM output power and the capacity constraints
            if index_time == 1
                rest_pmax_ESM1 = OPTIONS.E_Max(1) - OPTIONS.E_Min(1) ;
                rest_pmax_ESM2 = OPTIONS.E_Max(2) - OPTIONS.E_Min(2) ;
                rest_pmin_ESM1 = OPTIONS.E_Max(1) - OPTIONS.E_Max(1) ;
                rest_pmin_ESM2 = OPTIONS.E_Max(2) - OPTIONS.E_Max(2) ;
            else
                rest_pmax_ESM1 = E(1,index_time-1) - OPTIONS.E_Min(1) ;
                rest_pmax_ESM2 = E(2,index_time-1) - OPTIONS.E_Min(2) ;
                rest_pmin_ESM1 = E(1,index_time-1) - OPTIONS.E_Max(1) ;
                rest_pmin_ESM2 = E(2,index_time-1) - OPTIONS.E_Max(2) ;
            end
            if sum(delta(1,index_time))>=1
                upper_bound_ESM1P = roundn( min( varphi_sub * Delta_PL_island1, rest_pmax_ESM1 ), -4);
                upper_bound_ESM1N = roundn( max( varphi_sub * Delta_PL_island1, rest_pmin_ESM1 ), -4);
%                 upper_bound_ESM2P = roundn( min( varphi_sub * Delta_PL_island2, rest_pmax_ESM2 ), -4);
%                 upper_bound_ESM2N = roundn( max( varphi_sub * Delta_PL_island2, rest_pmin_ESM2 ), -4);
            else
                upper_bound_ESM1P = roundn( rest_pmax_ESM1, -4);
                upper_bound_ESM1N = roundn( rest_pmin_ESM1, -4);
%                 upper_bound_ESM2P = roundn( min( varphi_sub * Delta_PL_island2, rest_pmax_ESM2 ), -4);
%                 upper_bound_ESM2N = roundn( max( varphi_sub * Delta_PL_island2, rest_pmin_ESM2 ), -4);
            end
            if sum(delta(2,index_time))>=1
%                 upper_bound_ESM1P = roundn( min( varphi_sub * Delta_PL_island1, rest_pmax_ESM1 ), -4);
%                 upper_bound_ESM1N = roundn( max( varphi_sub * Delta_PL_island1, rest_pmin_ESM1 ), -4);
                upper_bound_ESM2P = roundn( min( varphi_sub * Delta_PL_island2, rest_pmax_ESM2 ), -4);
                upper_bound_ESM2N = roundn( max( varphi_sub * Delta_PL_island2, rest_pmin_ESM2 ), -4);
            else
%                 upper_bound_ESM1P = roundn( min( varphi_sub * Delta_PL_island1, rest_pmax_ESM1 ), -4);
%                 upper_bound_ESM1N = roundn( max( varphi_sub * Delta_PL_island1, rest_pmin_ESM1 ), -4);
                upper_bound_ESM2P = roundn( rest_pmax_ESM2, -4);
                upper_bound_ESM2N = roundn( rest_pmin_ESM2, -4);
            end
        end

cvx_begin quiet
    variable Pg_on(OPTIONS.N_g) nonnegative
    variable Pb_on(2,1)
    variable E_on(2,1) nonnegative
    variable Ppr_on(1) nonnegative
    variable delta_s(OPTIONS.N_g)
    variable redundent_sw_s(4,1)
    variable load_shedding_on(2, 1) nonnegative
    variable reduced_distance_on nonnegative
    variable Pc(OPTIONS.N_g) nonnegative
    variable Pd(OPTIONS.N_g) nonnegative
    dual variable temp_dual_delta_g1
    dual variable temp_dual_delta_g2
    dual variable temp_dual_Sp
    dual variable temp_dual_Ss
    minimize( sum(OPTIONS.G(1,1) * power(Pg_on(1,1).',2) ,1) ...
            + sum(OPTIONS.G(2,1) * power(Pg_on(2,1).',2) ,1) ...
            + sum(OPTIONS.G(1,2) * power(Pg_on(1,1).',1) ,1) ...
            + sum(OPTIONS.G(2,2) * power(Pg_on(2,1).',1) ,1) ...
            + OPTIONS.Xi_E * sum(OPTIONS.E(1:2,1).'* power(Pb_on(1:OPTIONS.N_e,1),2) ,2) ...
            + OPTIONS.Xi_E * sum([1 1] * ones(OPTIONS.N_e,1) ,2)...
            + sum(OPTIONS.Penalty_L * sum(load_shedding_on(1:2,1),2) ,1) ...
            + 1.02 * OPTIONS.Penalty_D * reduced_distance_on ...
            + sum(10 * OPTIONS.Penalty_D * Pc(1).' ,1) ...
            + sum(10 * OPTIONS.Penalty_D * Pc(2).' ,1) ...
            + sum(10 * OPTIONS.Penalty_D * Pd(1).' ,1) ...
            + sum(10 * OPTIONS.Penalty_D * Pd(2).' ,1))

    subject to
        % the range constraints of all the variables
        Pg_on(1,1) <= delta_s(1,1) * OPTIONS.Pg_Max(1) + Pc(1)
        Pg_on(1,1) >= delta_s(1,1) * OPTIONS.Pg_Min(1) - Pd(1)
        Pg_on(2,1) <= delta_s(2,1) * OPTIONS.Pg_Max(2) + Pc(2)
        Pg_on(2,1) >= delta_s(2,1) * OPTIONS.Pg_Min(2) - Pd(2)

        temp_dual_delta_g1 : delta_s(1,1) == delta(1,index_time)
        temp_dual_delta_g2 : delta_s(2,1) == delta(2,index_time)

        % ramping rate power of generators
        switch index_time
            case 1
                Pg_on(1,1) <= OPTIONS.R_G(1) 
                Pg_on(2,1) <= OPTIONS.R_G(2) 
            case OPTIONS.N_t
                Pg_on(1,1) <= OPTIONS.R_G(1) 
                Pg_on(2,1) <= OPTIONS.R_G(2)
            otherwise
                Pg_on(1,1) - Pg(1,index_time-1) <= OPTIONS.R_G(1) 
                Pg_on(1,1) - Pg(1,index_time-1) >= -OPTIONS.R_G(1) 
                Pg_on(2,1) - Pg(2,index_time-1) <= OPTIONS.R_G(2) 
                Pg_on(2,1) - Pg(2,index_time-1) >= -OPTIONS.R_G(2) 
        end
    %                 if index_time == 1
    %                     Pg_on(1,1) <= OPTIONS.R_G 
    %                     Pg_on(2,1) <= OPTIONS.R_G 
    %                 else
    %                     Pg_on(1,1) - Pg(1,index_time-1) <= OPTIONS.R_G 
    %                     Pg_on(1,1) - Pg(1,index_time-1) >= -OPTIONS.R_G 
    %                     Pg_on(2,1) - Pg(2,index_time-1) <= OPTIONS.R_G 
    %                     Pg_on(2,1) - Pg(2,index_time-1) >= -OPTIONS.R_G 
    %                 end

        % propulsion power limitation
        Ppr_on(1,1) <= OPTIONS.Ppr_Max

        % ESM limitation constraints
        if operation_mode <= 3
            if Delta_PL >= 0
                % the scaler factor is related with load power in each
                % island part and the adjusting factor.
                Pb_on(1,1) + Pb_on(2,1) <= 2 * OPTIONS.Pb_Max
                Pb_on(1,1) + Pb_on(2,1) <= upper_bound_ESM_P
                Pb_on(1,1) + Pb_on(2,1) >= 0
            elseif Delta_PL < 0
                Pb_on(1,1) + Pb_on(2,1) <= 0
                Pb_on(1,1) + Pb_on(2,1) <= upper_bound_ESM_N
                Pb_on(1,1) + Pb_on(2,1) >= 2 * OPTIONS.Pb_Min
            end
        elseif operation_mode <= 7
            if Delta_PL_island1 >=0
                % the scaler factor is related with load power in each
                % island part and the adjusting factor.
                Pb_on(1,1) <= OPTIONS.Pb_Max
%                 Pb_on(1,1) <= upper_bound_ESM1P
                Pb_on(1,1) >= 0
            elseif Delta_PL_island1 <0
                % the scaler factor is related with load power in each
                % island part and the adjusting factor.
                Pb_on(1,1) <= 0
%                 Pb_on(1,1) <= upper_bound_ESM1N
                Pb_on(1,1) >= OPTIONS.Pb_Min 
            end

            if Delta_PL_island2 >=0
                % the scaler factor is related with load power in each
                % island part and the adjusting factor.
                Pb_on(2,1) <= OPTIONS.Pb_Max
                Pb_on(2,1) <= upper_bound_ESM2P
                Pb_on(2,1) >= 0
            elseif Delta_PL_island2 <0
                % the scaler factor is related with load power in each
                % island part and the adjusting factor.
                Pb_on(2,1) <= 0
                Pb_on(2,1) <= upper_bound_ESM2N
                Pb_on(2,1) >= OPTIONS.Pb_Min
            end
        end


        % charging and discharging
        E_on(1,1) <=  OPTIONS.E_Max(1)
        E_on(2,1) <=  OPTIONS.E_Max(2)
        E_on(1,1) >= 0
        E_on(2,1) >= 0

        % ESM output power and the capacity constraints
        if index_time == 1
            OPTIONS.E_Max(1) - Pb_on(1,1) == E_on(1,1)
            OPTIONS.E_Max(2) - Pb_on(2,1) == E_on(2,1)
        else
            E(1,index_time-1) - Pb_on(1,1) == E_on(1,1)
            E(2,index_time-1) - Pb_on(2,1) == E_on(2,1)
        end

        % system power balance & load shedding
        if operation_mode <= 3
             OPTIONS.P_L_TIME_off(1,index_time) - sum(load_shedding_on(1,1)) + Ppr_on(1,1) == sum( Pg_on(1:OPTIONS.N_g,1) ) + sum(Pb_on(1:OPTIONS.N_e,1))
             load_shedding_on(2,1) == 0

            % load shedding range
            load_shedding_on(1) == zeros(1,1)
            load_shedding_on(2) == zeros(1,1)

        elseif operation_mode <= 7
            for t_index = index_time
                redundent_sw_s(1:2,1).'*OPTIONS.Coupled_load(:,index_time) + OPTIONS.island1_load(1, t_index) - load_shedding_on(1,1) + Ppr_on(1,1) == (Pg_on(1,1)) + Pb_on(1,1) 
                redundent_sw_s(3:4,1).'*OPTIONS.Coupled_load(:,index_time) + OPTIONS.island2_load(1, t_index) - load_shedding_on(2,1) == (Pg_on(2,1)) + Pb_on(2,1)
            end

            % load shedding range
            load_shedding_on(1) <= OPTIONS.island1_non_load(1,index_time)
            load_shedding_on(2) <= OPTIONS.island2_non_load(1,index_time)
        end

        temp_dual_Sp : redundent_sw_s(1:2,1) == redundent_sw(1:2,index_time)
        temp_dual_Ss : redundent_sw_s(3:4,1) == redundent_sw(3:4,index_time)

        % voyage planning
        switch operation_mode
            case 0
                Ppr_on(1,1) == OPTIONS.P_pr_avg;
                Pb_on(1:2,1) == 0;
            case 1
                (Ppr_on(1,1)./2.2e-3).^(1/3) >= Distance_slot_obj - reduced_distance_on;
                Pb_on(1:2,1) == 0;
            case 2
                Ppr_on(1,1) == OPTIONS.P_pr_avg;
            case 3
                (Ppr_on(1,1)./2.2e-3).^(1/3) >= Distance_slot_obj - reduced_distance_on; 
            case 4
                Ppr_on(1,1) == OPTIONS.P_pr_avg;
                Pb(1:2,1) == 0;
            case 5
                sum((Ppr_on(1,1)./2.2e-3).^(1/3)) >= Distance_slot_obj - reduced_distance_on;
                Pb_on(1:2,1) == 0;
            case 6
                Ppr_on(1,1) == OPTIONS.P_pr_avg;
            case 7
                sum((Ppr_on(1,1)./2.2e-3).^(1/3)) >= Distance_slot_obj - reduced_distance_on; 
        end

cvx_end

if strcmp(cvx_status, 'Infeasible')
    disp('stop');
end