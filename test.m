clear;

OPTIONS.Pg_Max = [8 5];
OPTIONS.Pg_Min = [1 1];

master_optval= zeros(4,4);

% %% total 
% for index_u = 1:4
%     if index_u == 1     
%         u1 = 0;
%         u2 = 1;
%     elseif index_u == 2     
%         u1 = 1;
%         u2 = 0;
%     elseif index_u == 3     
%         u1 = 1;
%         u2 = 1;
%     elseif index_u == 4    
%         u1 = 0;
%         u2 = 0;
%     end 
%     cvx_begin
%     % cvx_solver SeDuMi
%     cvx_begin quiet
%         variable Pg(2) nonnegative
%         variable P_c(2) nonnegative
%         variable P_d(2) nonnegative
%     %     dual variable y1
%     %     dual variable y2
%     %     variable u1 nonnegative
%     %     variable u2 nonnegative
%         minimize( 10*u1*power(Pg(1),2) + 5*u1*Pg(1) + u1*14 + 8*u2*power(Pg(2),2) + 4*u2*Pg(2) + u2*19 + 165*P_c(1) + 165*P_c(2) + 165*P_d(1) + 165*P_d(2)  )
%         subject to
%             % the range constraints of all the variables
%             Pg(1) <= u1 * OPTIONS.Pg_Max(1) + P_c(1)
%             Pg(2) <= u2 * OPTIONS.Pg_Max(2) + P_c(2)
%             Pg(1) >= u1 * OPTIONS.Pg_Min(1) - P_d(1)
%             Pg(2) >= u2 * OPTIONS.Pg_Min(2) - P_d(2)
% 
%     %         y1 : u1==1
%     %         y2 : u2==0
% 
%             Pg(1) + Pg(2) == 5
%     cvx_end
% 
%     disp(cvx_optval)
%     cvx_optval_total(index_u) = cvx_optval;
% end

% u1=0;
% u2=1;
u1_m = 0;
u2_m = 1;
for index_benders = 1:10
    cvx_begin
    % cvx_solver SeDuMi
    cvx_begin quiet
        variable Pg(2) nonnegative
        variable P_c(2) nonnegative
        variable P_d(2) nonnegative
        dual variable y1
        dual variable y2
        dual variable y3
        dual variable y4
        dual variable d1
        dual variable d2
        variable u1 nonnegative
        variable u2 nonnegative
        minimize( 10*power(Pg(1),2) + 5*Pg(1) + 8*power(Pg(2),2) + 4*Pg(2) + 165*P_c(1) + 165*P_c(2) + 165*P_d(1) + 165*P_d(2)  )
        subject to
            % the range constraints of all the variables
            y1 : Pg(1) <= u1 * OPTIONS.Pg_Max(1) + P_c(1)
            y2 : Pg(1) >= u1 * OPTIONS.Pg_Min(1) - P_d(1)
            y3 : Pg(2) <= u2 * OPTIONS.Pg_Max(2) + P_c(2)
            y4 : Pg(2) >= u2 * OPTIONS.Pg_Min(2) - P_d(2)

            d1 : u1==u1_m
            d2 : u2==u2_m

            Pg(1) + Pg(2) == 5
    cvx_end

    disp(cvx_optval)

%     sub_optval = cvx_optval;
    sub_optval(index_benders) = cvx_optval;
    u1_s(index_benders)=u1_m;
    u2_s(index_benders)=u2_m;
    upper_bound(index_benders) = sub_optval(index_benders) + u1_s(index_benders)*14 + u2_s(index_benders)*19;
    
    % dual variable and upperbound
%     lambda_u1(index_benders) = -y1*OPTIONS.Pg_Max(1) + y2*OPTIONS.Pg_Min(1);
%     lambda_u2(index_benders) = -y3*OPTIONS.Pg_Max(2) + y4*OPTIONS.Pg_Min(2);
    lambda_u1(index_benders) = d1;
    lambda_u2(index_benders) = d2;

    cvx_begin
    cvx_solver Mosek
    % cvx_begin quiet
        variable u1_m binary
        variable u2_m binary
        variable benders_cut
        minimize( u1_m*14 + u2_m*19 + benders_cut )

        subject to
    %         % startup detect
    %         startup(1,1:OPTIONS.N_t-1) >= (delta_master(1,2:OPTIONS.N_t) - delta_master(1,1:OPTIONS.N_t-1))
    %         startup(2,1:OPTIONS.N_t-1) >= (delta_master(2,2:OPTIONS.N_t) - delta_master(2,1:OPTIONS.N_t-1))

             % benders cuts
             for index = 1:index_benders
                benders_cut >= sub_optval(index) + lambda_u1(index)*(u1_m - u1_s(index)) + lambda_u2(index)*(u2_m - u2_s(index))
             end
    cvx_end

    disp('optimal');
    disp(cvx_optval);

    benders_cut_st(index_benders) = benders_cut;
    lower_bound(index_benders) = cvx_optval;
    error(index_benders) = (upper_bound(index_benders) - lower_bound(index_benders)) ;
    if (upper_bound(index_benders) - lower_bound(index_benders)) < 1e-1
        break;
    end
%     for index_temp = 1:index_benders
% %         for index_cal = 1:index_benders
% %         for index = 1:4
% %             if index == 1
%             u1_m = 1; u2_m = 1;
%             master(index_benders).optval(index_temp,index) = u1_m*14 + u2_m*19 + sub_optval(index_temp) + lambda_u1(index_temp)*(u1_m - u1_s(index_temp)) + lambda_u2(index_temp)*(u2_m - u2_s(index_temp));
% %             elseif index == 2
%             u1_m = 0; u2_m = 1;
%             master(index_benders).optval(index_temp,index) = u1_m*14 + u2_m*19 + sub_optval(index_temp) + lambda_u1(index_temp)*(u1_m - u1_s(index_temp)) + lambda_u2(index_temp)*(u2_m - u2_s(index_temp));
% %             elseif index == 3
%             u1_m = 1; u2_m = 0;
%             master(index_benders).optval(index_temp,index) = u1_m*14 + u2_m*19 + sub_optval(index_temp) + lambda_u1(index_temp)*(u1_m - u1_s(index_temp)) + lambda_u2(index_temp)*(u2_m - u2_s(index_temp));
% %             else
%             u1_m = 0; u2_m = 0;
% %             end
%             master(index_benders).optval(index_temp,index) = u1_m*14 + u2_m*19 + sub_optval(index_temp) + lambda_u1(index_temp)*(u1_m - u1_s(index_temp)) + lambda_u2(index_temp)*(u2_m - u2_s(index_temp));
% %         end
%     end
    
end    

for index_u = 1:4
    if index_u == 1
        u1_m = 0;
        u2_m = 1;
    elseif index_u == 2     
        u1_m = 1;
        u2_m = 0;
    elseif index_u == 3     
        u1_m = 1;
        u2_m = 1;
    elseif index_u == 4    
        u1_m = 0;
        u2_m = 0;
    end
        
    cvx_begin
    % cvx_solver Mosek
    % cvx_begin quiet
    %     variable u1_m binary
    %     variable u2_m binary
        variable benders_cut
        minimize( u1_m*14 + u2_m*19 + benders_cut )

        subject to
    %         % startup detect
    %         startup(1,1:OPTIONS.N_t-1) >= (delta_master(1,2:OPTIONS.N_t) - delta_master(1,1:OPTIONS.N_t-1))
    %         startup(2,1:OPTIONS.N_t-1) >= (delta_master(2,2:OPTIONS.N_t) - delta_master(2,1:OPTIONS.N_t-1))

             % benders cuts
            benders_cut >= sub_optval + lambda_u1*(u1_m - u1) + lambda_u2*(u2_m - u2)

    cvx_end

    disp('optimal');
    disp(cvx_optval);
    
    cvx_optval_master(index_u) = cvx_optval;

end


% cvx_begin
% % cvx_solver SeDuMi
% cvx_begin quiet
% %     variable x
%     variable x
%     variable y
%     dual variable d1
%     minimize(  -x -y   )
%     subject to
%         % the range constraints of all the variables
%         0.5*exp(2*y) - x <= 0.25
%         y >= 0
%         y <= 0.5
%         d1 : x==1
% cvx_end
