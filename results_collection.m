for index_distance = 1:6
    cost_comparison(index_distance) = data_1th(index_distance).optimal_cost_related(2,end);
    RD(index_distance) = roundn(data_1th(index_distance).RD, -2);
    LS(index_distance) = roundn(sum(sum(data_1th(index_distance).LS,1), 2), -2);
    No_sw(index_distance) = roundn(sum(sum(data_1th(index_distance).No_sw,1), 2), -2);
end