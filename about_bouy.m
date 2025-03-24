%% step1
b_Date = b_Date(~isnan(b_Hs));
b_Hs = b_Hs(~isnan(b_Hs));
b_Wind = b_Wind(~isnan(b_Hs));
%% step2
b_Hs(b_Hs == 0) = NaN;
b_Hs = fillmissing(b_Hs, 'linear');