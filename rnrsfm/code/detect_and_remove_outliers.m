function [out1,out2] = detect_and_remove_outliers(q,u,v,Pg,P,list,num,idd)
parfor i = 1: length(list)
    i
    [out_det{i},outlier{i}] = evaluate_outliers(q,u,v,Pg,P,list{i},num,idd);
end
out1 = [];
out2 = [];
for i = 1:length(list)
    out1 = [out1;out_det{i}];
    out2 = [out2;outlier{i}];
end
