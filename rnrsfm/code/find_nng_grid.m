function [nng] = find_nng_grid(q,qg,num)


for i = 1: length(q(1,:))
    for k = 1: length(qg(1,:))
        if q(1,i) ~= 0
            dist_q(i,k) = sum(abs(q(1:2,i)-qg(1:2,k)));
        else
            dist_q(i,k) = 1000;
        end
        
    end
    [val,id2]=sort(dist_q(i,:));
    idx = id2(1:num);
    nng{i} = idx;

    
end
