function P2=median_shape_all(N_a,u_all,v_all,list,ref,par,pts,th)
N_med = zeros(3*size(u_all,1),pts);

for k = 1: length(list)
    Ne = N_a(:,(k-1)*length(ref{k})*pts+1:(k-1)*length(ref{k})*pts+length(ref{k})*pts);
    s = sum(abs(Ne'));
    idx = find(s==0);
    if ~isempty(idx)
        Ne(idx,:)=[];
    end
    % normal similarity
    for i = 1: pts
        for j = 1:length(list{k}) % reconstructing jth view
            N_in = Ne(3*j-2:3*j,pts*(0:length(ref{k})-1)+i);
            idd = find(sum(abs(N_in))==0);
            if ~isempty(idd)
                N_in(:,idd) = [];
            end
            stop = 0; res_old = 1e5; nc = [0;0;0];
            if ~isempty(N_in)
                while ~stop
                    sim = [];
                    for ii = 1:  size(N_in,2)
                        for jj = ii+1: size(N_in,2)
                            sim(ii,jj) = min(acosd(dot(N_in(:,ii),N_in(:,jj))),acosd(dot(-N_in(:,ii),N_in(:,jj))));
                            sim(jj,ii) = sim(ii,jj);
                        end
                    end
                    sim(sim==0)=[]; sim=reshape(sim,size(N_in,2),size(N_in,2)-1);
                    med_sim = median(sim);[val, idx] = min(med_sim);
                    if val < res_old
                        if size(N_in,2)<=length(ref{k})-2
                            N_f = [median(nonzeros(N_in(1,:)));median(nonzeros(N_in(2,:)));median(nonzeros(N_in(3,:)))];
                            N_f = pick_closest(N_f,N_in); 
                            stop =1;
                        else
                            res_old = val;
                            %nc = N_in(:,idx);
                            nc = [median(nonzeros(N_in(1,:)));median(nonzeros(N_in(2,:)));median(nonzeros(N_in(3,:)))];
                            nc = pick_closest(nc,N_in);
                            [~, idx] = max(med_sim);N_in(:,idx) = [];
                            if val < th
                                N_f = nc; stop = 1;
                            end   
                        end
                    else
                        N_f = nc; stop = 1;
                    end  
                end
            end
            N_med(3*(list{k}(j)-1)+1:3*(list{k}(j)-1)+3,i) = N_f;
        end
    end
end

% integrate normals
P2=calculate_depth(N_med,u_all,v_all,par);


end

function Nc = pick_closest(N,Ni)
t=acosd(dot(repmat(N,1,size(Ni,2)),Ni));
[val, idx]= min(t);
Nc = Ni(:,idx);
end



