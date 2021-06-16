function [in_det,inlier2d] = evaluate_outliers(q,u,v,Pg,P,list,num,idd)

inlier2d = ones(length(list),length(P(1,:)));
id = list;
q_id = [2*(id-1)+1;2*(id-1)+2];
q_id = q_id(:);
P_id = [3*(id-1)+1;3*(id-1)+2;3*(id-1)+3];
P_id = P_id(:);
qg= zeros(size(q,1),size(u,2));
qg(1:2:end,:)=u;
qg(2:2:end,:)=v;
[u,v,u1,v1,~] = create_grid(qg(q_id,:),ones(length(id),size(u,2)),40); 
q_g= zeros(2*size(u,1)+2,size(u,2));
q_g(1:2:end,:)=[u(1,:);u1];
q_g(2:2:end,:)=[v(1,:);v1];
P_g = interpolate_Pgth(Pg(P_id,:),qg(q_id(1:2:end),:),qg(q_id(2:2:end),:),q_g);

% compute nearest neighbors
nng = find_nng_grid(q(q_id,:),q_g,num);

%find inliers-outliers
in_det = detect_outliers(P(P_id,:),P_g,nng);
% remove from the 2d points 
for j = 1: length(id)
    l = idd{id(j)};
    if ~isempty(l)
        for k = 1:length(l)
            inlier2d(j,l(k))=0;
        end
    end
end

