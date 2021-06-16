function P_grid=calculate_depth(N_res,u,v,par)
nC=20;
lambdas = par*ones(nC-3, nC-3);
parfor i=1:size(u,1)
    P{i} = calc_depth(i,u,v,nC,lambdas,N_res);
end
P_grid=[];
for i=1:size(u,1)
    P_grid = [P_grid;P{i}];
end
end

function Pg = calc_depth(i,u,v,nC,lambdas,N_res)
idx = find(u(i,:)~=0  & v(i,:)~=0 );
Pg = zeros(3,size(u,2));
    if ~isempty(idx)
        umin=min(u(i,:))-0.1;umax=max(u(i,:))+0.1;
        vmin=min(v(i,:))-0.1;vmax=max(v(i,:))+0.1;
        bbsd = bbs_create(umin, umax, nC, vmin, vmax, nC, 1);
        colocd = bbs_coloc(bbsd, u(i,idx), v(i,idx));
        bendingd = bbs_bending(bbsd, lambdas);
        [ctrlpts3Dn]=ShapeFromNormals(bbsd,colocd,bendingd,[u(i,idx);v(i,idx);ones(1,length(u(i,idx)))],N_res(3*(i-1)+1:3*(i-1)+3,idx));
        mu=bbs_eval(bbsd, ctrlpts3Dn,u(i,idx)',v(i,idx)',0,0);
        Pg(:,idx) = [u(i,idx);v(i,idx);ones(1,length(u(i,idx)))].*[mu;mu;mu];
    end
end

% for i=1:size(u,1)
%     idx = find(u(i,:)~=0  & v(i,:)~=0);
%     if ~isempty(idx)
%         umin=min(u(i,idx))-0.1;umax=max(u(i,idx))+0.1;
%         vmin=min(v(i,idx))-0.1;vmax=max(v(i,idx))+0.1;
%         bbsd = bbs_create(umin, umax, nC, vmin, vmax, nC, 1);
%         colocd = bbs_coloc(bbsd, u(i,idx), v(i,idx));
%         bendingd = bbs_bending(bbsd, lambdas);
%         [ctrlpts3Dn]=ShapeFromNormals(bbsd,colocd,bendingd,[u(i,idx);v(i,idx);ones(1,length(u(i,idx)))],N_res(3*(i-1)+1:3*(i-1)+3,idx));
%         mu=bbs_eval(bbsd, ctrlpts3Dn,u(i,idx)',v(i,idx)',0,0);
%         P_grid(3*(i-1)+1:3*(i-1)+3,idx) = [u(i,idx);v(i,idx);ones(1,length(u(i,idx)))].*[mu;mu;mu];
%     end
% end