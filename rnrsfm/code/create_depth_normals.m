function [P2,N] = create_depth_normals(res,a,b,c,d,t1,t2,visb,I1u,I2u,I1v,I2v,q_n)
num=size(visb,1);
% recover first order derivatives on rest of the surfaces
k1_all = [res(:,1)';a.*repmat(res(:,1)',num-1,1) + b.*repmat(res(:,2)',num-1,1) + t1];
k2_all = [res(:,2)';c.*repmat(res(:,1)',num-1,1) + d.*repmat(res(:,2)',num-1,1) + t2];

idx = find(visb(1,:)==0);
for i = 1: length(idx)
    id = find(visb(1:end,idx(i))>0);
    I2u(id(1)-1,idx(i)) = I1u(1,idx(i)); I2v(id(1)-1,idx(i)) = I1v(1,idx(i)); I1u(:,idx(i)) = 0; I1v(:,idx(i)) = 0;
    k1_all(id(1),idx(i)) = k1_all(1,idx(i)); k2_all(id(1),idx(i)) = k2_all(1,idx(i)); k1_all(1,idx(i)) = 0; k2_all(1,idx(i)) = 0;
end

u_all = [I1u(1,:);I2u]; v_all = [I1v(1,:);I2v];

% find normals on all surfaces N= [N1;N2;N3]
N1 = k1_all; N2 = k2_all; N3 = 1-u_all.*k1_all-v_all.*k2_all;
n = sqrt(N1.^2+N2.^2+N3.^2);
N1 = N1./n ; N2 = N2./n; N3 = N3./n;

N = [N1(:),N2(:),N3(:)]';
N_res = reshape(N(:),3*num,size(u_all,2));

% find indices with no solution
idx = find(res(:,1)==0);
N_res(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = [];

% Integrate normals to find depth
if ~isempty(N_res)
P_grid=calculate_depth(N_res,u_all,v_all,1e0);

P2 = interpolate_Pgth(P_grid,u_all,v_all,q_n);
N = obtain_normals(P2,q_n);
else
    P2 = []; N=[];
end