function [I1u,I1v,I2u,I2v,visb2] = create_grid(q_n,visb,p)

n = size(visb,1);
er = 1e-4;
t= 1e-3;
nC = 20;
idx = visb(1,:)==1;
q1 = q_n(1:2,idx);  
umin = min(q1(1,:))-t; umax = max(q1(1,:))+t;
vmin = min(q1(2,:))-t; vmax = max(q1(2,:))+t;
bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);

[xv,yv]=meshgrid(linspace(bbs.umin,bbs.umax,p),linspace(bbs.vmin,bbs.vmax,p));
I1u = repmat(xv(:)',n-1,1);
I1v = repmat(yv(:)',n-1,1);
visb2 = ones(n,length(I1u(1,:)));
for i = 2:n
    idx = find(visb(1,:)==1 & visb(i,:)==1);
    if ~isempty(idx)
        q1 = q_n(1:2,idx);
        q2 = q_n(2*(i-1)+1:2*(i-1)+2,idx);
        umin = min([q1(1,:),I1u(1,:)])-t; umax = max([q1(1,:),I1u(1,:)])+t;
        vmin = min([q1(2,:),I1v(1,:)])-t; vmax = max([q1(2,:),I1v(1,:)])+t;
        bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
        lambdas = er*ones(nC-3, nC-3);
        coloc = bbs_coloc(bbs, q_n(1,idx), q_n(2,idx));
        bending = bbs_bending(bbs, lambdas);
        
        cpts = (coloc'*coloc + bending) \ (coloc'*q2');
        ctrlpts = cpts';
        qw =  bbs_eval(bbs,ctrlpts,q_n(1,idx)',q_n(2,idx)',0,0);
        error=sqrt(mean((qw(1,:)-q_n(2*(i-1)+1,idx)).^2+(qw(2,:)-q_n(2*(i-1)+2,idx)).^2));
        q =  bbs_eval(bbs,ctrlpts,I1u(1,:)',I1v(1,:)',0,0);
        I2u(i-1,:) = q(1,:);
        I2v(i-1,:) = q(2,:);
    end
end

