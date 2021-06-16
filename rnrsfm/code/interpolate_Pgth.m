function P2 = interpolate_Pgth(P,u,v,q)
er = 1e-4;
nC = 20;
P2 = zeros(3*length(u(:,1)),length(q(1,:)));
for i = 1: size(u,1)
    
    idx = find(u(i,:)~=0 & v(i,:)~=0);
    idx2 = q(2*(i-1)+1,:)~=0;
    q1 = q(2*(i-1)+1:2*(i-1)+2,idx2);
    umin=min([u(i,idx),q1(1,:)])-0.001;umax=max([u(i,idx),q1(1,:)])+0.001;
    vmin=min([v(i,idx),q1(2,:)])-0.001;vmax=max([v(i,idx),q1(2,:)])+0.001;
    
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, u(i,idx), v(i,idx));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*P(3*(i-1)+1:3*(i-1)+3,idx)');
    ctrlpts = cpts';
    qw = bbs_eval(bbs, ctrlpts, q1(1,:)',q1(2,:)',0,0);
    P2(3*(i-1)+1:3*(i-1)+3,idx2) = qw;
end