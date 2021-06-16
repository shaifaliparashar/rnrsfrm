function [J21a,J21b,J21c,J21d,H21uva,H21uvb] = create_warps(I1u,I1v,I2u,I2v)

er = 1e-5;
t= 1e-3;
nC = 20;
J21a = zeros(size(I2u)); J21b = J21a; J21c=J21a; J21d = J21a; H21uva = J21a; H21uvb = J21a;
% calculate Eta_21 derivatives using schwarzian warps
parfor i = 1:size(I2u,1)
    q1 = [I1u(1,:);I1v(1,:)];  q2 = [I2u(i,:);I2v(i,:)];
    umin = min([q2(1,:),q1(1,:)])-t; umax = max([q2(1,:),q1(1,:)])+t; 
    vmin = min([q2(2,:),q1(2,:)])-t; vmax = max([q2(2,:),q1(2,:)])+t;
    
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
    coloc = bbs_coloc(bbs, q2(1,:), q2(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*q1');
    ctrlpts = cpts';
      dqu = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',1,0); 
    dqv = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',0,1);
    [xv,yv]=meshgrid(linspace(umin,umax,30),linspace(vmin,vmax,30));
    ctrlpts = optimPanalSchwarz(bbs,ctrlpts,q2',q1',[xv(:),yv(:)],1e-3);
    ctrlpts=ctrlpts';
       
    qw2 = bbs_eval(bbs,ctrlpts,q2(1,:)',q2(2,:)',0,0);
    error=sqrt(mean((qw2(1,:)-q1(1,:)).^2+(qw2(2,:)-q1(2,:)).^2));
    dqu = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',1,0); 
    dqv = bbs_eval(bbs, ctrlpts, q2(1,:)',q2(2,:)',0,1);
    dquv = bbs_eval(bbs,ctrlpts,q2(1,:)',q2(2,:)',1,1);
    J21a(i,:) = dqu(1,:); J21b(i,:) = dqu(2,:);
    J21c(i,:) = dqv(1,:); J21d(i,:) = dqv(2,:);
    H21uva(i,:) = dquv(1,:); H21uvb(i,:) = dquv(2,:);
end
