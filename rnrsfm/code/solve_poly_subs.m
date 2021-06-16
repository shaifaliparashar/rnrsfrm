function res = solve_poly_subs(a,b,c,d,t1,t2,e2,e1,u2,u1,v2,v1)

eq121 = e1.*v2.*-2.0+c.*e2.*u1.*2.0+e1.*e2.*t2.*2.0+d.*e2.*v1.*2.0;
eq112 = e1.*u2.*2.0-a.*e2.*u1.*2.0-b.*e2.*v1.*2.0-e1.*e2.*t1.*2.0;
eq120 = e1-c.^2.*e2-d.^2.*e2-e1.*t2.*v2.*2.0+e1.*e2.*t2.^2;
eq111 = a.*u1.*v2.*4.0-c.*u1.*u2.*4.0+b.*v1.*v2.*4.0-d.*u2.*v1.*4.0-a.*e2.*t2.*u1.*4.0+c.*e2.*t1.*u1.*4.0-b.*e2.*t2.*v1.*4.0+d.*e2.*t1.*v1.*4.0;
eq102 = -e1+a.^2.*e2+b.^2.*e2+e1.*t1.*u2.*2.0-e1.*e2.*t1.^2;
eq110 = a.*u1.*-2.0-b.*v1.*2.0+c.^2.*u2.*2.0+d.^2.*u2.*2.0-c.^2.*e2.*t1.*2.0-d.^2.*e2.*t1.*2.0+a.*t2.*u1.*v2.*4.0+b.*t2.*v1.*v2.*4.0-a.*e2.*t2.^2.*u1.*2.0-b.*e2.*t2.^2.*v1.*2.0;
eq101 = c.*u1.*2.0+d.*v1.*2.0-a.^2.*v2.*2.0-b.^2.*v2.*2.0+a.^2.*e2.*t2.*2.0+b.^2.*e2.*t2.*2.0-c.*t1.*u1.*u2.*4.0-d.*t1.*u2.*v1.*4.0+c.*e2.*t1.^2.*u1.*2.0+d.*e2.*t1.^2.*v1.*2.0;
eq100 = a.^2+b.^2-c.^2-d.^2-a.^2.*t2.*v2.*2.0+c.^2.*t1.*u2.*2.0-b.^2.*t2.*v2.*2.0+d.^2.*t1.*u2.*2.0+a.^2.*e2.*t2.^2+b.^2.*e2.*t2.^2-c.^2.*e2.*t1.^2-d.^2.*e2.*t1.^2;

eqa(:,:,1)=eq121; eqa(:,:,2)=eq112; eqa(:,:,3)=eq120; eqa(:,:,4)=eq111; eqa(:,:,5)=eq102; eqa(:,:,6)=eq110; eqa(:,:,7)=eq101; eqa(:,:,8)=eq100;

eq212 = -e1.*v2+c.*e2.*u1+e1.*e2.*t2+d.*e2.*v1;
eq203 = e1.*u2-a.*e2.*u1-b.*e2.*v1-e1.*e2.*t1;
eq211 = e1-c.^2.*e2-d.^2.*e2-e1.*t2.*v2.*2.0+e1.*e2.*t2.^2;
eq202 = a.*c.*e2+b.*d.*e2+a.*u1.*v2.*2.0-c.*u1.*u2.*2.0+b.*v1.*v2.*2.0+e1.*t2.*u2-d.*u2.*v1.*2.0+e1.*t1.*v2-a.*e2.*t2.*u1.*2.0+c.*e2.*t1.*u1.*2.0-b.*e2.*t2.*v1.*2.0-e1.*e2.*t1.*t2+d.*e2.*t1.*v1.*2.0;
eq210 = -c.*u1-d.*v1+c.^2.*v2+d.^2.*v2-c.^2.*e2.*t2-d.^2.*e2.*t2+c.*t2.*u1.*v2.*2.0+d.*t2.*v1.*v2.*2.0-c.*e2.*t2.^2.*u1-d.*e2.*t2.^2.*v1;
eq201 = -a.*u1-b.*v1+c.^2.*u2+d.^2.*u2-a.*c.*v2.*2.0-b.*d.*v2.*2.0-c.^2.*e2.*t1-d.^2.*e2.*t1+a.*t2.*u1.*v2.*2.0-c.*t2.*u1.*u2.*2.0-c.*t1.*u1.*v2.*2.0+b.*t2.*v1.*v2.*2.0-d.*t2.*u2.*v1.*2.0-d.*t1.*v1.*v2.*2.0-a.*e2.*t2.^2.*u1-b.*e2.*t2.^2.*v1+a.*c.*e2.*t2.*2.0+b.*d.*e2.*t2.*2.0+c.*e2.*t1.*t2.*u1.*2.0+d.*e2.*t1.*t2.*v1.*2.0;
eq200 = a.*c+b.*d+c.^2.*t2.*u2+c.^2.*t1.*v2+d.^2.*t2.*u2+d.^2.*t1.*v2-a.*c.*t2.*v2.*2.0-b.*d.*t2.*v2.*2.0+a.*c.*e2.*t2.^2+b.*d.*e2.*t2.^2-c.^2.*e2.*t1.*t2-d.^2.*e2.*t1.*t2;

eqb(:,:,1)=eq212; eqb(:,:,2)=eq203; eqb(:,:,3)=eq211; eqb(:,:,4)=eq202; eqb(:,:,5)=eq210; eqb(:,:,6)=eq201; eqb(:,:,7)=eq200;
DD = eq211.^2-4.*eq212.*eq210;

%eq7 = eq121.*eq203.^2-eq112.*eq203.*eq212;
eq6 = eq102.*eq212.^2+eq120.*eq203.^2-eq111.*eq203.*eq212-eq112.*eq202.*eq212-eq112.*eq203.*eq211+eq121.*eq202.*eq203.*2.0;
eq5 = eq101.*eq212.^2+eq121.*eq202.^2+eq102.*eq211.*eq212.*2.0-eq110.*eq203.*eq212-eq111.*eq202.*eq212-eq111.*eq203.*eq211-eq112.*eq201.*eq212-eq112.*eq202.*eq211-eq112.*eq203.*eq210+eq120.*eq202.*eq203.*2.0+eq121.*eq201.*eq203.*2.0;
eq4 = eq100.*eq212.^2+eq102.*eq211.^2+eq120.*eq202.^2+eq101.*eq211.*eq212.*2.0+eq102.*eq210.*eq212.*2.0-eq110.*eq202.*eq212-eq110.*eq203.*eq211-eq111.*eq201.*eq212-eq111.*eq202.*eq211-eq111.*eq203.*eq210-eq112.*eq200.*eq212-eq112.*eq201.*eq211-eq112.*eq202.*eq210+eq120.*eq201.*eq203.*2.0+eq121.*eq200.*eq203.*2.0+eq121.*eq201.*eq202.*2.0;
eq3 = eq101.*eq211.^2+eq121.*eq201.^2+eq100.*eq211.*eq212.*2.0+eq101.*eq210.*eq212.*2.0+eq102.*eq210.*eq211.*2.0-eq110.*eq201.*eq212-eq110.*eq202.*eq211-eq110.*eq203.*eq210-eq111.*eq200.*eq212-eq111.*eq201.*eq211-eq111.*eq202.*eq210-eq112.*eq200.*eq211-eq112.*eq201.*eq210+eq120.*eq200.*eq203.*2.0+eq120.*eq201.*eq202.*2.0+eq121.*eq200.*eq202.*2.0;
eq2 = eq100.*eq211.^2+eq102.*eq210.^2+eq120.*eq201.^2+eq100.*eq210.*eq212.*2.0+eq101.*eq210.*eq211.*2.0-eq110.*eq200.*eq212-eq110.*eq201.*eq211-eq110.*eq202.*eq210-eq111.*eq200.*eq211-eq111.*eq201.*eq210-eq112.*eq200.*eq210+eq120.*eq200.*eq202.*2.0+eq121.*eq200.*eq201.*2.0;
eq1 = eq101.*eq210.^2+eq121.*eq200.^2+eq100.*eq210.*eq211.*2.0-eq110.*eq200.*eq211-eq110.*eq201.*eq210-eq111.*eq200.*eq210+eq120.*eq200.*eq201.*2.0;
eq0 = eq100.*eq210.^2+eq120.*eq200.^2-eq110.*eq200.*eq210;

eq(:,:,1)=eq6; eq(:,:,2)=eq5; eq(:,:,3)=eq4; eq(:,:,4)=eq3; eq(:,:,5)=eq2; eq(:,:,6)=eq1; eq(:,:,7)=eq0;
options = optimoptions('lsqnonlin','Display','none','UseParallel',true);
parfor i =1:size(DD,2)
    res(i,:)= solve_eqn(DD(:,i),eq(:,i,:),eqa(:,i,:),eqb(:,i,:),a(:,i),b(:,i),c(:,i),d(:,i),options);   
end
end

function soln = solve_eqn(DD,eq,eqa,eqb,a,b,c,d,options)
eqa = reshape(eqa,size(eqa,1),size(eqa,3));
eqb = reshape(eqb,size(eqb,1),size(eqb,3));
idx = find(abs(DD)>1e-5);
id = randi(length(idx),1,max(10,ceil(0.1*length(idx))));
t= [];
for j= 1:length(id)
    eqn = [eq(id(j),:,1),eq(id(j),:,2),eq(id(j),:,3),eq(id(j),:,4),eq(id(j),:,5),eq(id(j),:,6),eq(id(j),:,7)];
    k2 = roots(eqn);
    k2 = k2(imag(k2)==0); % real roots to k2
    k2 = k2(abs(k2)<5);
    eq2 = repmat(eqb(id(j),:),length(k2),1);
    k1 = -(eq2(:,2).*k2.^3 + eq2(:,4).*k2.^2+eq2(:,6).*k2 + eq2(:,7))./(eq2(:,1).*k2.^2+eq2(:,3).*k2 + eq2(:,5));
    k1k2=[k1(:),k2(:)]';
    Jtinv = inv([a(id(j)),b(id(j));c(id(j)),d(id(j))]);
    x1x2 = Jtinv*k1k2;
    x1x2 = x1x2(:,abs(x1x2(1,:))< 5);
    x1x2 = x1x2(:,abs(x1x2(2,:))< 5);
    t = [t,x1x2];
end 
if isempty(t)
    x_init = [0,0];
    idx = 1:size(eqa,1);
else
[x_init, idx] = pick_k1k2(t,eqa,eqb,a,b,c,d);
end
soln = lsqnonlin(@(x)cost(x,eqa(idx,:),eqb(idx,:),a(idx),b(idx),c(idx),d(idx)),x_init,[],[],options);
end

function F = cost(x,eqa,eqb,a,b,c,d)
k1 = a.*repmat(x(1),length(a),1) + b.*repmat(x(2),length(a),1);
k2 = c.*repmat(x(1),length(a),1) + d.*repmat(x(2),length(a),1);
c1 = [k1.^2.*k2,k1.*k2.^2,k1.^2,k1.*k2,k2.^2,k1,k2,ones(size(k1))];
c2 = [k1.*k2.^2,k2.^3,k1.*k2,k2.^2,k1,k2,ones(size(k1))];
E1 = sum(eqa.*c1,2);
E2 = sum(eqb.*c2,2);
F = [E1(:);E2(:)];
end

function [x_init, idx] = pick_k1k2(t,eqa,eqb,a,b,c,d)
for j=1:size(t,2)
    k1 = a.*repmat(t(1,j),length(a),1) + b.*repmat(t(2,j),length(a),1);
    k2 = c.*repmat(t(1,j),length(a),1) + d.*repmat(t(2,j),length(a),1);
    c1 = [k1.^2.*k2,k1.*k2.^2,k1.^2,k1.*k2,k2.^2,k1,k2,ones(size(k1))];
    c2 = [k1.*k2.^2,k2.^3,k1.*k2,k2.^2,k1,k2,ones(size(k1))];
    res(j,:) = abs(sum(eqa.*c1,2)) + abs(sum(eqb.*c2,2));
    resj = res(j,:);
    id{j}= find(resj<10*median(resj));
    m(j)=median(res(j,id{j}));
end
[~, idx] = min(m);
x_init = t(:,idx)';
idx = id{idx};
end