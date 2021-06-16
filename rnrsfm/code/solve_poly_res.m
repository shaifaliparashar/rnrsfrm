function res = solve_poly_res(a,b,c,d,t1,t2,e,e1,u,u1,v,v1)

g1120 = a.^2.*e;
g1111 = a.*b.*e.*2;
g1102 = b.^2.*e;
g1110 = (a.*e.*t1.*2) -(a.*u.*2);
g1101 = (b.*e.*t1.*2) -(b.*u.*2);
g1100 = (e.*t1.^2) -(t1.*u.*2) + 1;

g1220 = a.*c.*e;
g1211 = (a.*d.*e) + (b.*c.*e);
g1202 = b.*d.*e;
g1210 = -(a.*v)-(c.*u)+(a.*e.*t2)+(c.*e.*t1);
g1201 = -(b.*v)-(d.*u)+(b.*e.*t2)+(d.*e.*t1);
g1200 = -(t2.*u)-(t1.*v)+(e.*t1.*t2);

g2220 = c.^2.*e;
g2211 = c.*d.*e.*2;
g2202 = d.^2.*e;
g2210 = (c.*e.*t2.*2) -(c.*v.*2);
g2201 = (d.*e.*t2.*2) -(d.*v.*2);
g2200 = (e.*t2.^2) -(t2.*v.*2) + 1;

% pullback metric tensor at image i: g11b, g12b, g22b

g11b20 = a.^2.*e1;
g11b11 = a.*b.*e1.*2;
g11b02 = b.^2.*e1;
g11b10 = -(a.^2.*u1.*2)-(a.*b.*v1.*2);
g11b01 = -(b.^2.*v1.*2)-(a.*b.*u1.*2);
g11b00 = a.^2+b.^2;

g12b20 = a.*c.*e1;
g12b11 = (a.*d.*e1) + (b.*c.*e1);
g12b02 = b.*d.*e1;
g12b10 = -(v1.*((a.*d) + (b.*c)))-(a.*c.*u1.*2);
g12b01 = -(u1.*((a.*d) + (b.*c)))-(b.*d.*v1.*2);
g12b00 = a.*c+b.*d;

g22b20 = c.^2.*e1;
g22b11 = c.*d.*e1.*2;
g22b02 = d.^2.*e1;
g22b10 = -(c.^2.*u1.*2)-(c.*d.*v1.*2);
g22b01 = -(d.^2.*v1.*2)-(c.*d.*u1.*2);
g22b00 = c.^2+d.^2;

%eq1 = g11*g22b - g11b*g22
eq130 = g1110.*g22b20 + g1120.*g22b10 - g11b10.*g2220 - g11b20.*g2210;
eq121 = g1101.*g22b20 + g1110.*g22b11 + g1111.*g22b10 + g1120.*g22b01 - g11b01.*g2220 - g11b10.*g2211 - g11b11.*g2210 - g11b20.*g2201;
eq112 = g1101.*g22b11 + g1102.*g22b10 + g1110.*g22b02 + g1111.*g22b01 - g11b01.*g2211 - g11b02.*g2210 - g11b10.*g2202 - g11b11.*g2201;
eq103 = g1101.*g22b02 + g1102.*g22b01 - g11b01.*g2202 - g11b02.*g2201;

eq120 = g1100.*g22b20 + g1110.*g22b10 + g1120.*g22b00 - g11b00.*g2220 - g11b10.*g2210 - g11b20.*g2200;
eq111 = g1100.*g22b11 + g1101.*g22b10 + g1110.*g22b01 + g1111.*g22b00 - g11b00.*g2211 - g11b01.*g2210 - g11b10.*g2201 - g11b11.*g2200;
eq102 = g1100.*g22b02 + g1101.*g22b01 + g1102.*g22b00 - g11b00.*g2202 - g11b01.*g2201 - g11b02.*g2200;

eq110 = g1100.*g22b10 + g1110.*g22b00 - g11b00.*g2210 - g11b10.*g2200;
eq101 = g1100.*g22b01 + g1101.*g22b00 - g11b00.*g2201 - g11b01.*g2200;

eq100 = g1100.*g22b00 - g11b00.*g2200;

%eq2 = g12*g22b - g12b*g22
eq230 = g1210.*g22b20 + g1220.*g22b10 - g12b10.*g2220 - g12b20.*g2210;
eq221 = g1201.*g22b20 + g1210.*g22b11 + g1211.*g22b10 + g1220.*g22b01 - g12b01.*g2220 - g12b10.*g2211 - g12b11.*g2210 - g12b20.*g2201;
eq212 = g1201.*g22b11 + g1202.*g22b10 + g1210.*g22b02 + g1211.*g22b01 - g12b01.*g2211 - g12b02.*g2210 - g12b10.*g2202 - g12b11.*g2201;
eq203 = g1201.*g22b02 + g1202.*g22b01 - g12b01.*g2202 - g12b02.*g2201;

eq220 = g1200.*g22b20 + g1210.*g22b10 + g1220.*g22b00 - g12b00.*g2220 - g12b10.*g2210 - g12b20.*g2200;
eq211 = g1200.*g22b11 + g1201.*g22b10 + g1210.*g22b01 + g1211.*g22b00 - g12b00.*g2211 - g12b01.*g2210 - g12b10.*g2201 - g12b11.*g2200;
eq202 = g1200.*g22b02 + g1201.*g22b01 + g1202.*g22b00 - g12b00.*g2202 - g12b01.*g2201 - g12b02.*g2200;

eq210 = g1200.*g22b10 + g1210.*g22b00 - g12b00.*g2210 - g12b10.*g2200;
eq201 = g1200.*g22b01 + g1201.*g22b00 - g12b00.*g2201 - g12b01.*g2200;

eq200 = g1200.*g22b00 - g12b00.*g2200;

eq1(:,:,1)=eq130;eq1(:,:,2)=eq121;eq1(:,:,3)=eq112;eq1(:,:,4)=eq103;eq1(:,:,5)=eq120;
eq1(:,:,6)=eq111;eq1(:,:,7)=eq102;eq1(:,:,8)=eq110;eq1(:,:,9)=eq101;eq1(:,:,10)=eq100;
eq2(:,:,1)=eq230;eq2(:,:,2)=eq221;eq2(:,:,3)=eq212;eq2(:,:,4)=eq203;eq2(:,:,5)=eq220;
eq2(:,:,6)=eq211;eq2(:,:,7)=eq202;eq2(:,:,8)=eq210;eq2(:,:,9)=eq201;eq2(:,:,10)=eq200;
th1 = 0.01*median(abs(eq1),3);th2 = 0.01*median(abs(eq2),3);
eqa = eq1; eqb = eq2; eq1(abs(eq1)<th1)=0; eq2(abs(eq2)<th2)=0;
options = optimoptions('lsqnonlin','Display','none','UseParallel',true);
parfor i=1:size(eq1,2)
    res(i,:)=compute_k1k2(eq1(:,i,:),eq2(:,i,:),eqa(:,i,:),eqb(:,i,:),options);
end
end

function F = cost2(x,eq1,eq2)
x1 = repmat(x(1),size(eq1,1),1); x2 = repmat(x(2),size(eq1,1),1);
cc = [x1.^3,x1.^2.*x2,x1.*x2.^2,x2.^3,x1.^2,x1.*x2,x2.^2,x1,x2,ones(size(x1))];
ea = reshape(eq1,size(eq1,1),size(eq1,3));eb = reshape(eq2,size(eq1,1),size(eq1,3));
E1 = sum(ea.*cc,2);
E2 = sum(eb.*cc,2);
F= [E1(:);E2(:)];
end

function k1k2 = compute_k1k2(eq1,eq2,eqa,eqb,options)
%tic
id = randi(size(eq1,1),1,max(10,ceil(0.1*length(size(eq1,1)))));
x= [];
for j= 1:length(id)
    det_syl = solve_recon(eq1(id(j),:,:),eq2(id(j),:,:));
    k1 = roots(det_syl);
    k1 = k1(imag(k1)==0); % real roots to k1
    k1 = k1(abs(k1)< 5);
    ea = reshape(eqa(id(j),:,:),1,size(eqa,3));
    eb = reshape(eqb(id(j),:,:),1,size(eqb,3));
    for k = 1:length(k1)
       k2 = roots([ea(4),ea(3)*k1(k)+ea(7),ea(2)*k1(k)^2+ea(6)*k1(k)+ea(9),ea(10)+ea(8)*k1(k)+ea(5)*k1(k)^2+ea(1)*k1(k)^3]);
       k2 = k2(imag(k2)==0); % real roots to k1
       k2 = k2(abs(k2)< 2);
       s = [repmat(k1(k),length(k2),1),k2];
       x=[x;s];
    end
     %t=[t;k1(:)];
end
%toc
if isempty(x)
    x_init = [0,0];
    idd = 1:size(eq1,1);
else
[x_init, idd] = pick_k1k2(x,eqa,eqb);
end
 k1k2 = lsqnonlin(@(y)cost2(y,eqa(idd,:,:),eqb(idd,:,:)),[x_init(1),x_init(2)],[],[],options);
end

function [x_init, idx] = pick_k1k2(t,eqa,eqb)
eqa = reshape(eqa,size(eqa,1),size(eqa,3));eqb = reshape(eqb,size(eqb,1),size(eqb,3));
for j=1:size(t,1)
    k1 = repmat(t(j,1),size(eqa,1),1);
    k2 = repmat(t(j,1),size(eqa,1),1);
    cc = [k1.^3,k1.^2.*k2,k1.*k2.^2,k2.^3,k1.^2,k1.*k2,k2.^2,k1,k2,ones(size(k1))];
    res(j,:) = abs(sum(eqa.*cc,2)) + abs(sum(eqb.*cc,2));
    resj = res(j,:);
    id{j}= find(resj<10*median(resj));
    m(j)=median(res(j,id{j}));
end
[~, idx] = min(m);
x_init = t(idx,:);
idx = id{idx};
end


