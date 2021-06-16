%INPUT
% P : feature points
% C : Control Center
% L : Warp Parameters
% lambda: Inernal Regulrization Parameter

%OUTPUT
% I,J,M,N : The differential invariants.


% function [I J M N] = schwarzian(P, L, C, lambda,choice,MM00,MM10,MM01,MM20,MM02,MM11)
function [I J M N jI jJ jM jN] = schwarzian(bbs,ctrlpts,X2,Y2,choice)

DAdetabydu1  = bbs_eval(bbs, ctrlpts, X2, Y2, 1,0)';
DAdetabydv1  = bbs_eval(bbs, ctrlpts, X2, Y2, 0,1)';
DAdetabydu2  = bbs_eval(bbs, ctrlpts, X2, Y2, 2,0)';
DAdetabydv2  = bbs_eval(bbs, ctrlpts, X2, Y2, 0,2)';
DAdetabydudv = bbs_eval(bbs, ctrlpts, X2, Y2, 1,1)';
coloc_du = bbs_coloc_deriv(bbs, X2, Y2,1,0);
coloc_dv = bbs_coloc_deriv(bbs, X2, Y2,0,1);
coloc_duu= bbs_coloc_deriv(bbs, X2, Y2,2,0);
coloc_duv= bbs_coloc_deriv(bbs, X2, Y2,1,1);
coloc_dvv= bbs_coloc_deriv(bbs, X2, Y2,0,2);
if(strcmp(choice,'den'))
    
    I = (DAdetabydu2(:,1).*DAdetabydu1(:,2) - DAdetabydu2(:,2).*DAdetabydu1(:,1)) ...
        ./ (DAdetabydu1(:,1).*DAdetabydv1(:,2) - DAdetabydu1(:,2).*DAdetabydv1(:,1));
    
    J = (DAdetabydv2(:,2).*DAdetabydv1(:,1) - DAdetabydv2(:,1).*DAdetabydv1(:,2)) ...
        ./ (DAdetabydu1(:,1).*DAdetabydv1(:,2) - DAdetabydu1(:,2).*DAdetabydv1(:,1));
    
    M = (DAdetabydu2(:,1).*DAdetabydv1(:,2) - DAdetabydu2(:,2).*DAdetabydv1(:,1) + ...
        2.*(DAdetabydudv(:,1).*DAdetabydu1(:,2) - DAdetabydudv(:,2).*DAdetabydu1(:,1))) ...
        ./ 3.*((DAdetabydu1(:,1).*DAdetabydv1(:,2) - DAdetabydu1(:,2).*DAdetabydv1(:,1)));
    
    N = (DAdetabydv2(:,2).*DAdetabydu1(:,1) - DAdetabydv2(:,1).*DAdetabydu1(:,2) + ...
        2.*(DAdetabydudv(:,2).*DAdetabydv1(:,1) - DAdetabydudv(:,1).*DAdetabydv1(:,2))) ...
        ./ 3.*((DAdetabydu1(:,1).*DAdetabydv1(:,2) - DAdetabydu1(:,2).*DAdetabydv1(:,1)));
    
    

elseif (strcmp(choice, 'noden'))
    
    I = (DAdetabydu2(:,1).*DAdetabydu1(:,2) - DAdetabydu2(:,2).*DAdetabydu1(:,1));
    
    J = (DAdetabydv2(:,2).*DAdetabydv1(:,1) - DAdetabydv2(:,1).*DAdetabydv1(:,2));
    
    M = (DAdetabydu2(:,1).*DAdetabydv1(:,2) - DAdetabydu2(:,2).*DAdetabydv1(:,1) + ...
        2.*(DAdetabydudv(:,1).*DAdetabydu1(:,2) - DAdetabydudv(:,2).*DAdetabydu1(:,1)));
    
    N = (DAdetabydv2(:,2).*DAdetabydu1(:,1) - DAdetabydv2(:,1).*DAdetabydu1(:,2) + ...
        2.*(DAdetabydudv(:,2).*DAdetabydv1(:,1) - DAdetabydudv(:,1).*DAdetabydv1(:,2)));
    
    D12u=DAdetabydu1';
    D12v=DAdetabydv1';
    D12uv=DAdetabydudv';
    D12uu=DAdetabydu2';
    D12vv=DAdetabydv2';
    nparam=size(ctrlpts,2);
    zerosm=zeros(size(coloc_du));
    r = size(zerosm,1); c = size(zerosm,2);

    A = bsxfun(@times,[D12u(1,:)';D12u(2,:)';D12v(1,:)';D12v(2,:)'],[coloc_duu,coloc_duv,coloc_dvv;coloc_duu,coloc_duv,coloc_dvv;coloc_duu,coloc_duv,coloc_dvv;coloc_duu,coloc_duv,coloc_dvv]);
    B = bsxfun(@times,[D12uu(1,:)';D12uu(2,:)';D12uv(1,:)';D12uv(2,:)';D12vv(1,:)';D12vv(2,:)'],[coloc_du,coloc_dv;coloc_du,coloc_dv;coloc_du,coloc_dv;coloc_du,coloc_dv;coloc_du,coloc_dv;coloc_du,coloc_dv]);

    
    jI= [A(r+1:2*r,1:c)-B(r+1:2*r,1:c),B(1:r,1:c)-A(1:r,1:c)];
    jJ= [B(5*r+1:6*r,c+1:2*c)-A(3*r+1:4*r,2*c+1:3*c),A(2*r+1:3*r,2*c+1:3*c)-B(4*r+1:5*r,c+1:2*c)];
    
    jM= [A(3*r+1:4*r,1:c)-B(r+1:2*r,c+1:2*c)+2*A(r+1:2*r,c+1:2*c)-2*B(3*r+1:4*r,1:c),...
         B(1:r,c+1:2*c)-A(2*r+1:3*r,1:c)+2*B(2*r+1:3*r,1:c)-2*A(1:r,c+1:2*c)];
    jN=[B(5*r+1:6*r,1:c)-A(r+1:2*r,2*c+1:3*c)-2*A(3*r+1:4*r,c+1:2*c)+2*B(3*r+1:4*r,c+1:2*c),...
        A(1:r,2*c+1:3*c)-B(4*r+1:5*r,1:c)-2*B(2*r+1:3*r,c+1:2*c)+2*A(2*r+1:3*r,c+1:2*c)];
    
    
    
    
end
end


