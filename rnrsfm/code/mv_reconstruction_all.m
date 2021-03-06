function [Nall1,Nall2,uall,vall]=mv_reconstruction_all(list,ref,q,visb,p)
i1 = length(list{1});
e= length(list);
i2 = length(list{e});

% compute  the normals for  each image subset with each reference in the 
Nall1 = []; Nall2 = []; uall = []; vall = [];
for k = 1:length(list) % for each subset of images
    k
    t1 = 2*(list{k}-1);
    tt = [t1+1;t1+2];
    q_n = q(tt(:),:);
    vis = visb(list{k},:);
    [u,v,u1,v1,visb2] = create_grid(q_n,vis,p);
    ua = [u(1,:);u1];
    va = [v(1,:);v1];
    N_res1 = []; N_res2 = []; uall = [uall;ua]; vall = [vall;va];
    for i = 1:length(ref{k})% for each reference in the subset
        a = 1:length(list{k});
        t = ref{k};
        temp = a(t(i)); a(t(i))=a(1); a(1)= temp;
        ua2 = ua(a,:); va2 = va(a,:);
        visb3 = visb2(a,:);
        [N_r1,N_r2]=mv_reconstruction(repmat(ua2(1,:),length(list{k})-1,1),...
            repmat(va2(1,:),length(list{k})-1,1),ua2(2:end,:),va2(2:end,:),visb3);
        N_res1=[N_res1,N_r1];N_res2=[N_res2,N_r2];
    end
    if i2>i1 && k~= length(list)
        N = [N_res1;zeros(3*(i2-i1),size(N_res1,2))]; 
        Nall1 = [Nall1,N];
        N = [N_res2;zeros(3*(i2-i1),size(N_res2,2))]; 
        Nall2 = [Nall2,N];
    elseif i2<i1 && k==length(list)
        N = [N_res1;zeros(3*(i1-i2),size(N_res1,2))]; 
        Nall1 = [Nall1,N]; 
        N = [N_res2;zeros(3*(i1-i2),size(N_res2,2))]; 
        Nall2 = [Nall2,N]; 
    else
        Nall1 = [Nall1,N_res1]; 
        Nall2 = [Nall2,N_res2]; 
    end
    
end