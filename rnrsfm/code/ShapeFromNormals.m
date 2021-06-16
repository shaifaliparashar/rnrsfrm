function [ctrpts3Dn]=ShapeFromNormals(bbs,coloc,bending,m1,n1)

% warp derivatives
%wu=bbs_eval(bbs, ctrlpts, m1(1,:)', m1(2,:)',1,0);
%wv=bbs_eval(bbs, ctrlpts, m1(1,:)', m1(2,:)',0,1);
coloc_du = bbs_coloc_deriv(bbs, m1(1,:)', m1(2,:)', 1, 0);
coloc_dv = bbs_coloc_deriv(bbs, m1(1,:)', m1(2,:)', 0, 1);


npts=size(m1,2);
nparams=size(coloc,2);

M=sparse(2*npts+1,nparams);
B=sparse(2*npts,1);

for i=1:npts
n=n1(:,i);
etat=[m1(:,i)];
etatu=[1;0;0];
etatv=[0;1;0];
crow=coloc(i,:);
crow_du=coloc_du(i,:);
crow_dv=coloc_dv(i,:);
% normal crow_du * ctr *etat + crow *etatu
M(i,1:nparams)=(n'*(etat*crow_du+ etatu*crow));%./(etatu'*etatu+etat'*etat);    
M(npts+i,1:nparams)=n'*(etat*crow_dv+ etatv*crow);%./(etatv'*etatv+etat'*etat);    
end

B(2*npts+1)=1;
M(2*npts+1,:)=ones(1,nparams);
%[U,D,V]=svds((M'*M+bending),1,0);
X=(M'*M+bending)\(M'*B);
%X=V(:,end);

ctrpts3Dn=reshape(X,nparams,1)';


end


