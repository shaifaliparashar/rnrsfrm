function [N1,N2]=mv_reconstruction(I1u,I1v,I2u,I2v,visb)
% [N,u_all,v_all,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uva,H21uvb]
u_all = [I1u(1,:);I2u];
v_all = [I1v(1,:);I2v];
[J21a,J21b,J21c,J21d,H21uva,H21uvb] = create_warps(I1u,I1v,I2u,I2v);
d = J21a.*J21d - J21b.*J21c;
J12a = J21d./d; J12b = -J21b./d; J12c = -J21c./d; J12d = J21a./d;


a = J21a; b = J21b; c = J21c; d = J21d; 
t1 = -(J12b.*H21uva + J12d.*H21uvb); t2 = -(J12a.*H21uva + J12c.*H21uvb);

e = 1+ I2u.^2 + I2v.^2; u = I2u; v = I2v;
e1 = 1+ I1u.^2 + I1v.^2; u1 = I1u; v1 = I1v;

qg = zeros(2*size(u_all,1),size(u_all,2));
qg(1:2:end,:)=u_all;
qg(2:2:end,:)=v_all;
% FAST SUBSTITUTION METHOD
res = solve_poly_subs(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
% obtain normals and depth
[P, N1] = create_depth_normals(res,a,b,c,d,t1,t2,visb,I1u,I2u,I1v,I2v,qg);


% FAST RESULTANT METHOD
res = solve_poly_res(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
% obtain normals and depth
[P, N2] = create_depth_normals(res,a,b,c,d,t1,t2,visb,I1u,I2u,I1v,I2v,qg);


