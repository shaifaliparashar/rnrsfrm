% example script
clear all
close all
addpath(genpath('code'));


load Kinect_paper.mat

grid_size = 20; % size of grid

[I1u,I1v,I2u,I2v,~] = create_grid(qg,visb',grid_size);
[J21a,J21b,J21c,J21d,H21uva,H21uvb] = create_warps(I1u,I1v,I2u,I2v);
d = J21a.*J21d - J21b.*J21c;
J12a = J21d./d; J12b = -J21b./d; J12c = -J21c./d; J12d = J21a./d;


a = J21a; b = J21b; c = J21c; d = J21d; 
t1 = -(J12b.*H21uva + J12d.*H21uvb); t2 = -(J12a.*H21uva + J12c.*H21uvb);

e = 1+ I2u.^2 + I2v.^2; u = I2u; v = I2v;
e1 = 1+ I1u.^2 + I1v.^2; u1 = I1u; v1 = I1v;

% FAST SUBSTITUTION METHOD
res = solve_poly_subs(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
% obtain normals and depth
[P, N] = create_depth_normals(res,a,b,c,d,t1,t2,visb',I1u,I2u,I1v,I2v,qg);
% compute erros
[err_p, err_n] = compute_errors(P,N,Pgth,Ngth);


% FAST RESULTANT METHOD
res1 = solve_poly_res(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
% obtain normals and depth
[P1, N1] = create_depth_normals(res1,a,b,c,d,t1,t2,visb',I1u,I2u,I1v,I2v,qg);
% compute erros
[err_p1, err_n1] = compute_errors(P1,N1,Pgth,Ngth);


mean(mean(err_n'))
mean(mean(err_n1'))
350*mean(mean(err_p'))
350*mean(mean(err_p1'))


