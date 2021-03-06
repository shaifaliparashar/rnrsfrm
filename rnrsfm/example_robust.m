% example script
clear all
close all
addpath(genpath('code'));


load Kinect_paper.mat

% partition images in to subsets of 35 images with 7 references
[list,ref] = create_list_of_images(35,qg,7);

 % 1. remove MAD correspondences
[qpn,qpn2,idd,rmv] = fix_mad_correspondences(qp,visb',list); % idd are flagged points
% mad corrected points
qna=zeros(size(qg));
visb2 = visb;
for i = 1: length(qpn)
    qn{i} = qpn{i};
    visb2(rmv{i},i)=0;
    idx = find(visb2(:,i)~=0);
    temp = inv(K)*qn{i}(:,idx);
    qna(2*i-1:2*i,idx) = temp(1:2,:);
    %visb2(idd{i},i)=0;
%     figure(1)
%     clf
%     hold on
%     plot(qg(2*i-1,:),qg(2*i,:),'*r')
%     plot(qna(2*i-1,:),qna(2*i,:),'og')
%     hold off
%     pause(0.1)
end


%  MAD + Multiple view reconstruction: 1 (Fast) , 2 (Resultant)
[Nalln1,Nalln2,ualln,valln] = mv_reconstruction_all(list,ref,qna,visb2,20);
P_mvn1 = median_shape_all(Nalln1,ualln,valln,list,ref,1e0,400,5);
P_mvn2 = median_shape_all(Nalln2,ualln,valln,list,ref,1e0,400,5);
P21 = interpolate_Pgth(P_mvn1,ualln,valln,qg);
P22 = interpolate_Pgth(P_mvn2,ualln,valln,qg);
[P1,err_p1] = compute_errors_P(P21,Pgth);
[P2,err_p2] = compute_errors_P(P22,Pgth);

350*mean(mean(err_p1'))
350*mean(mean(err_p2'))

% Outlier rejection
%for a very noisy dataset 
[in_det,in_pts] = detect_and_remove_outliers(qg,ualln,valln,P_mvn1,P21,list,20,idd);
sum(sum(in_det))/(1503*191)
sum(sum(in_pts))/(1503*191)

%inliers
in1 = in_det+in_pts;
in1(in1<2)=0; in1(in1>0)=1;

[in_det2,in_pts2] = detect_and_remove_outliers(qg,ualln,valln,P_mvn2,P22,list,20,idd);
sum(sum(in_det2))/(1503*191)
sum(sum(in_pts2))/(1503*191)

in2 = in_det2+in_pts2;
in2(in2<2)=0; in2(in2>0)=1;

350*mean(mean((err_p1.*in1)'))
350*mean(mean((err_p2.*in2)'))



