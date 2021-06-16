function [q_new, q_new2, ids,rmv] = fix_mad_correspondences(q,visb,list)
% q_new : uncorrected  points
% q_new2 : points corrected using a warp (use them if the registration is bad)
% ids : flagged points
% rmv : removed points
th = .001*sqrt(640^2+480^2);% set this according to the image size
er = 1e-4;
t = 5;
nC = 20;
q_new = q;
q_new2 = q;
for i = 1:length(list)
    s = list{i}(1);
    for j = 2: length(list{i})
        
        e = list{i}(j);
        idx = find(visb(s,:)+visb(e,:)==2);
        if ~isempty(idx)
            q1 = q{s}(1:2,idx);  q2 = q{e}(1:2,idx);
            umin = min(q{s}(1,:))-t; umax = max(q{s}(1,:))+t;
            vmin = min(q{s}(2,:))-t; vmax = max(q{s}(2,:))+t;
            bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
            coloc = bbs_coloc(bbs, q{s}(1,idx), q{s}(2,idx));
            lambdas = er*ones(nC-3, nC-3);
            bending = bbs_bending(bbs, lambdas);
            cpts = (coloc'*coloc + bending) \ (coloc'*q{e}(1:2,idx)');
            cpts = cpts';
            
            qw2 = bbs_eval(bbs,cpts,q1(1,:)',q1(2,:)',0,0);
            error = sqrt(mean((qw2(1,:)-q2(1,:)).^2+(qw2(2,:)-q2(2,:)).^2));
            m = sum(abs(qw2 - q2(1:2,:)));
            mad_score = 4.5*median(m);
            median_old = 1e5;
            id = find(m < mad_score);
            idd = find(m >= mad_score);
            rmv{e}=[];
            while  abs(median_old - median(m))>th
                median_old = median(m);
                coloc = bbs_coloc(bbs, q1(1,id), q1(2,id));
                lambdas = er*ones(nC-3, nC-3);
                bending = bbs_bending(bbs, lambdas);
                cpts = (coloc'*coloc + bending) \ (coloc'*q2(1:2,id)');
                cpts = cpts';
                qw2 = bbs_eval(bbs,cpts,q1(1,id)',q1(2,id)',0,0);
                error = sqrt(mean((qw2(1,:)-q2(1,id)).^2+(qw2(2,:)-q2(2,id)).^2));
                m = sum(abs(qw2 - q2(1:2,id)));
                mad_score = 4.5*median(m);
                if median(m) < median_old
                    id2 = find(m < mad_score);
                    idd2 = find(m >= 4.5*median(m));
                    idd = [idd,id(idd2)];
                    id = id(id2);
                    idd = unique(idd);      
                end
                qw2 = bbs_eval(bbs,cpts,q1(1,:)',q1(2,:)',0,0);
            end
            q_new{e}(1:2,idx) = q2;
            q_new{e}(1:2,idd) = qw2(:,idd);
            q_new2{e}(1:2,idx) = qw2;
            ids{e}= idx(idd);
            m = sum(abs(qw2(:,idd) - q2(1:2,idd)));
            rm = find(m > 4.5*median(m));
            rmv{e} = [rmv{e},idd(rm)];
            rmv{e} = idx(unique(rmv{e}));
        else
            q_new{e}(1:2,:) =  q{e}(1:2,:);
            q_new2{e}(1:2,:) =  q{e}(1:2,:);
            ids{e}= [];
            rmv{e} = [];
        end
    end
end
