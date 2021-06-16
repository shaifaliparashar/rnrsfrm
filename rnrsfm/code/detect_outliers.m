function inlier = detect_outliers(Pr,P_grid,nng)
inlier = ones(size(Pr,1)/3,size(Pr,2));
% find relative scale
for i = 1: size(Pr,2)
    if Pr(1,i)~=0 && ~isnan(Pr(1,i)) && ~isinf(Pr(1,i))
        P1 = [Pr(1:3,i),P_grid(1:3,nng{i})];
        dist1 = [];
        for j = 2: size(P1,2)
            dist1 = [dist1,sqrt(sum((P1(:,j)-P1(:,1)).^2))];
        end
        % recscale P1 such that sum(dist1)=1
        P1_sc = P1./sum(dist1);
        dist1 = [];
        for j = 2: length(P1)
            dist1 = [dist1,sqrt(sum((P1_sc(:,j)-P1_sc(:,1)).^2))];
        end
        for t = 2: size(Pr,1)/3
            if Pr(3*(t-1)+1,i)~=0 && ~isnan(Pr(3*(t-1)+1,i)) && ~isinf(Pr(3*(t-1)+1,i))
                if ~sum(sum(isnan(P_grid(3*(t-1)+1:3*(t-1)+3,nng{i}))))
                    P2 = [Pr(3*(t-1)+1:3*(t-1)+3,i),P_grid(3*(t-1)+1:3*(t-1)+3,nng{i})];
                    [~,P2_sc,~] = absor(P2,P1_sc,'doScale',true);
                    dist2 =[];
                    for j = 2: length(P2_sc)
                        dist2 = [dist2,sqrt(sum((P2_sc(:,j)-P2_sc(:,1)).^2))];
                    end
                    dist(t-1,:)=dist2./dist1;
                end
            end
        end
        dist_med = median(dist'); % alpha_i
        dist_f = dist1;
        for t=2:size(Pr,1)/3
            Qj = dist_med(t-1)*P1_sc;
            P2 = [Pr(3*(t-1)+1:3*(t-1)+3,i),P_grid(3*(t-1)+1:3*(t-1)+3,nng{i})];
            [~,P2e,~] = absor(P2,Qj,'doScale',true);
            dist2 =[];
            for j = 2: length(P2_sc)
                dist2 = [dist2,sqrt(sum((P2e(:,j)-Qj(:,1)).^2))];
            end
            dist_f(t,:)=dist2;
        end
        distfn = dist_f(:,2:end);
        distf1 = dist_f(:,1);
        th = 0.1*mean(mean(distfn));
        %th = median(distfn(:));
        idf=find(distf1<=th);
        if length(idf)< round(length(distf1(:))/2)
           inlier(t,i) =0;
        end
    end
end



