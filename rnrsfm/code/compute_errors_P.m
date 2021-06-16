function [P,err_p] = compute_errors_P(P,Pgth)

for i=1:size(P,1)/3
    [~,P(3*i-2:3*i,:),~]= absor(P(3*i-2:3*i,:),Pgth(3*i-2:3*i,:),'doScale',true);
    scale = max(max(Pgth(3*i-2:3*i,:)')-min(Pgth(3*i-2:3*i,:)'));
    err_p(i,:) = sqrt(mean((Pgth(3*i-2:3*i,:)-P(3*i-2:3*i,:)).^2))/scale;
end