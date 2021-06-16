function [err_p, err_n] = compute_errors(P,N,Pgth,Ngth)

for i=1:size(P,1)/3
    [~,P(3*i-2:3*i,:),~]= absor(P(3*i-2:3*i,:),Pgth(3*i-2:3*i,:),'doScale',true);
    scale = max(max(Pgth(3*i-2:3*i,:)')-min(Pgth(3*i-2:3*i,:)'));
    err_p(i,:) = sqrt(mean((Pgth(3*i-2:3*i,:)-P(3*i-2:3*i,:)).^2))/scale;
    err = acosd(sum(N(3*i-2:3*i,:).*Ngth(3*i-2:3*i,:)));
    err(err>=90) = 180 - err(err>=90);
    err_n(i,:) = err;
end