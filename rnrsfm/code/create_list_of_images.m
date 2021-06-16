function [img_final,ref] = create_list_of_images(n,qgth,n1)
window = 1:n:size(qgth,1)/2;
for i= 1:length(window)
    img_list{i} = window(i):1:window(i)+n-1;
    
end
i = length(window);
if sum(img_list{i}> size(qgth,1)/2) >0 
    img_list{i}(img_list{i}>size(qgth,1)/2)=[];
end
if length(img_list{i}) <= n1
    img_list{i-1}= [img_list{i-1},img_list{i}];
    img_list{i}=[];
    for j = 1: length(window)-1
        img_final{j} = img_list{j};
    end
    
else
    img_final = img_list;
end

for i = 1:length(img_final)
    t = img_final{i};
    ref{i}= round(linspace(1,length(t),round(n/5)));
end