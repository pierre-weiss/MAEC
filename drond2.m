function d=drond2(im)
d=zeros(size(im));
d(:,1:end-1)=im(:,2:end)-im(:,1:end-1);
