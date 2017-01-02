
function d=drond1T(im)
d=zeros(size(im));
d(2:end-1,:)=im(1:end-2,:)-im(2:end-1,:);
d(1,:)=-im(1,:);
d(end,:)=im(end-1,:);

