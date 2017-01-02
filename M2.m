function x2=M2(a)
x2=cumsum(a(end:-1:1,:),1);
x2=x2(end:-1:1,:);