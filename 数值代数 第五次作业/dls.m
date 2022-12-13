function[a,s]=dls(x,y)%x,y是横坐标
[~,col]=size(x);
G=zeros(3,3);
c=zeros(3,1);
a=zeros(3,1);

for i=1:3
    for j=1:3
        if (j==1&&i==1)
            G(1,1)=col;
        else
            G(i,j)=sum(x.^(i-1).*x.^(j-1));
        end
    end
end
for i=1:3
    c(i,1)=sum(x.^(i-1).*y);
end
a=G\c;
s=cond(G,2);
