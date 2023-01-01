x=-1:0.02:1;
y=1./(1+25.*x.*x);
plot(A6(:,1),A6(:,2),A11(:,1),A11(:,2),A21(:,1),A21(:,2),A41(:,1),A41(:,2),A81(:,1),A81(:,2),x,y,'g');
plot(Anatural6(:,1),Anatural6(:,2),Anatural11(:,1),Anatural11(:,2),Anatural21(:,1),Anatural21(:,2),Anatural41(:,1),Anatural41(:,2),Anatural81(:,1),Anatural81(:,2),x,y,'g');
plot(Anot6(:,1),Anot6(:,2),Anot11(:,1),Anot11(:,2),Anot21(:,1),Anot21(:,2),Anot41(:,1),Anot41(:,2),Anot81(:,1),Anot81(:,2),x,y,'g');
