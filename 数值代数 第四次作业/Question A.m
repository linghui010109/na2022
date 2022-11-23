x=linspace(0.99,1.01,101);
f=x.^8-8*x.^7+28*x.^6-56*x.^5+70*x.^4-56*x.^3+28*x.^2-8*x+1;
g=(((((((x-8).*x+28).*x-56).*x+70).*x-56).*x+28).*x-8).*x+1;
h=(x-1).^8;
plot(x,f,LineWidth=1);
hold on;
plot(x,g,LineWidth=1);
hold on;
plot(x,h,LineWidth=1);
