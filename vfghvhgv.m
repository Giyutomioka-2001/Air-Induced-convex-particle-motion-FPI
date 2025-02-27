x = [0 0.2 0.6 0.8 1.0]
y = [20 30 40 60 80]

n= length(x)- 1
xp = 35;
sm =0;
for i = 1:n+1
    pr = 1;
    for j=1:n+1
        if j~=i
            pr = pr*(xp-x(j))/(x(i)-x(j));
        end
    end
    sm= sm+y(i)*pr;
end
yp=sm


