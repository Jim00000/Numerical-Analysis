function w = poisson(xl, xr, yb, yt, M, N)
%myFun - Description
%
% Syntax: w = poisson(xl, xr, yb, yt, M, N)
%
% Long description
    
    f=@(x,y) 0;
    g1=@(x) log(x.^2+1);
    g2=@(x) log(x.^2+4);
    g3=@(y) 2 * log(y);
    g4=@(y) log(y.^2+1);
    m = M + 1;n = N + 1; mn = m * n;
    h = (xr - xl) / M; h2 = h^2; k=(yt - yb) / N; k2 = k^2;
    x = xl + (0:M)*h;
    y = yb+(0:N)*k;
    A = zeros(mn,mn);b = zeros(mn,1);

    for i=2:m-1
        for j = 2:n-1
            A(i+(j-1)*m,i-1+(j-1)*m) = 1 / h2;
            A(i+(j-1)*m,i+1+(j-1)*m) = 1 / h2;
            A(i+(j-1)*m,i+(j-1)*m) = - 2 / h2 - 2 / k2;
            A(i+(j-1)*m,i+(j-2)*m) = 1 / k2;
            A(i+(j-1)*m,i+j*m) = 1 / k2;
            b(i+(j-1)*m) = f(x(i),y(j));
        end
    end

    for i = 1:m
        j = 1;
        A(i+(j-1)*m, i+(j-1)*m) = 1;
        b(i+(j-1)*m) = g1(x(i));
        j = n;
        A(i+(j-1)*m, i+(j-1)*m) = 1;
        b(i+(j-1)*m) = g2(x(i));
    end

    for j = 2:n-1
        i = 1;
        A(i+(j-1)*m, i+(j-1)*m) = 1;
        b(i+(j-1)*m) = g3(y(j));
        i = m;
        A(i+(j-1)*m, i+(j-1)*m) = 1;
        b(i+(j-1)*m) = g4(y(j));
    end

    v = A\b;
    w = reshape(v(1:mn),m,n);

    mesh(x,y,w');

end
