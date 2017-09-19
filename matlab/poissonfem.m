function w = poissonfem(xl,xr,yb,yt,M,N)
%myFun - Description
%
% Syntax: w = poissonfem(xl,xr,yb,yt,M,N)
%
% Long description
    
    f = @(x, y) 0;  % input function
    r = @(x, y) 0;
    g1 = @(x) log(x.^2 + 1);    % bottom
    g2 = @(x) log(x.^2 + 4);    % top
    g3 = @(y) 2 * log(y);       % left
    g4 = @(y) log(y.^2 + 1);    % right
    m = M + 1;n = N + 1; mn = m * n;
    h = (xr - xl) / M; h2 = h^2; k=(yt - yb) / N; k2 = k^2; hk = h * k;
    x = xl + (0:M)*h;
    y = yb+(0:N)*k;
    A = zeros(mn,mn);
    b = zeros(mn,1);

    B1x = @(i) x(i) - 2 * h / 3;
    B1y = @(j) y(j) - k / 3;
    B2x = @(i) x(i) - h / 3;
    B2y = @(j) y(j) - 2 * k / 3;
    B3x = @(i) x(i) + h / 3;
    B3y = @(j) y(j) - k / 3;
    B4x = @(i) x(i) + 2 * h / 3;
    B4y = @(j) y(j) + k / 3;
    B5x = @(i) x(i) + h / 3;
    B5y = @(j) y(j) + 2 * k / 3;
    B6x = @(i) x(i) - h / 3;
    B6y = @(j) y(j) + k / 3;


    % interior points
    for i=2:m-1 
        for j=2:n-1
            rsum =  r(B1x(i),B1y(j)) + r(B2x(i),B2y(j)) + r(B3x(i),B3y(j)) + ...
                    r(B4x(i),B4y(j)) + r(B5x(i),B5y(j)) + r(B6x(i),B6y(j));
            fsum =  f(B1x(i),B1y(j)) + f(B2x(i),B2y(j)) + f(B3x(i),B3y(j)) + ...
                    f(B4x(i),B4y(j)) + f(B5x(i),B5y(j)) + f(B6x(i),B6y(j));
            A(i + (j - 1) * m, i + (j - 1) * m) = 2 * (h2 + k2) / (hk) - hk * rsum / 18;
            A(i + (j - 1) * m, i - 1 + (j - 1) * m) = - k / h - hk * (r(B1x(i),B1y(j)) + r(B6x(i),B6y(j))) / 18;
            A(i + (j - 1) * m, i - 1 + (j - 2) * m) = - hk * (r(B1x(i),B1y(j)) + r(B2x(i),B2y(j))) / 18;
            A(i + (j - 1) * m, i + (j - 2) * m) = - h / k - hk * (r(B2x(i),B2y(j)) + r(B3x(i),B3y(j))) / 18;
            A(i + (j - 1) * m, i + 1 + (j - 1) * m) = - k / h - hk * (r(B3x(i),B3y(j)) + r(B4x(i),B4y(j))) / 18;
            A(i + (j - 1) * m, i + 1 + j * m) = - hk * (r(B4x(i),B4y(j)) + r(B5x(i),B5y(j))) / 18;
            A(i + (j - 1) * m, i + j * m) = - h / k - hk * (r(B5x(i),B5y(j)) + r(B6x(i),B6y(j))) / 18;
            b(i + (j - 1) * m) = - h * k * fsum / 6; 
        end
    end

    % boundary points
    for i = 1:m
        j = 1;
        A(i+(j-1)*m, i+(j-1)*m) = 1;
        b(i+(j-1)*m) = g1(x(i));
        j = n;
        A(i+(j-1)*m, i+(j-1)*m) = 1;
        b(i+(j-1)*m) = g2(x(i));
    end

    % boundary points
    for j = 2:n-1
        i = 1;
        A(i+(j-1)*m, i+(j-1)*m) = 1;
        b(i+(j-1)*m) = g3(y(j));
        i = m;
        A(i+(j-1)*m, i+(j-1)*m) = 1;
        b(i+(j-1)*m) = g4(y(j));
    end

    v = A \ b;
    w = reshape(v(1:mn), m, n);
    mesh(x, y, w');

        
end