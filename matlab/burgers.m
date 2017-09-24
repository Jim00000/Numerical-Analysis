function w = burgers(xl,xr,tb,te,M,N)
%myFun - Description
%
% Syntax: w = burgers(xl,xr,tb,te,M,N)
%
% Long description
    alf = 5;
    bet = 4;
    D = 0.05;
    f = @(x) 2 * D * bet * pi * sin(pi * x) ./ (alf + bet * cos(pi * x));
    l = @(t) 0 * t;
    r = @(t) 0 * t;
    h = (xr - xl) / M; k = (te - tb) / N; m = M + 1; n = N;
    sigma = D * k / (h * h);

    % initial conditions
    w(:,1) = f(xl + (0:M) * h);
    w1 = w;

    for j=1:n
        % Newton iteration
        for it = 1:3
            DF1 = diag(1 + 2 * sigma * ones(m, 1)) + diag(-sigma * ones(m-1, 1), 1) + diag(-sigma * ones(m-1, 1), -1);
            DF2 = diag([0;k*w1(3:m)/(2*h);0]) - diag([0;k * w1(1:(m - 2)) / (2 * h);0])...
                + diag([0;k * w1(2:m - 1) / (2 * h)], 1) - diag([k * w1(2:m - 1) / (2 * h);0], -1);
            DF = DF1 + DF2;
            F = -w(:,j) + (DF1 + DF2 / 2) * w1;
            % Dirichlet conditions for DF
            DF(1,:) = [1 zeros(1,m-1)];
            % Dirichlet conditions for F
            F(1) = w1(1) - l(j); F(m) = w1(m) - r(j);
            w1 = w1 -DF\F;
        end
        w(:, j + 1) = w1;
    end

    x = xl + (0:M) * h; t = tb + (0:n) * k;
    % 3-D plot of solution w
    mesh(x,t,w');

end