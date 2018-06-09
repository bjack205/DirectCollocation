function xdot = invertedPendulumDynamics(X,U,p)
    Xdot = zeros(size(X));
    for i = 1:size(X,2)
        x = X(1,i);
        v = X(2,i);
        th = X(3,i);
        om = X(4,i);
        u = U(:,i);

        A = [p.M+p.m         , p.m*p.ell*cos(th);
             p.m*p.ell*cos(th), p.I+p.m*p.ell^2];
        c = [-p.b*v + p.m*p.ell*om^2*sin(th) + u, -p.m*p.g*p.ell*sin(th)];
        Ainv = inv(A);
        xdot(:,i) = [v;  Ainv(1,1)*c(1)+Ainv(1,2)*c(2);
                om; Ainv(2,1)*c(1)+Ainv(2,2)*c(2)];
    end
end
