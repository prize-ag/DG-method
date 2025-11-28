function xi = newton_burgers(F, dF, xi0, maxIter, tol)
    xi = xi0;
    for k = 1:maxIter
        xi_new = xi - F(xi) / dF(xi);
        if abs(xi_new - xi) < tol
            xi = xi_new;
            return
        end
        xi = xi_new;
    end
end
