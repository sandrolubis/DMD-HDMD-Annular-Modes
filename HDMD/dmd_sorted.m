function [phi_n Lambda_n] = dmd_sorted(phi, Lambda)

[Lambda_n, I] = sort(abs(diag(Lambda)), 'descend');

phi_n = phi(:, I);

Lambda_n = diag(Lambda);
Lambda_n = Lambda_n(I);
Lambda_n = diag(Lambda_n);

end


