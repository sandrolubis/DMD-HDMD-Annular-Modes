f = 'data/u.anom.erainterim.nc';

u = ncread(f, 'u');
lat = ncread(f, 'lat');
lev = ncread(f, 'lev');

lat_index = find(lat >= -80 & lat <= -20);
lev_index = find(lev >= 100 & lev <= 1000);

u_subset = squeeze(u(1, lat_index, lev_index, :));
lat_subset = lat(lat_index);
lev_subset = lev(lev_index);

clear u;
clear lat;
clear lev;
clear lat_index;
clear lev_index;

wgt = sqrt(cos(deg2rad(lat_subset)));

u_subset = u_subset .* wgt;

X = reshape(u_subset, size(u_subset, 1) * size(u_subset, 2), size(u_subset, 3));
X(isnan(X)) = 0.0;

clear wgt;

n = size(X, 1);
m = size(X, 2);
d = 2; %input variable

H = zeros(n * d, m - d + 1);

for i = 1:d
    H(n * (i - 1) + 1:n * i, :) = X(:, i:m - d + i);
end

clear X;

P = 11; %input variable

X = H(:, 1:end - P);
Y = H(:, 1 + P:end);

[U, Sigma, V] = svd(X, 'econ');

r = 100; %input variable

U = U(:, 1:r);
Sigma = Sigma(1:r, 1:r);
V = V(:, 1:r);

A_tilde = U' * Y * V / Sigma;

[W, Lambda] = eig(A_tilde);

Phi = zeros(n * d, r);

for i = 1:r
    Phi(:, i) = (1.0 / Lambda(i, i)) * Y * V / Sigma * W(:, i);
end

[Phi, Lambda] = dmd_sorted(Phi, Lambda);
Phi = my_compress(Phi, d);

Phi_real = real(Phi);
Phi_imag = imag(Phi);

dmd_real = reshape(Phi_real, size(u_subset, 1), size(u_subset, 2), size(Phi_real, 2));
dmd_imag = reshape(Phi_imag, size(u_subset, 1), size(u_subset, 2), size(Phi_imag, 2));

lambda = diag(Lambda);

lambda_real = real(lambda);
lambda_imag = imag(lambda);

clear H;
clear X;
clear Y;
clear U;
clear Sigma;
clear V;
clear A_tilde;
clear W;
clear Lambda;
clear Phi;
clear Phi_real;
clear Phi_imag;
clear lambda;

ndmd = 1:size(dmd_real, 3);
ndmd = ndmd';

g = ['data/hdmd_u_erainterim_d_', num2str(d), '_p_', num2str(P), '_r_', num2str(r), '.nc'];

nccreate(g, 'dmd_real', 'Dimensions', {'lat', size(dmd_real, 1), 'lev', size(dmd_real, 2), 'ndmd', size(dmd_real, 3)}, 'Format', 'netcdf4', 'ChunkSize', [size(dmd_real, 1) size(dmd_real, 2) size(dmd_real, 3)]);
nccreate(g, 'dmd_imag', 'Dimensions', {'lat', size(dmd_imag, 1), 'lev', size(dmd_imag, 2), 'ndmd', size(dmd_imag, 3)}, 'Format', 'netcdf4', 'ChunkSize', [size(dmd_imag, 1) size(dmd_imag, 2) size(dmd_imag, 3)]);
nccreate(g, 'lambda_real', 'Dimensions', {'ndmd', size(lambda_real, 1)}, 'Format', 'netcdf4', 'ChunkSize', size(lambda_real, 1));
nccreate(g, 'lambda_imag', 'Dimensions', {'ndmd', size(lambda_imag, 1)}, 'Format', 'netcdf4', 'ChunkSize', size(lambda_imag, 1));
nccreate(g, 'lat', 'Dimensions', {'lat', size(lat_subset, 1)}, 'Format', 'netcdf4', 'ChunkSize', size(lat_subset, 1));
nccreate(g, 'lev', 'Dimensions', {'lev', size(lev_subset, 1)}, 'Format', 'netcdf4', 'ChunkSize', size(lev_subset, 1));
nccreate(g, 'ndmd', 'Dimensions', {'ndmd', size(ndmd, 1)}, 'Format', 'netcdf4', 'ChunkSize', size(ndmd, 1));

ncwrite(g, 'dmd_real', dmd_real);
ncwrite(g, 'dmd_imag', dmd_imag);
ncwrite(g, 'lambda_real', lambda_real);
ncwrite(g, 'lambda_imag', lambda_imag);
ncwrite(g, 'lat', lat_subset);
ncwrite(g, 'lev', lev_subset);
ncwrite(g, 'ndmd', ndmd);