function [] = scarf3d(fcorr, dir, skip, npts, dh, cl, sigma, acf);



fid = fopen(fcorr, 'rb');

tmp = fread(fid, [npts(1) * npts(2) * npts(3)], 'double');

corr = reshape(tmp, [npts(1) npts(2) npts(3)]);

% acf is along x-direction
if dir == 1
    
    acftitle  = 'ACF along X';
    psdftitle = 'PSDF along X';
    
    % lags vector (only positive lags)
    lag = [0:npts(1)-1] * dh;
    
    kvec = [1:skip:npts(3)];
    jvec = [1:skip:npts(2)];
    
    n = length(kvec) * length(jvec);
    
    var = zeros(npts(1), n);
    
    n = 0;
    
    for k = 1:length(kvec)
        for j = 1:length(jvec)
            n = n + 1;
            var(:, n)  = corr(:, jvec(j), kvec(k));
        end
    end    
    
% acf is along y-direction    
elseif dir == 2
    
    acftitle  = 'ACF along Y';
    psdftitle = 'PSDF along Y';
    
    % lags vector (only positive lags)
    lag = [0:npts(2)-1] * dh;
    
    kvec = [1:skip:npts(3)];
    ivec = [1:skip:npts(1)];
    
    n = length(kvec) * length(ivec);
    
    var = zeros(npts(2), n);
    
    n = 0;
    
    for k = 1:length(kvec)
        for i = 1:length(ivec)
            n = n + 1;
            var(:, n) = corr(ivec(i), :, kvec(k));
        end
    end       
    
% acf is along z-direction    
elseif dir == 3    
    
    acftitle  = 'ACF along Z';
    psdftitle = 'PSDF along Z';    
   
    % lags vector (only positive lags)
    lag = [0:npts(3)-1] * dh;
    
    jvec = [1:skip:npts(2)];
    ivec = [1:skip:npts(1)];
   
    n = length(jvec) * length(ivec);
    
    var = zeros(npts(3), n);
    
    n = 0;
    
    for j = 1:length(jvec)
        for i = 1:length(ivec)
            n = n + 1;
            var(:, n) = corr(ivec(i), jvec(j), :);
        end
    end       
    
end    
    
% compute theoretical ACF
if acf == 'gs'
    fun = sigma^2 * exp(-lag.^2 / cl^2);
elseif acf == 'vk'
    x      = lag / cl(1); 
    fun    = sigma^2 * 2^(1-cl(2)) / gamma(cl(2)) * (x).^cl(2) .* besselk(cl(2), x); 
    fun(1) = sigma^2; 
end

% plot ACF
figure; hold on;

title(acftitle);

for i = 1:size(var, 2)
    plot(lag, var(:, i), 'Color', [0.75 0.75 0.75]);
end

h1 = plot(lag, mean(var, 2), 'b', 'LineWidth', 2);
h2 = plot(lag, fun, 'r', 'LineWidth', 2);
     
legend([h1 h2], {'Average', 'Theoretical'});

xlabel('lag (m)');
ylabel('Biased ACF');

grid on; axis tight;

% compute PSDF
n = size(var, 1);

psdf = zeros(n, size(var, 2));

% compute PSDF
for i = 1:size(var, 2)
    psdf(:, i) = pwr(var(:, i));
end

% compute theoretical PSDF
ks = 1 / dh;
k  = 2 * pi * ks * (0:n-1) / (2*n - 1);
kn = ks * pi / 2;

if acf == 'gs'
    fun = sigma^2 * sqrt(pi) * cl * exp(-cl^2 * k.^2 / 4);
    fun = fun / dh / n;
elseif acf == 'vk'
    x   = (1 + cl(1)^2 * k.^2).^(cl(2) + 0.5); 
    fun = 2 * sqrt(pi) * gamma(cl(2) + 0.5) / gamma(cl(2)) * sigma^2 * cl(1) ./ x;
    fun = fun / dh / n;
end

% plot PSDF
figure; hold on;

title(psdftitle);

for i = 1:size(psdf, 2)
    plot(k*cl(1), psdf(:, i), 'Color', [0.75 0.75 0.75]);
end

h1 = plot(k*cl(1), mean(psdf, 2), 'b', 'LineWidth', 2);
h2 = plot(k*cl(1), fun, 'r', 'LineWidth', 2);

set(gca, 'Xscale', 'log', 'Yscale', 'log');

yl = ylim;

plot([kn*cl(1) kn*cl(1)], yl, '--g');

legend([h1 h2], {'Average', 'Theoretical'});

xlabel('ka');
ylabel('PSDF');

grid on; axis tight;

end


function X = pwr(x);

    % length of one-sided acr
    n = length(x);
    
    % length of two sided: we will have always odd number of points
    npts = 2*n - 1;
    
    y = zeros(npts, 1);
    
    % negative lags
    y(1:n - 1) = x(end:-1:2);
    % positive lags
    y(n:end)   = x;
        
    Y = fft(y, npts);
    
    % one-sided amplitude spectrum: since npts is always odd, we never
    % cover nyquist wavenumber
    X        = abs(Y(1:n)) / npts;
    X(2:end) = 2 * X(2:end);
    
end

    