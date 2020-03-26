function [] = scarf2d(prec, fcorr, dir, skip, npts, dh, cl, sigma, acf);
%function [] = scarf2d(corr, dir, skip, npts, dh, cl, sigma, acf);


fid = fopen(fcorr, 'rb');

tmp = fread(fid, [npts(1) * npts(2)], prec);

corr = reshape(tmp, [npts(1) npts(2)]);

% acf is along i-direction
if dir == 1
    
    acftitle  = 'ACF along i';
    psdftitle = 'PSDF along i';
    
    % lags vector (only positive lags)
    lag = [0:npts(1)-1] * dh;
    
    jvec = [1:skip:npts(2)];
    
    n = length(jvec);
    
    var = zeros(npts(1), n);
    
    n = 0;
    
    for j = 1:length(jvec)
       n = n + 1;
       var(:, n)  = autocorr(corr(:, jvec(j)));
    end  
    
% acf is along j-direction    
elseif dir == 2
    
    acftitle  = 'ACF along j';
    psdftitle = 'PSDF along j';
    
    % lags vector (only positive lags)
    lag = [0:npts(2)-1] * dh;
    
    ivec = [1:skip:npts(1)];
    
    n = length(ivec);
    
    var = zeros(npts(2), n);
    
    n = 0;
    
    for i = 1:length(ivec)
        n = n + 1;
        var(:, n) = autocorr(corr(ivec(i), :));
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

% xx = pwr(fun);
% xx(1), xx(2)

% plot ACF
figure;

handle = subplot(1, 2, 1);

hold on;

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
%kn = pi / dh;
dk = 2 * pi / (2*n - 1) / dh;

if acf == 'gs'
    fun = sigma^2 * sqrt(pi) * cl * exp(-cl^2 * k.^2 / 4);
    fun = fun / dh / n;
elseif acf == 'vk'
    x   = (1 + cl(1)^2 * k.^2).^(cl(2) + 0.5); 
    fun = 2 * sqrt(pi) * gamma(cl(2) + 0.5) / gamma(cl(2)) * sigma^2 * cl(1) ./ x;
    fun = fun / dh / n;
end

% fun(1), fun(2)

xx = ipwr(fun);
h3 = plot(lag, xx, '--k', 'LineWidth', 2);

legend([h1 h2 h3], {'Average', 'Theoretical', 'From PSDF'});

%fun = pwr(xx);

['Std.Dev from continuous PSD: ', num2str(sqrt( 1/ pi * sum(fun(2:end) * dh * n) * dk))]

sig = mean(psdf, 2);

['Average Std.Dev. from discrete PSD: ', num2str(sqrt( sum(sig(2:end))))]


% plot PSDF
%figure; hold on;

handle = subplot(1, 2, 2);

hold on;

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


function x = autocorr(y);

    n = length(y);

    z = xcorr(y, 'biased');
    %z = xcorr(y, 'none');
    
    x = z(n:end);
    
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

function x = ipwr(X);
        
    n = length(X);
    
    npts = 2*n - 1;

    Y = zeros(npts, 1);
    
    Y(1) = X(1);
    
    Y(2:n) = X(2:n) / 2;
    
    y = ifft(Y, 'symmetric');
    
    x = y(1:n) * npts;
    
end
    
