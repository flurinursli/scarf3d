function [] = scarf3d(prec, fcorr, skip, npts, dh, cl, sigma, acf);

% simple script to compare continuous and discrete (synthetic) ACFs & PSDFs.
% input random field must be isotropic.
% e.g.: scarf3d('single', 'fim_unstruct_whole_3d', 20, [500 500 500], 50, [2000 0.25], 0.05, 'vk');
% e.g.: scarf3d('single', 'fim_unstruct_xslice', 20, [500 500], 50, [2000 0.25], 0.05, 'vk');

fprintf("%s\n", '');
fprintf("%s\n", 'This script compares continuous and discrete ACF&PSD along each direction');

ndim = length(npts);

if ndim == 2; npts = [npts(1) npts(2) 1]; end;

fid = fopen(fcorr, 'rb');

tmp = fread(fid, [prod(npts)], prec);

field = reshape(tmp, [npts(1) npts(2) npts(3)]);

figure;

for dir = 1:ndim

    % acf is along x-direction
    if dir == 1

        direction = 'X';
        acftitle  = 'ACF X-AXIS';
        psdftitle = 'PSDF X-AXIS';

        % lags vector (only positive lags)
        lag = [0:npts(1)-1] * dh;

        kvec = [1:skip:npts(3)];
        jvec = [1:skip:npts(2)];

        n = length(kvec) * length(jvec);

        corr = zeros(npts(1), n);

        n = 0;

        for k = 1:length(kvec)
            for j = 1:length(jvec)
                n = n + 1;
                corr(:, n)  = autocorr(field(:, jvec(j), kvec(k)));
            end
        end

    % acf is along y-direction
    elseif dir == 2

        direction = 'Y';
        acftitle  = 'ACF Y-AXIS';
        psdftitle = 'PSDF Y-AXIS';

        % lags vector (only positive lags)
        lag = [0:npts(2)-1] * dh;

        kvec = [1:skip:npts(3)];
        ivec = [1:skip:npts(1)];

        n = length(kvec) * length(ivec);

        corr = zeros(npts(2), n);

        n = 0;

        for k = 1:length(kvec)
            for i = 1:length(ivec)
                n = n + 1;
                corr(:, n) = autocorr(field(ivec(i), :, kvec(k)));
            end
        end

    % acf is along z-direction
    elseif dir == 3

        direction = 'Z';
        acftitle  = 'ACF Z-AXIS';
        psdftitle = 'PSDF Z-AXIS';

        % lags vector (only positive lags)
        lag = [0:npts(3)-1] * dh;

        jvec = [1:skip:npts(2)];
        ivec = [1:skip:npts(1)];

        n = length(jvec) * length(ivec);

        corr = zeros(npts(3), n);

        n = 0;

        for j = 1:length(jvec)
            for i = 1:length(ivec)
                n = n + 1;
                corr(:, n) = autocorr(field(ivec(i), jvec(j), :));
            end
        end

    end

    % theoretical 1D ACF
    if acf == 'gs'
        fun = sigma^2 * exp(-lag.^2 / cl^2);
    elseif acf == 'vk'
        x      = lag / cl(1);
        fun    = sigma^2 * 2^(1-cl(2)) / gamma(cl(2)) * (x).^cl(2) .* besselk(cl(2), x);
        fun(1) = sigma^2;
    end

    % plot ACF
    handle = subplot(2, ndim, dir);

    hold on;

    title(acftitle);

    for i = 1:size(corr, 2)
        plot(lag, corr(:, i), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.25);
    end

    h1 = plot(lag, mean(corr, 2), 'b', 'LineWidth', 2);
    h2 = plot(lag, fun, 'r', 'LineWidth', 2);

    if (dir == ndim)
        legend([h1 h2], {'Average', 'Theoretical'});
    end

    xlabel('Lag');

    if (dir == 1)
        ylabel('ACF');
    end

    grid on; axis tight;

    % compute discrete PSD
    n = size(corr, 1);

    psdf = zeros(n, size(corr, 2));

    for i = 1:size(corr, 2)
        psdf(:, i) = pwr(corr(:, i));
    end

    % theoretical 1D PSD using discrete resolution ('PWR' double number of points)
    ks = 1 / dh;
    k  = 2 * pi * ks * (0:n-1) / (2*n - 1);
    dk = 2 * pi / (2*n - 1) / dh;

    if acf == 'gs'
        fun = sigma^2 * sqrt(pi) * cl * exp(-cl^2 * k.^2 / 4);
    elseif acf == 'vk'
        x   = (1 + cl(1)^2 * k.^2).^(cl(2) + 0.5);
        fun = 2 * sqrt(pi) * gamma(cl(2) + 0.5) / gamma(cl(2)) * sigma^2 * cl(1) ./ x;
    end

    % continuous std.dev. (assume zero-mean), 1D-equivalent of eqs. 5-6 in Frenje & Juhlin
    csig = sqrt(1/pi * sum(fun(2:end)) * dk);

    % average discrete PSDF
    dPSD = mean(psdf, 2);

    fprintf("%s%s%s %12.5f%s %12.5f%s\n", 'Std.Dev. ', direction, '-direction:', csig, '(cont.)', sqrt(sum(dPSD(2:end))), '(disc.)');

    % normalise continuous PSD to compare with discrete PSD
    fun = fun / dh / n;

    % note that sqrt(sum(fun(2:end))) = csig above

    % now plot PSDFs

    handle = subplot(2, ndim, dir + ndim);

    hold on;

    title(psdftitle);

    for i = 1:size(psdf, 2)
        plot(k*cl(1), psdf(:, i), 'Color', [0.75 0.75 0.75], 'LineWidth', 0.25);
    end

    h1 = plot(k*cl(1), dPSD, 'b', 'LineWidth', 2);
    h2 = plot(k*cl(1), fun, 'r', 'LineWidth', 2);

    set(gca, 'Xscale', 'log', 'Yscale', 'log');

    if (dir == ndim)
        legend([h1 h2], {'Average', 'Theoretical'});
    end

    xlabel('Wavenumber * CL');

    if (dir == 1)
        ylabel('PSDF');
    end

    grid on; axis tight;

end

end

%-------------------------------------------------------------------------------

function x = autocorr(y);
% one-sided autocorrelation of input sequence

    n = length(y);

    z = xcorr(y, 'biased');

    x = z(n:end);

end

%-------------------------------------------------------------------------------

function X = pwr(x);
% PSD from input autocorrelation

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
