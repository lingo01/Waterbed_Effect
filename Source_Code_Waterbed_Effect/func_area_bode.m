function [area, area_array] = func_area_bode(wout, mag)
    area_array = zeros(size(wout));

    % x_mag = real(mag);
    % y_mag = imag(mag);
    % ratio_mag = y_mag./x_mag;
    % 
    % log10_abs_mag = log10(abs(x_mag)) + 0.5 * log10(1+ratio_mag.^2);

    log10_abs_mag = mag;

    for ii = 2:length(wout)
        area_array(ii) = area_array(ii-1) + log10_abs_mag(ii) * (wout(ii)-wout(ii-1));
    end

    area = area_array(end);
end