% A very fast function to estimate the Gamma function with complex variable
% @see https://doi.org/10.1007/s00013-010-0146-9 Equation (4.1) and Euler's reflection formula
function y = gamma_complex(z)
    % seperate matrix z to real(z) >= 1 and real z < 1 
    z_larger_1  = z(real(z) >= 1);
    z_smaller_1 = z(real(z) < 1);
    
    % When real(z) > 1
    y_with_z_larger_1 = gamma_complex_0(z_larger_1);
    
    % When real(z) < 1, use Euler's reflection formula
    y_with_z_smaller_1 = pi ./ sin(pi .* z_smaller_1) ...
            ./ gamma_complex_0(1 - z_smaller_1);
    
    % Combine two part and return
    y = zeros(size(z));
    y(real(z) >= 1) = y_with_z_larger_1;
    y(real(z) < 1)  = y_with_z_smaller_1;
end

% The Gamma function for real(z) > 1
% @see https://doi.org/10.1007/s00013-010-0146-9 Equation (4.1)
function y = gamma_complex_0(z)
    y = (2*pi ./z ).^0.5 .* ...
            (1/exp(1) .* (z + 1./ (12.*z - 1/10./z) ) ) .^ z;
end

% for test
function gamma_test()
    z_r = linspace(-6, 4, 601);
    z_i = linspace(-5, 5, 601);
    [Z_real, Z_imag] = meshgrid(z_r, z_i);
    
    Z = complex(Z_real, Z_imag);
    W_gamma = gamma_complex( Z);
    
    W_ang = angle(W_gamma);
    W_abs = abs(W_gamma);
    
    figure('Position', [40 40, 500, 500]);
    imagesc([-6, 4], [-5, 5], W_ang);
    colormap(hsv(256));
    
    figure('Position', [40 40, 500, 500]);
    imagesc([-6, 4], [-5, 5], log(W_abs), [-3, 5]);
    colormap(parula(256));
end
