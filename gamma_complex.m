% A very fast function to estimate the Gamma function with complex variable
% @see https://doi.org/10.1007/s00013-010-0146-9 Equation (4.1) and Euler's reflection formula
% Seperate the Gamma function in real(z) >= 1, real(z) < 1 and z ~ 1.0 to calculate.
% Warning: the result may no continous in connection part.
function y = gamma_complex(z)
    % seperate matrix z to real(z) >= 1, real(z) < 1 and z ~ 1.0 
    z_larger_1  = z(real(z) >= 1);
    z_smaller_1 = z(real(z) < 1 & abs(z-1) >= 0.35);
    z_smaller_1_near_1 = z(real(z) < 1 & abs(z-1) < 0.35);
    
    % When real(z) > 1
    y_with_z_larger_1 = gamma_complex_0(z_larger_1);
    
    % When real(z) < 1, use Euler's reflection formula
    y_with_z_smaller_1 = pi ./ sin(pi .* z_smaller_1) ...
            ./ gamma_complex_0(1 - z_smaller_1);
    
    % When z ~ 1, use Taylor Expansion at 1.0
    y_with_z_smaller_1_near_1 = gamma_complex_near_1(z_smaller_1_near_1);
    
    % Combine three part and return
    y = zeros(size(z));
    y(real(z) >= 1) = y_with_z_larger_1;
    y(real(z) < 1 & abs(z-1) >= 0.35)  = y_with_z_smaller_1;
    y(real(z) < 1 & abs(z-1) < 0.35)  = y_with_z_smaller_1_near_1;
end

% The Gamma function for real(z) > 1
% @see https://doi.org/10.1007/s00013-010-0146-9 Equation (4.1)
function y = gamma_complex_0(z)
    y = (2*pi ./z ).^0.5 .* ...
            (1/exp(1) .* (z + 1./ (12.*z - 1/10./z) ) ) .^ z;
end

% The Gamma function for abs(z-1) ~ 0
% The Taylor expansion of Gamma(x=1)
% const_2 = (6*eulergamma^2+pi^2)/12 = 0.9890 ...
function y = gamma_complex_near_1(z)
    the_eulergamma = ...
        0.577215664901532860606512090082402431042159335939923598805;
    const_2 = ...
        0.989055995327972555395395651500634707939183520728214090443;
    
    y = 1 - the_eulergamma .* (z-1) + const_2 .* (z-1).^2;
end

% for test
function gamma_test()
    z_r = linspace(-6, 4, 601);
    z_i = linspace(-5, 5, 601);
    [Z_real, Z_imag] = meshgrid(z_r, z_i);
    
    Z = complex(Z_real, Z_imag);
    W_gamma = gamma_complex( Z);
    
    W_arg = angle(W_gamma);
    W_abs = abs(W_gamma);
    
    figure('Position', [40 40, 500, 500]);
    imagesc([-6, 4], [-5, 5], W_arg);
    colormap(hsv(256));
    title('Arg \Gamma(z)');
    xlabel('Re z');
    ylabel('Im z');
    
    figure('Position', [40 40, 500, 500]);
    imagesc([-6, 4], [-5, 5], log(W_abs), [-3, 5]);
    colormap(parula(256));
    title('Abs \Gamma(z)');
    xlabel('Re z');
    ylabel('Im z');
end
