function g = Green(z, a1, a2)
    
    function H = infProduct(z) 
        H = 1;
        for n=1:1000
            H = H * (1 - exp(2*pi*1i*(n*(a2/a1)+(z/a1)))) ...
                * (1 - exp(2*pi*1i*(n*(a2/a1)-(z/a1))));
        end
    
        H = log(abs(H));
    end

    Lambda = det([real(a1) real(a2) ; imag(a1) imag(a2)]);
    
    infProd = infProduct(z);

    g = (abs(z)^2 / (abs(Lambda) * 4)) - (1 / (2 * pi)) * ...
        log(abs(exp(2*pi*1i*((z^2 * conj(a1) / (4i * ...
        abs(Lambda) * a1)) - (z / (2 * a1)) + (a2 / (12 * a1)))) * ...
        (1 - exp(2*pi*1i*(mod(real(z), 1) + 1i*imag(z))/a1)))) - ...
        (1 / (2 * pi)) * infProd;
end