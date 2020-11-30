function minEnergy = MinEnergy(K, a1, a2, point, G)

    % Local term
    function P = Pfunc(K, pho) 
        P = K * 2 * pi * pho;
    end

    % Nonlocal term
    function I = Ifunc(K, pho, a1, a2, G)
        % Consider only lattices of size 1
        Lambda = 1.;
    
        % Value of the harmonic function at 0
        function h = hZero()
        
            infSum = 0;
            for n=1:1000
                infSum = infSum + log(abs((1 - exp(2*pi*1i*n*(a2/a1)))^2));
            end

            h = (-1/(2*pi)) * (log(abs(sqrt(imag(a2/a1)) * exp(2*pi*1i*(a2/a1)/12))) + infSum);
        end
        
        if K > 1
            Gterm = (K - 1) * pi^2 * pho^4 * G(point, a1, a2);
            
        else
            Gterm = 0.;
        end
        
        I = Gterm + K * ((pi^2 * pho^4 * hZero()) + ...
            (pi * pho^4 / 8) - (pi * pho^4 * ...
            log(2*pi*pho/sqrt(abs(Lambda))) / 2)) + ...
            (K^2 * pi^2 * pho^6 / (4 * abs(Lambda))); 
    end

    % omega must be in (0, 1)
    omega = 0.01;  
    pho = sqrt(omega / (K * pi));
    
    % Computes the K value 
    minEnergy = (Pfunc(K, pho) ^ (2/3)) * ...
                (Ifunc(K, pho, a1, a2, G) ^ (1/3));
end