function g = GreenTheta(z, a1, a2)
    tau = a2 / a1;
    
    % Calculate the constant in the Green's function
    function C = GreenConstant()

        C = Square();

        function Ctau = Square() 
            Gfunc = @(x, y) IntegrateG(x, y);
    
            % The constant is the negative of the integral over the cell
            Ctau = -integral2(Gfunc, 0, 1, 0, 1);
            
            function G = IntegrateG(w1, w2)
        
                % Change coordinates to integrate over [0,1]x[0,1] 
    
                % Transformation matrix 
                A = [real(a1) real(a2); imag(a1) imag(a2)];
    
                % integral2 sends in a matrix of values, 
                % this code multiplies each value in the matrix by A 
                x = zeros(size(w1));
                y = zeros(size(w2));
                for i = 1 : size(w1, 1)
                    for j = 1 : size(w1, 2)
                        Z = A * [w1(i, j); w2(i, j)];
                        x(i, j) = Z(1);
                        y(i, j) = Z(2);
                    end
                end
         
                G = (-1 / (2 * pi)) * log(abs(theta(complex(x, y)))) ...
                    + (y.^2 / (2 * imag(tau)));
            end  
        end
    end

    % Theta function, based on Lin and Wang 2010 
    function t = theta(z)
        q = exp(1i * pi * tau);
        t = 0;
    
        for n = 0:50
            t = t + ((-1)^n * q^((n+0.5)^2) * sin((2*n+1)*pi*z));
        end
    
        t = t * 2;
    end
    
    g = (-1 / (2 * pi)) * log(abs(theta(z))) ...
        + (imag(z) ^ 2 / (2 * imag(tau))) + GreenConstant();
end