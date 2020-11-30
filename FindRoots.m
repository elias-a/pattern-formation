function roots = FindRoots(a1, a2)

    tau = a2 / a1;

    % When searching for critical points, a1 must equal 1 and 
    % a2 must equal tau. To then get the critical points for the normalized 
    % basis vectors, multiply the critical points by the factor
    % 1. / sqrt(imag(tau)) 

    warning('off', 'MATLAB:singularMatrix')
    warning('off', 'MATLAB:nearlySingularMatrix')
    
    cpSearch(a1, a2);


    function cpSearch(a1, a2)
        a1Re = real(a1);
        a1Im = imag(a1);
        a2Re = real(a2);
        a2Im = imag(a2);
        
        % Preallocate space for array of roots 
        roots = zeros(1, 500);
        index = 1;
        
        % Set step sizes 
        stepReal = real(a1 + a2) / 10.;
        stepImag = imag(a1 + a2) / 10.;
        
        % The for loops iterate over a rectangular region.
        % However, the `inLattice` function used below prevents points 
        % from outside the parallelogram from being considered. 
        for j = min([0, a1Re, a2Re]):stepReal:max([(a1Re+a2Re), a1Im, a2Im])
            for k = min([0, a1Im, a2Im]):stepImag:max([(a1Im+a2Im), a1Im, a2Im])
                z = complex(j, k);
                
                % If the theta function equals 0, the gradient of the
                % Green's function will be infinite, so in this case,
                % skip the current loop iteration. 
                if round(theta(z), 5) == 0
                    continue
                end
                
                    lastwarn('') 
                    
                    % Uses `cxroot` program available on MathWorks File Exchange 
                    [root, ssq, ~] = cxroot(@gradG, z, 'FunTol', 1e-10, ...
                                 'XTol', 1e-10, 'MaxIter', 30);
                    
                    root = round(root, 4);
                    
                    % Ignores results from singular matrices   
                    warnMsg = lastwarn;
                    if ~isempty(warnMsg)
                        continue 
                    end
                    
                    % If root has not already been found, 
                    % add root to list.
                    if isempty(find(roots == root, 1))
                        
                        % Checks that the root that has been found
                        % is inside the parallelogram cell. 
                        if inLattice(root, a1, a2)
                            
                            % A large value of the sum of squares
                            % indicates that the point returned by
                            % `cxroot` is not actually a root. 
                            if ssq > 1e-15
                                continue
                            end
                            
                            roots(index) = root;
                            index = index + 1;
                        end             
                    end
            end
        end
        
        roots = roots(1:index-1);
    end
    
    % Checks if point `z` is in the parallelogram cell
    function in = inLattice(z, a1, a2)
    
        x1 = [0, 0];
        x2 = [real(a1), imag(a1)];
        x3 = [real(a1+a2), imag(a1+a2)];
        x4 = [real(a2), imag(a2)];
    
        x = [x1(1), x2(1), x3(1), x4(1)];
        y = [x1(2), x2(2), x3(2), x4(2)];
    
        in = inpolygon(real(z), imag(z), x, y); 
    
        % Function checks if `point` is on the line between
        % points `p1` and `p2`. 
        function on = onLine(p1, p2, point)
            on = false;
        
            % Calculate normal along the line connecting p1 and p2
            p12 = p2 - p1;
            l12 = sqrt(p12 * p12');
            n12 = p12 / l12;
        
            % Calculate normal along the line connecting p1 and point
            p1p = point - p1;
            l1p = sqrt(p1p * p1p');
            n1p = p1p / l1p;
        
            % If the normals are parallel, then the point is on the line
            if round(dot(n12, n1p), 3) == 1
                on = true;
            end
        end
    
        % Note that points on the right and top edges of the cell
        % are considered to be in other cells. 

        % Checks whether `z` is on the right edge of the cell.
        if onLine(x2, x3, [real(z), imag(z)])
            in = false;
        end
    
        % Checks whether `z` is on the top edge of the cell. 
        if onLine(x4, x3, [real(z), imag(z)])
            in = false;
        end
    end
    
    % Gradient of the Green's function 
    function dGdz = gradG(z)    
        dGdz = (thetaPrime(z) / theta(z)) ...
            + (2 * pi * 1i * imag(z) / imag(tau));
    end
    
    % Theta function, based on Lin and Wang 2010 
    function t = theta(z)
        q = exp(1i * pi * tau);
        t = 0;
    
        for n = 0:30
            t = t + ((-1)^n * q^((n+0.5)^2) * sin((2*n+1)*pi*z));
        end
    
        t = t * 2;
    end
    
    % Derivative of the theta function 
    function tp = thetaPrime(z)
        q = exp(1i * pi * tau);
        tp = 0;
    
        for n = 0:30
            tp = tp + ((-1)^n * q^((n+0.5)^2) * (2*n+1)*pi * cos((2*n+1)*pi*z));
        end
    
        tp = tp * 2;
    end
end