function MinLattice(K, filePath, G, cases) 

    % The parameter `cases` can assume one of three values
    % (all referring to the two disc configuration): 
    %    1) 'half' - place second disc at the minimal half period
    %    2) 'five' - consider only lattices with 5 critical points 
    %    3) 'all'  - consider all lattices and place the second disc
    %                at the minimal critical point

    
    file = fopen(filePath, 'w');
    
    % Iterate over the fundamental domain in the tau plane
    for j = -0.495:0.005:0.5
        for k = 0.5:0.005:1      
            tau = complex(j, k);
            
            % Skip points outside the fundamental domain
            if abs(tau) < 1
                continue
            end
            
            % Normalize the basis vectors
            a1 = 1 / sqrt(imag(tau));
            a2 = tau / sqrt(imag(tau));
            
            switch K
                case 1
                    % The `point` parameter of MinEnergy is not used for K=1.
                    % 0 is passed for simplicity but any value can be used. 
                    energy = MinEnergy(K, a1, a2, 0, G);
                    
                case 2
                    switch cases
                        case 'half'  
                            % Determine the energy for each of the half periods,
                            % and choose the minimum. 
                            energy1 = MinEnergy(K, a1, a2, a1/2, G);
                            energy2 = MinEnergy(K, a1, a2, a2/2, G);
                            energy12 = MinEnergy(K, a1, a2, (a1+a2)/2, G);
                            
                            energy = min([energy1, energy2, energy12]);
                            
                        case 'five'
                            % When using the FindRoots routine, we require that
                            % the basis vectors be 1 and tau.
                            w1 = 1.;
                            w2 = tau;
                    
                            % Search for roots and normalize them 
                            roots = FindRoots(w1, w2) / sqrt(imag(tau));
                            
                            %  Skip lattices with 3 critical points 
                            if length(roots) < 5
                                continue
                            end
                            
                            gList = zeros(1, length(roots));
                            index = 1;
                            for root = roots
                                gList(index) = Green(root, a1, a2);
                                index = index + 1;
                            end
                    
                            % Determine the critical point that minimizes
                            % the Green's function (there should be two
                            % minima, so take the smaller one). 
                            point = min(roots(find(gList == min(gList))));
                    
                            energy = MinEnergy(K, a1, a2, point, G);
                        
                        case 'all'
                            % When using the FindRoots routine, we require that
                            % the basis vectors be 1 and tau.
                            w1 = 1.;
                            w2 = tau;
                    
                            % Search for roots and normalize them 
                            roots = FindRoots(w1, w2) / sqrt(imag(tau));
                            
                            gList = zeros(1, length(roots));
                            index = 1;
                            for root = roots
                                gList(index) = Green(root, a1, a2);
                                index = index + 1;
                            end
                    
                            % Determine the critical point that minimizes
                            % the Green's function (there should be two
                            % minima, so take the smaller one). 
                            point = min(roots(find(gList == min(gList))));
                    
                            energy = MinEnergy(K, a1, a2, point, G);
                    end
            end
            
            [j, k, energy] 
            
            fprintf(file, '%.3f %.3f %.10f\n', j, k, energy);
        end
    end

    fclose(file);
end