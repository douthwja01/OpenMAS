    for j=1:n
        pj = p_all(:,j);                                                   % Position of the second agents from the global set 
        if j~=i && A(i,j)~=0                                               % If the reference agent ID and eval agent ID do not match or to not neighbour
            %% POSITION BASED FORMATION CONTROL LAW 
            %pi_star=p_star_all(:,i);
            %pj_star=p_star_all(:,j);
            %fi=fi+(pj-pi-(pj_star-pi_star));
            
            %% DISTANCE BASED FORMATION CONTROL LAW
            ell_ij = dis_star(i,j);                                        % Get the desired distance between agents i & j
            fi = fi + (norm(pi - pj)^2 - ell_ij^2)*(pj - pi);              % Get the formation control feedback for that agent (fi)              
        end
    end