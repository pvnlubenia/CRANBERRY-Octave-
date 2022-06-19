% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                             %
%    network_numbers                                                          %
%                                                                             %
%                                                                             %
% OUTPUT: Returns the values of the following network numbers of a chemical   %
%            reaction network:                                                %
%               - Species (m)                                                 % 
%               - Complexes (n)                                               % 
%               - Reactant complexes (n_r)                                    % 
%               - Reversible reactions (r_rev)                                %
%               - Irreversible reactions (r_irrev)                            %
%               - Reactions (r) = r_irrev + 2 r_rev                           %
%               - Linkage classes (l)                                         %
%               - Strong linkage classes (sl)                                 %
%               - Terminal linkage classes (t)                                %
%               - Rank (s)                                                    %
%               - Reactant rank (q)                                           %
%               - Deficiency (delta) = n - l - s                              %
%               - Reactant deficiency (delta_p) = n_r - q                     %
%         The output variable 'model' allows the user to view the complete    %
%            network with all the species listed in the 'species' field of    %
%            the structure 'model'.                                           %
% INPUT: model: a structure, representing the CRN, with the following fields  %
%           (the kinetics of the network is not needed):                      %
%           - id: name of the model                                           %
%           - species: a list of all species in the network; this is left     %
%                blank since incorporated into the function is a step which   %
%                compiles all species used in the model                       %
%           - reaction: a list of all reactions in the network, each with the %
%                following subfields:                                         %
%                   - id: a string representing the reaction                  %
%                   - reactant: has the following further subfields:          %
%                        - species: a list of strings representing the        %
%                             species in the reactant complex                 %
%                        - stoichiometry: a list of numbers representing the  %
%                             stoichiometric coefficient of each species in   %
%                             the reactant complex (listed in the same order  %
%                             of the species)                                 %
%                   - product: has the following further subfields:           %
%                        - species: a list of strings representing the        %
%                             species in the product complex                  %
%                        - stoichiometry: a list of numbers representing the  %
%                             stoichiometric coefficient of each species in   %
%                             the product complex (listed in the same order   %
%                             of the species)                                 %
%                   - reversible: has the value true or false indicating if   %
%                        the reaction is reversible or not, respectively      %
%                                                                             %
% References                                                                  %
%    [1] Arceo C, Jose E, Lao A, Mendoza E (2017) Reactant subspaces and      %
%           kinetics of chemical reaction networks. J Math Chem 56:395–422.   %
%           http://doi.org/10.1007/s10910-017-0809-x                          %
%    [2] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical reaction %
%           network theory. Bioinform 25(21):2853–2854.                       %
%           https://doi.org/10.1093/bioinformatics/btp513                     %
%                                                                             %
% Created: 13 June 2021                                                       %
% Last Modified: 19 June 2022                                                 %
%                                                                             %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function model = network_numbers(model)
    
    %
    % Add to 'model.species' all species indicated in the reactions
    %
    
    % Get all species from reactants
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).reactant)
            model.species{end+1} = model.reaction(i).reactant(j).species;
        end
    end
    
    % Get species from products
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).product)
            model.species{end+1} = model.reaction(i).product(j).species;
        end
    end
    
    % Get only unique species
    model.species = unique(model.species);
    
    
    
    %
    % Species
    %
    
    % Count the number of species
    m = numel(model.species);
    
    
    
    %
    % Complexes
    %
    
    % Initialize the matrix of reactant complexes
    reactant_complex = [ ];
    
    % Initialize the matrix of product complexes
    product_complex = [ ];
    
    % Initialize the stoichiometric matrix
    N = [ ];
    
    % For each reaction in the model
    for i = 1:numel(model.reaction)
      
        % Initialize the vector for the reaction's reactant complex
        reactant_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the reactant complex
        for j = 1:numel(model.reaction(i).reactant)
            reactant_complex(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
        end
        
        % Initialize the vector for the reaction's product complex
        product_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the product complex
        for j = 1:numel(model.reaction(i).product)
            product_complex(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
        end
        
        % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
        N(:, end+1) = product_complex(:, end) - reactant_complex(:, end);
        
        % If the reaction is reversible
        if model.reaction(i).reversible
          
            % Insert a new vector for the reactant complex: make it same as the product complex
            reactant_complex(:, end+1) = product_complex(:, end);
            
            % Insert a new vector for the product complex: make it the same as the reactant complex
            product_complex(:, end+1) = reactant_complex(:, end-1);
            
            % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
            N(:, end+1) = -N(:, end);
        end
    end
    
    % Get just the unique complexes
    % ind2(i) is the index in Y of the reactant complex in reaction i
    [Y, ~, ind2] = unique([reactant_complex, product_complex]', 'rows');
    
    % Construct the matrix of complexes
    Y = Y';
    
    % Count the number of complexes
    n = size(Y, 2);
    
    
    
    %
    % Reactant complexes
    %
    
    % Get just the unique reactant complexes
    reactant_complexes_unique = unique(reactant_complex', 'rows')';
    
    % Count the number of unique reactant complexes
    n_r = size(reactant_complexes_unique, 2);
    
    
    
    %
    % Reversible, irreversible, and total reactions
    %
    
    % Count the number of reversible, irreversible, and total reactions
    r_rev = 0;
    r_irrev = 0;
    for i = 1:numel(model.reaction)
        if (model.reaction(i).reversible == 1)
            r_rev = r_rev + 1;
         else
            r_irrev = r_irrev + 1;
        end
    end
    r = r_irrev + 2*r_rev;
    
    
    
    %
    % Linkage classes
    %
    
    % Initialize an undirected graph g
    g = init_graph();

    % Go through each column of Y (a complex)
    for i = 1:n
        
        % For the zero complex
        if numel(find(Y(:, i))) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(find(Y(:, i)))
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [model.species{find(Y(:, i))(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [complex, '+', model.species{find(Y(:, i))(j)}];
                    else
                        complex = [complex, '+', num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                end
            end 
        end
        
        % Add this complex in the list of vertices of g
        g = add_vertex(g, complex);
    end
    
    % Add edges to g: Ci -> Cj forms an edge
    % ~ suppresses the original output
    for i = 1:r
        g = add_edge(g, g.vertices{[~, loc] = ismember(reactant_complex(:, i)', Y', 'rows')}, g.vertices{[~, loc] = ismember(product_complex(:, i)', Y', 'rows')});
    end
    
    % Initialize the vector which will indicate in which linkage class number a vertex (i.e., complex) belongs to
    linkage_class = zeros(numel(g.vertices), 1);
    
    % Initialize the linkage class number tracker
    linkage_class_num = 0;
    
    % Go to each vertex
    for i = 1:numel(g.vertices)
        
        % Pay attention only to a vertex which has no linkage class number yet
        if linkage_class(i) == 0
            
            % This vertex goes to the next linkage class number
            linkage_class_num += 1;
            
            % Assign the linkage class number to the vertex
            linkage_class(i) = linkage_class_num;
            
            % Take note of the vertex that needs to be checked for edges
            to_check = [i];
            
            % Continue assigning a linkage class number to vertices that get into the check list
            while ~isempty(to_check)
                
                % Get the vertex in the check list
                v1 = to_check(end);
                
                % Remove the vertex from the check list (since we now know which vertex to focus on)
                to_check(end) = [ ];
                
                % Check to which vertices the vertex is connected to
                for j = 1:numel(g.edges{v1})
                    
                    % Take note of the vertex it is connected to
                    v2 = g.edges{v1}(j).vertex;
                    
                    % Pay attention to this vertex if it has no linkage class number yet
                    if linkage_class(v2) == 0
                        
                        % Assign this vertex with the same linkage class number as the vertex it is connected to
                        linkage_class(v2) = linkage_class_num;
                        
                        % Add this vertex to our check list: in the next round, we'll check to which other vertices it is connected to
                        to_check(end+1) = v2;
                    end
                end
            end
        end
    end
    
    % Count the number of linkage classes
    l = max(linkage_class);
    
    
    
    %
    % Strong linkage classes
    %
    
    % Initialize a directed graph g
    g = init_graph();

    % Go through each column of Y (a complex)
    for i = 1:n
        
        % For the zero complex
        if numel(find(Y(:, i))) == 0
            complex = '0';
        
        % Otherwise
        else
            
            % Check which species appear in the complex
            for j = 1:numel(find(Y(:, i)))
                
                % For the first species
                if j == 1
                    
                    % Don't show the stoichiometry if it's 1
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [model.species{find(Y(:, i))(j)}];
                    
                    % Otherwise, include it
                    else
                        complex = [num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                
                % We need the + sign for succeeding species
                else
                    if Y(:, i)(find(Y(:, i))(j)) == 1
                        complex = [complex, '+', model.species{find(Y(:, i))(j)}];
                    else
                        complex = [complex, '+', num2str(Y(:, i)(find(Y(:, i))(j))), model.species{find(Y(:, i))(j)}];
                    end
                end
            end 
        end
        
        % Add this complex in the list of vertices of G
        g = add_vertex(g, complex);
    end
    
    % Add a directed edge to g: Ci -> Cj forms an edge
    for i = 1:r
        g = add_path(g, g.vertices{[~, loc] = ismember(reactant_complex(:, i)', Y', 'rows')}, g.vertices{[~, loc] = ismember(product_complex(:, i)', Y', 'rows')});
    end
    
    % Define function which visits each complex (i.e., vertex) v and the other vertices connected to it
    % This has to be placed here (this was taken from [2])
    function visit(v)
        
        % Set the discovery time of the complex as the current time
        discovery_time(v) = time;
        
        % Set the discovery time of the strong linkage class as the current time
        slc_discovery_time(v) = time;
        
        % Move the time forward for the next complex
        time = time + 1;
        
        % Add the complex in the list of complexes in the same strong linkage class
        stack(end+1) = v;
        
        % Note that the complex is already listed in the strong linkage class
        on_stack(v) = true;
        
        % Go through each edge connected to the vertex (i.e., complex)
        for j = 1:numel(g.edges{v})
            
            % Take the vertex connected to it
            v2 = g.edges{v}(j).vertex;
            
            % If the vertex is not yet visited
            if discovery_time(v2) == 0
                
                % Apply the visit function to this vertex
                visit(v2);
                
                % Set the discovery time of the strong linkage class
                % slc_discovery_time(v2) < slc_discovery_time(v) iff a vertex in the stack before v is reachable from v2 (and so they are all in the same strong linkage class)
                slc_discovery_time(v) = min(slc_discovery_time(v), slc_discovery_time(v2));
            
            % If v2 was visited before v
            elseif on_stack(v2)
                
                % So v and v2 are in the same component, and they must have the same slc_discovery_time
                slc_discovery_time(v) = min(slc_discovery_time(v), slc_discovery_time(v2));
            end
        end
        
        % If v is the first visited node of its strong linkage class, all the other vertices of the strong linkage class follow it on the stack
        if slc_discovery_time(v) == discovery_time(v) 
            while true
                v2 = stack(end);
                stack(end) = [ ];
                on_stack(v2) = false;
                strong_linkage_class(v2) = slc_num;
                if v2 == v
                    break
                end
            end
            
            % Add 1 to the strong linkage class number
            slc_num = slc_num + 1;
        end
    end
    
    % Initialize the discovery time of the vertices
    discovery_time = zeros(numel(g.vertices), 1);
    
    % Initialize the discovery time of the [first visited complex of the] strong linkage class of the vertices
    slc_discovery_time = zeros(numel(g.vertices), 1);
    
    % Start the timer at 1
    time = 1;
    
    % Initialize the list of complexes in the strong linkage class
    stack = [ ];
    
    % Initialize that no vertex is on a strong linkage class
    on_stack = false(numel(g.vertices), 1);
    
    % Initialize vector of strong linkage class numbers
    strong_linkage_class = zeros(numel(g.vertices), 1);
    
    % Start numbering the strong linkage class at 1
    slc_num = 1;
    
    % Go through each complex (i.e., vertex)
    for i = 1:numel(g.vertices)
        
        % Use the visit function if it has no strong linkage classs number
        if strong_linkage_class(i) == 0
            visit(i);
        end
    end
    
    % Count the number of strong linkage classes
    sl = max(strong_linkage_class);
    
    
    
    %
    % Terminal linkage classes
    %
    
    % Initialize a matrix (complexes x complexes) for the reacts_to relation
    % This is for testing reversibility of the network
    reacts_to = false(n, n);
    
    % Fill out the entries of the matrix
    for i = 1:r
        
        % reacts_to(i, j) = true iff there is a reaction r: y_i -> y_j
        reacts_to(ind2(i), ind2(i + r)) = true;
    end
    
    % Count the number of terminal strong linkage classes
    % Initialize the count
    t = 0;
    is_nontrivial_terminal_slc = false(sl, 1);
    is_terminal_complex = false(n, 1);
    for i = 1:sl
        
        % Locate the indexes in Y of the complexes present in strong-linkage class i
        complexes_i = find(strong_linkage_class == i);
        
        % Locate the indexes in Y of the complexes which in some reactions are products of complexes_i
        products_of_complexes_i = find(any(reacts_to(complexes_i, :), 1));
        
        % Products_of_complexes_i is a subset of complexes_i, so the strong-linkage class i is terminal
        if all(ismember(products_of_complexes_i, complexes_i))
            t = t + 1;
            is_terminal_complex(complexes_i) = true;
            if numel(complexes_i) > 1
                is_nontrivial_terminal_slc(i) = true;
            end
        end
    end
    
    
    
    %
    % Rank
    %
    
    % Get the rank of the reaction network
    % S = Im N
    % dim S = dim (Im N) = rank(N)
    % Note: We talk of "dimension of a linear transformation" and "rank of a matrix"
    s = rank(N);
    
    
    
    %
    % Reactant rank
    %
    
    % Initialize matrix (complexes x total reactions) for the reacts_in relation
    % This is the incidence matrix I_a
    reacts_in = zeros(n, r);
    
    % Fill out the entries of the matrix
    for i = 1:r
        
        % reacts_in(i, r) = -1 and reacts_in(j, r) = 1) iff there is a reaction r: y_i -> y_j
        reacts_in(ind2(i), i) = -1;
        reacts_in(ind2(i+r), i) = 1;
    end
    
    % Construct the incidence matrix
    % We can decompose this into I_a = I_a^+ - I_a^-
    I_a = reacts_in;
    
    % Construct I_a^-: for reactant complexes only
    I_a_minus = I_a;
    I_a_minus(I_a_minus > 0) = 0;
    I_a_minus(I_a_minus < 0) = 1;
    
    % Construct N^-: "reactant subspace" matrix
    N_minus = Y*I_a_minus;
    
    % Get the reactant rank
    % R = Im N^-
    % dim R = dim (Im N^-) = rank(N^-)
    q = rank(N_minus);
    
    
    %
    % Deficiency
    %
    
    % Compute the deficiency of the reaction network
    delta = n - l - s;
    
    
    
    %
    % Reactant deficiency
    %
    
    % Compute the reactant deficiency
    delta_p = n_r - q;
    
    
    
    %
    % Display results
    %
    
    % Print header
    % Use 'fprintf' instead of 'disp' to interpret '\n' as 'newline'
    fprintf(['NETWORK NUMBERS - ', model.id '\n\n']);
    
    % Save each result as a string
    mR = ['Species (m): ', num2str(m)];
    nR = ['Complexes (n): ', num2str(n)];
    n_rR = ['Reactant complexes (n_r): ', num2str(n_r)];
    r_revR = ['Reversible reactions (r_rev): ', num2str(r_rev)];
    r_irrevR = ['Irreversible reactions (r_irrev): ', num2str(r_irrev)];
    rR = ['Reactions (r): ', num2str(r)];
    lR = ['Linkage classes (l): ', num2str(l)];
    slR = ['Strong linkage classes (sl): ', num2str(sl)];
    tR = ['Terminal linkage classes (t): ', num2str(t)];
    sR = ['Rank (s): ', num2str(s)];
    qR = ['Reactant rank (q): ', num2str(q)];
    deltaR = ['Deficiency (delta): ', num2str(delta)];
    delta_pR = ['Reactant deficiency (delta_p): ', num2str(delta_p)];
    
    % Combine all results in a vector for display
    disp([mR;
          nR;
          n_rR;
          r_revR;
          r_irrevR;
          rR;
          lR;
          slR;
          tR;
          sR;
          qR;
          deltaR;
          delta_pR]);
    fprintf('\n\n')

end










% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% %                                                   % %
% % The following are functions used in the algorithm % %
% %                                                   % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Function 1 of 4: init_graph                                             %
%                                                                         %
%    - Purpose: To initialize an undirected graph                         %
%    - Input: none                                                        %
%    - Output                                                             %
%         - g: empty structure with subfields 'vertices' and 'edges'      %
%    - Used in network_numbers (linkage classes, strong linkage classes)  %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = init_graph()
    
    % Initialize a structure with empty fields
    g = struct();
    
    % Initialize 'vertices' field
    g.vertices = cell();
    
    % Initialize 'edges' field
    g.edges = cell();

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Function 2 of 4: add_vertex                                             %
%                                                                         %
%    - Purpose: To add a vertex to a graph g                              %
%    - Inputs                                                             %
%         - g: graph structure                                            %
%         - v: name of vertex to be added                                 %
%    - Output                                                             %
%         - g: structure with vertex added                                %
%    - Used in network_numbers (linkage classes, strong linkage classes)  %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = add_vertex(g, v)
    
    % Determine the index of 'v', if it already exists, in the 'vertices' field
    location = find(strcmp(v, g.vertices));
    
    % Case 1: 'location' is empty
    if isempty(location)
        
        % Add vertex 'v' in the list of vertices in 'g'
        g.vertices{end+1} = v;
        
        % Get the vertex number of the added vertex
        vertex_num = find(strcmp(v, g.vertices));
        
        % Initialize the place in 'g.edges' where the edges connecting v to other vertices will be indicated
        g.edges{vertex_num} = [ ];
    
    % Case 2: 'location' is not empty
    else
        disp(['Vertex ' v ' is already in the graph.']);
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
% Function 3 of 4: add_edge                             %
%                                                       %
%    - Purpose: To add an edge to a graph g             %
%    - Inputs                                           %
%         - g: graph structure with vertices v1 and v2  %
%         - v1: one vertex of edge to be added          %
%         - v2: another vertex of edge to be added      %
%    - Output                                           %
%         - g: structure with edge added                %
%    - Used in network_numbers (linkage classes)        %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = add_edge(g, v1, v2)
    
    % Make sure 'v1' and 'v2' are already in 'g'
    
    % Case 1: 'v1' or 'v2' is not in 'g'
    if isempty(find(strcmp(v1, g.vertices))) || isempty(find(strcmp(v2, g.vertices)))
        disp(['Make sure both ' v1 ' and ' v2 ' are in the graph.']);
        
    % Case 2: Both vertices are already in 'g'
    else
        
        % Case 2.1: The vertices are the same
        if strcmp(v1, v2) == 1
            disp(['Make sure the vertices are different.']);
            
        % Case 2.2: The vertices are different
        else
            
            % Get the index of 'v1' and 'v2' in 'g.vertices'
            v1_index = find(strcmp(v1, g.vertices));
            v2_index = find(strcmp(v2, g.vertices));
            
            % Check if an edge with 'v1' already exists
            try
                g.edges{v1_index};
            catch
                
                % If none, initialize 'g.edges' for 'v1'
                g.edges{v1_index} = cell();
            end
            
            % Check if an edge with 'v2' already exists
            try
                g.edges{v2_index};
            catch
                
                % If none, initialize 'g.edges' for 'v2'
                g.edges{v2_index} = cell();
            end
            
            % Check if an edge between 'v1' and 'v2' already exists
            for i = 1:numel(g.edges{v1_index})
                if g.edges{v1_index}(i).vertex == v2_index
%%                    disp(['Edge ' v1 '-' v2 ' already exists.']);
                    
                    % 'return' exits the function; we don't need to continue the code
                    % If we wanted to just get out of the loop, we use 'break'
                    return
                end
            end
            
            % As a control, check also if an edge between 'v2' and 'v1' already exists
            for i = 1:numel(g.edges{v2_index})
                if g.edges{v2_index}(i).vertex == v1_index
%%                    disp(['Edge ' v2 '-' v1 ' already exists.']);
                    return
                end
            end
            
            % After all the controls above have been implemented, add an edge in the 'edges' field in 'g' for both 'v1' and 'v2'
            g.edges{v1_index}(end+1) = struct('vertex', v2_index, 'label', [v1 '-' v2]);
            g.edges{v2_index}(end+1) = struct('vertex', v1_index, 'label', [v1 '-' v2]);
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                       %
% Function 4 of 4: add_path                             %
%                                                       %
%    - Purpose: To add a directed edge to a graph g     %
%    - Inputs                                           %
%         - g: graph structure with vertices v1 and v2  %
%         - v1: starting vertex of edge to be added     %
%         - v2: ending vertex of edge to be added       %
%    - Output                                           %
%         - g: structure with edge added                %
%    - Used in network_numbers (strong linkage classes) %
%                                                       %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = add_path(g, v1, v2)
    
    % Case 1: v1 or v2 is not in g
    if (isempty(find(strcmp(v1, g.vertices))) || isempty(find(strcmp(v2, g.vertices))))
        disp(['Make sure both ' v1 ' and ' v2 ' are in the graph.']);
        
    % Case 2: Both vertices are already in g
    else
        
        % Case 2.1: The vertices are the same
        if strcmp(v1, v2) == 1
            disp(['Make sure the vertices are different.']);
            
        % Case 2.2: The vertices are different
        else
            
            % Get the index of v1 and v2 in g.vertices
            v1_index = find(strcmp(v1, g.vertices));
            v2_index = find(strcmp(v2, g.vertices));
            
            % Check if an edge from v1 to v2 already exists
            for i = 1:numel(g.edges{v1_index})
                if g.edges{v1_index}(i).vertex == v2_index
%%                    disp(['Edge ' v1 '->' v2 ' already exists.']);
                    
                    % 'return' exits the function; we don't need to continue the code
                    % If we wanted to just get out of the loop, we use 'break'
                    return
                end
            end
            
            % Add a directed edge in the 'edges' field in g for v1
            g.edges{v1_index}(end+1) = struct('vertex', v2_index, 'label', [v1 '->' v2]);
        end
    end

end