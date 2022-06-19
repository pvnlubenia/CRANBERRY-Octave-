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



function [model] = network_numbers(model)
    
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
    reactant_complexes = [ ];
    
    % Initialize the matrix of product complexes
    product_complexes = [ ];
    
    % Initialize the stoichiometric matrix
    N = [ ];
    
    % For each reaction in the model
    for i = 1:numel(model.reaction)
        
        % Initialize the vector for the reaction's reactant complex
        reactant_complexes(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the reactant complex
        for j = 1:numel(model.reaction(i).reactant)
            reactant_complexes(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
        end
        
        % Initialize the vector for the reaction's product complex
        product_complexes(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the product complex
        for j = 1:numel(model.reaction(i).product)
            product_complexes(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
        end
        
        % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
        N(:, end+1) = product_complexes(:, end) - reactant_complexes(:, end); % N %
        
        % If the reaction is reversible
        if model.reaction(i).reversible
            
            % Insert a new vector for the reactant complex: make it same as the product complex
            reactant_complexes(:, end+1) = product_complexes(:, end);
            
            % Insert a new vector for the product complex: make it the same as the reactant complex
            product_complexes(:, end+1) = reactant_complexes(:, end-1);
            
            % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
            N(:, end+1) = -N(:, end);
        end
    end
    
    % Get just the unique complexes
    % ind2(i) is the index in Y of the reactant complex in reaction i
    [Y, ~, ind2] = unique([reactant_complexes product_complexes]', 'rows');
    
    % Construct the matrix of complexes
    Y = Y';
    
    % Count the number of complexes
    n = size(Y, 2);
    
    
    
    %
    % Reactant complexes
    %
    
    % Get just the unique reactant complexes
    reactant_complexes_unique = unique([reactant_complexes]', 'rows')';
    
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
    
    % Initialize a matrix (complexes x complexes) for the reacts_to relation
    % This is for testing reversibility of the network
    reacts_to = false(n, n);
    
    % Initialize matrix (complexes x total reactions) for the reacts_in relation
    % This is the incidence matrix I_a
    reacts_in = zeros(n, r);
    
    % Fill out the entries of the matrices
    for i = 1:r
        
        % reacts_to(i, j) = true iff there is a reaction r: y_i -> y_j
        reacts_to(ind2(i), ind2(i + r)) = true;
        
        % reacts_in(i, r) = -1 and reacts_in(j, r) = 1) iff there is a reaction r: y_i -> y_j
        reacts_in(ind2(i), i) = -1;
        reacts_in(ind2(i+r), i) = 1;
    end
    
    % Linkage classes
    linkage_class = connected_components(umultigraph(reacts_to | reacts_to'));
    
    % Count the number of linkage classes
    l = max(linkage_class);
    
    
    
    %
    % Strong linkage classes
    %
    
    % Check if the network is reversibile
    is_reversible = isequal(reacts_to, reacts_to');
    
    % Strong linkage classes
    if is_reversible
        strong_linkage_class = linkage_class;
    else
        strong_linkage_class = strongly_connected_components(multigraph(reacts_to));
    end
    
    % Count the number of strong linkage classes
    sl = max(strong_linkage_class);
    
    
    
    %
    % Terminal linkage classes
    %
    
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



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 1 of 5: umultigraph (taken from [2])                       %
%                                                                     %
%    - Purpose: To fix duplicate edges in an undirected graph         %
%    - Input                                                          %
%         - g: undirected graph                                       %
%    - Output                                                         %
%         - g: undirected graph where edges with the same endpoints   %
%              and the same label, is substituted with an edge with   %
%              weight equal to the sum of their weights               %
%    - Used in umultigraph                                            %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = fix_dup_edges(g)
    
    fixed_edges = 0;
    for i = 1:numel(g.vertices)
        j = 1;
        while j <= numel(g.edges{i})-1
            v2 = g.edges{i}(j).vertex;
            label = g.edges{i}(j).label;
            k = j + 1;
            while k <= numel(g.edges{i})
                if v2 == g.edges{i}(k).vertex && strcmp(label, g.edges{i}(k).label)
                    g.edges{i}(j).weight = g.edges{i}(j).weight + g.edges{i}(k).weight;
                    g.edges{i}(k) = [];
                    if i < v2
                        disp(['Duplicated edge: ' g.vertices{i} ' -(' label ')-> ' g.vertices{v2}])
                    end
                    fixed_edges = fixed_edges + 1;
                else
                    k = k + 1;
                end
            end
            j = j + 1;
        end
    end
    if fixed_edges > 0
        disp(['Summary: fixed ' int2str(fixed_edges / 2) ' duplicated edges'])
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 2 of 5: umultigraph (taken from [2])                       %
%                                                                     %
%    - Purpose: To initialize an undirected graph                     %
%    - Input: none                                                    %
%    - Output                                                         %
%         - g: empty undirected graph                                 %
%    - Used in network_numbers (linkage class computation)            %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function ug = umultigraph(varargin)
    
    if nargin == 0
        ug.vertices = {};
        ug.edges = {};
    elseif nargin == 1 && isa(varargin{1}, 'umultigraph')
        ug = varargin{1};
    elseif nargin == 1 && isreal(varargin{1}) && ~ischar(varargin{1}) && isequal(varargin{1}, varargin{1}')
        ug.vertices = cell(size(varargin{1}, 1), 1);
        for i = 1:size(varargin{1}, 1)
            ug.vertices{i} = num2str(i);
        end
        ug.edges = cell(numel(ug.vertices), 1);
        for i = 1:size(varargin{1}, 1)
            for j = i:size(varargin{1}, 2)
                if varargin{1}(i, j)
                    ug.edges{i}(end+1) = struct('vertex', j, 'label', '', 'weight', varargin{1}(i, j));
                    if i ~= j
                        ug.edges{j}(end+1) = struct('vertex', i, 'label', '', 'weight', varargin{1}(i, j));
                    end
                end
            end
        end
    elseif nargin == 2 && isvector(varargin{1}) && iscellstr(varargin{1}) && ndims(varargin{2}) == 2 && iscell(varargin{2}) && size(varargin{2}, 2) >= 3 && size(varargin{2}, 2) <= 4 && iscellstr(varargin{2}(:, 1:3))
        if size(varargin{2}, 2) == 4 && (~isreal(cell2mat(varargin{2}(:, 4))) || ischar(cell2mat(varargin{2}(:, 4))))
            error('Wrong arguments: the edge weights should be real numbers.')
        end
        ug.vertices = varargin{1};
        ug.edges = cell(numel(ug.vertices), 1);
        for i = 1:size(varargin{2}, 1)
            v1 = find(strcmp(varargin{2}{i, 1}, ug.vertices), 1);
            v2 = find(strcmp(varargin{2}{i, 2}, ug.vertices), 1);
            if size(varargin{2}, 2) == 3
                ug.edges{v1}(end+1) = struct('vertex', v2, 'label', varargin{2}{i, 3}, 'weight', 1);
                if v1 ~= v2
                    ug.edges{v2}(end+1) = struct('vertex', v1, 'label', varargin{2}{i, 3}, 'weight', 1);
                end
            else
                ug.edges{v1}(end+1) = struct('vertex', v2, 'label', varargin{2}{i, 3}, 'weight', varargin{2}{i, 4});
                if v1 ~= v2
                    ug.edges{v2}(end+1) = struct('vertex', v1, 'label', varargin{2}{i, 3}, 'weight', varargin{2}{i, 4});
                end
            end
        end
        ug = fix_dup_edges(ug);
    else
        error('Wrong arguments.');
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 3 of 5: connected_components (taken from [2])              %
%                                                                     %
%    - Purpose: To determine to which component each vertex belongs   %
%         to                                                          %
%    - Input                                                          %
%         - g: graph with edges and vertices                          %
%    - Output                                                         %
%         - cc: list of component number to which each vertex belongs %
%         to                                                          %
%    - Used in network_numbers (linkage class computation)            %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function cc = connected_components(g)

    function visit(v)
        cc(v) = n_cc;
        to_finish = [v];
        while ~isempty(to_finish)
            v1 = to_finish(end);
            to_finish(end) = [ ];
            for j = 1:numel(g.edges{v1})
                v2 = g.edges{v1}(j).vertex;
                if cc(v2) == 0 % not yet visited
                    cc(v2) = n_cc;
                    to_finish(end+1) = v2;
                end
            end
        end
    end
    
    cc = zeros(numel(g.vertices), 1);
    n_cc = 0;
    for i = 1:numel(g.vertices)
        if cc(i) == 0 % not yet visited
            n_cc = n_cc + 1;
            visit(i);
        end
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 4 of 5: multigraph (taken from [2])                        %
%                                                                     %
%    - Purpose: To initialize a directed graph                        %
%    - Input: none                                                    %
%    - Output                                                         %
%         - g: empty directed graph                                   %
%    - Used in network_numbers (strong linkage class computation)     %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function g = multigraph(varargin)
    
    if nargin == 0
        g.vertices = {};
        g.edges = {};
    elseif nargin == 1 && isa(varargin{1}, 'multigraph')
        g = varargin{1};
    elseif nargin == 1 && isreal(varargin{1}) && ~ischar(varargin{1}) && ndims(varargin{1}) == 2 && size(varargin{1}, 1) == size(varargin{1}, 2)
        g.vertices = cell(size(varargin{1}, 1), 1);
        for i = 1:size(varargin{1}, 1)
            g.vertices{i} = num2str(i);
        end
        g.edges = cell(numel(g.vertices), 1);
        for i = 1:size(varargin{1}, 1)
            for j = 1:size(varargin{1}, 2)
                if varargin{1}(i, j)
                    g.edges{i}(end+1) = struct('vertex', j, 'label', '', 'weight', varargin{1}(i, j));
                end
            end
        end
    elseif nargin == 2 && isvector(varargin{1}) && iscellstr(varargin{1}) && ndims(varargin{2}) == 2 && iscell(varargin{2}) && size(varargin{2}, 2) >= 3 && size(varargin{2}, 2) <= 4 && iscellstr(varargin{2}(:, 1:3))
        if size(varargin{2}, 2) == 4 && (~isreal(cell2mat(varargin{2}(:, 4))) || ischar(cell2mat(varargin{2}(:, 4))))
            error('Wrong arguments: the edge weights should be real numbers.')
        end
        g.vertices = varargin{1};
        g.edges = cell(numel(g.vertices), 1);
        for i = 1:size(varargin{2}, 1)
            v1 = find(strcmp(varargin{2}{i, 1}, g.vertices), 1);
            v2 = find(strcmp(varargin{2}{i, 2}, g.vertices), 1);
            if size(varargin{2}, 2) == 3
                g.edges{v1}{end+1} = struct('vertex', v2, 'label', varargin{2}{i, 3}, 'weight', 1);
            else
                g.edges{v1}{end+1} = struct('vertex', v2, 'label', varargin{2}{i, 3}, 'weight', varargin{2}{i, 4});
            end
        end
    else
        error('Wrong arguments.');
    end

end



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                     %
% Function 5 of 5: strongly_connected_components (taken from [2])     %
%                                                                     %
%    - Purpose: To determine to which strongly connected component    %
%         each vertex belongs to                                      %
%    - Input                                                          %
%         - g: graph with edges and vertices                          %
%    - Output                                                         %
%         - scc: list of strongly connected component number to which %
%         each vertex belongs to                                      %
%    - Used in network_numbers (strong linkage class computation)     %
%                                                                     %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function scc = strongly_connected_components(g)
    
    function visit(v)
        dtime(v) = time;
        scc_dtime(v) = time;
        time = time + 1;
        stack(end+1) = v;
        on_stack(v) = true;
        for j = 1:numel(g.edges{v})
            v2 = g.edges{v}(j).vertex;
            if dtime(v2) == 0
                visit(v2);
                scc_dtime(v) = min(scc_dtime(v), scc_dtime(v2)); % scc_dtime(v2) < scc_dtime(v) iff a vertex in the stack before v is reachable from v2 (and so they are all in the same scc)
            elseif on_stack(v2) % v2 was visited before v
                scc_dtime(v) = min(scc_dtime(v), scc_dtime(v2)); % so v and v2 are in the same component, and they must have the same scc_dtime
            end
        end
        if scc_dtime(v) == dtime(v) % v is the first visited node of its scc, all the other vertices of the scc follow it on the stack
            while true
                v2 = stack(end);
                stack(end) = [];
                on_stack(v2) = false;
                scc(v2) = n_cc;
                if v2 == v
                    break
                end
            end
            n_cc = n_cc + 1;
        end
    end
    
    dtime = zeros(numel(g.vertices), 1); % vertex discovery times by visit()
    scc_dtime = zeros(numel(g.vertices), 1); % discovery time of (the first visited node of) the vertex scc
    time = 1;
    stack = [];
    on_stack = false(numel(g.vertices), 1);
    scc = zeros(numel(g.vertices), 1);
    n_cc = 1;
    for i = 1:numel(g.vertices)
        if scc(i) == 0
            visit(i);
        end
    end
 
end