function vaccineGameOfLife = VaccinatedGame(matrix_x, matrix_y, vaccinated_num, start_col, start_row)
    % inputs: 
        % matrix_x: size of the matrix in x dir (int) 
        % matrix_y: size of the matrix in y dir (int)
        % vaccinated_num: number of individuals vaccinated (int)
        % start_row: where to start infection in y (int) 
        % start_col: where to start infection in x (int)

    matrix = zeros(matrix_x, matrix_y); % initial size of matrix based on user input 
    matrix_len = height(matrix);  % get the length of the matrix 
    
    % find a starting spot on the matrix, random spot if not defined
    if nargin < 4 || isempty(start_col) || isempty(start_row)
        start_row = randi(matrix_len);    
        start_col = randi(matrix_len);
    end 
    
    %  mark the infected individual with a 1 
    matrix(start_col, start_row) = 1;
    
    % randomly vaccinate individuals
    vaccinated = 0;
    while vaccinated < vaccinated_num
        vac_row = randi(matrix_len);
        vac_col = randi(matrix_len);
        if matrix(vac_row, vac_col) == 0
            matrix(vac_row, vac_col) = -1; % denote vaccinated with -1
            vaccinated = vaccinated + 1;
        end
    end
    
    % store each iteration's matrix to visualize the spread
    spread_history = cell(1, matrix_x*matrix_y); % store matrices in a cell array
    
    % continue spreading until no new cells are infected in an iteration
    iteration = 1;
    while true
        disp(['Iteration: ' num2str(iteration)])
        disp(matrix)
        spread_history{iteration} = matrix; % store the current state for visualization
        
        new_infected = zeros(size(matrix)); % matrix to track newly infected cells
        
        [row, col] = find(matrix == 1);
        
        % check and update neighboring cells
        neighborhood = [-1, 0, 1];
        for i = 1:length(row) % loop through row 
            infected_count = 0; % track the number of infections caused by each infected cell
            for dr = neighborhood  
                for dc = neighborhood
                    if dr == 0 && dc == 0 % check if it's the current spot
                        continue; % skip the current cell
                    end
                    % update positions 
                    new_row = row(i) + dr;
                    new_col = col(i) + dc;
                    
                    % check if the neighboring cell is within bounds
                    if new_row >= 1 && new_row <= matrix_len && new_col >= 1 && new_col <= matrix_len
                        % update the neighboring cell if it's '0' and not vaccinated
                        if matrix(new_row, new_col) == 0 && infected_count < 2 && matrix(new_row, new_col) ~= -1
                            matrix(new_row, new_col) = 1;
                            new_infected(new_row, new_col) = 1; % Mark as newly infected
                            infected_count = infected_count + 1; % Increment the infection count
                        end
                    end
                end
            end
        end
        
        % mark previously infected cells as '2'
        matrix(matrix == 1) = 2;
        
        % check if there are no new infections in this iteration
        if nnz(new_infected) == 0
            break; % break the loop if no new infections occurred
        end
        
        % update only the newly infected cells for the next iteration
        matrix(new_infected == 1) = 1;
        
        iteration = iteration + 1;
        %disp(matrix)
    end
    % make visualization
    figure;
    colormap([0 1 0; 1 1 1; 1 0 0; 0.5 0.5 0.5]); % define colormap: green, white, red, grey

    for iter = 1:iteration  % 
        matrix = spread_history{iter};
        subplot(1, iteration, iter); % Create subplots for each iteration
        imagesc(matrix); % Display the matrix as an image
        title(['Iteration: ' num2str(iter)]);
        axis square;
        pause(0.1); % Pause to see each iteration, adjust as needed
    end
end 