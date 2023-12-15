function onlySomeSpread = GameOfLifeSome(matrix_x, matrix_y, infected_num)
    % purpose: simulate disease spread, but limit the number able to be
    % infected
    % inputs 
        % matrix_x: size of the matrix in x dir (int) 
        % matrix_y: size of the matrix in y dir (int)
        % infected_num: number of individuals the infected can spread to
        % (int) 

    matrix = zeros(matrix_x, matrix_y); % define size
    matrix_len = height(matrix);  % get the length of the matrix 
    
    % find a starting spot on the matrix
    start_row = randi(matrix_len);    
    start_col = randi(matrix_len);
    
    % mark the infected individual with a 1 
    matrix(start_col, start_row) = 1;
    
    % store each iteration's matrix to visualize the spread
    spread_history = cell(1, matrix_x*matrix_y); % cell array based on input
    
    % continue spreading until no new cells are infected in an iteration
    iteration = 1;
    while true
        disp(['Iteration: ' num2str(iteration)])
        disp(matrix) % different way to visualize, seeing actual matrix after each iteration
        spread_history{iteration} = matrix; % store the current state for visualization
        new_infected = zeros(size(matrix)); % matrix to track newly infected cells
        
        [row, col] = find(matrix == 1); % starting location
        
        % check and update neighboring cells
        neighborhood = [-1, 0, 1]; % allowed within a 1 cell limit 
        for i = 1:length(row)   % loop through the row
            infected_count = 0; % track the number of infections caused by each infected cell
            for dr = neighborhood
                for dc = neighborhood
                    if dr == 0 && dc == 0 % check if it's the same cell
                        continue; % skip the current cell
                    end
                    % update positions
                    new_row = row(i) + dr;
                    new_col = col(i) + dc;
                    
                    % check if the neighboring cell is within bounds
                    if new_row >= 1 && new_row <= matrix_len && new_col >= 1 && new_col <= matrix_len
                        % update the neighboring cell if it's 0
                        if matrix(new_row, new_col) == 0 && infected_count < infected_num % make sure 
                            matrix(new_row, new_col) = 1; 
                            new_infected(new_row, new_col) = 1; % mark as newly infected
                            infected_count = infected_count + 1; % increment the infection count
                        end
                    end
                    
                    if infected_count >= infected_num
                        break; % stop infecting more cells if the limit is reached
                    end
                end
                
                if infected_count >= infected_num
                    break; % stop infecting
                end
            end
        end
        
        % mark previously infected cells as '2' for recovered
        matrix(matrix == 1) = 2;
        
        % xheck if there are no new infections in this iteration
        if nnz(new_infected) == 0
            break; % break the loop 
        end
        
        % update only the newly infected cells for the next iteration
        matrix(new_infected == 1) = 1;   
        iteration = iteration + 1; % update iteration count
    end
   % add visualization 
   figure;
   colormap([1 1 1; 1 0 0; 0.5 0.5 0.5]); % define colormap: red, white, grey

   for iter = 1:iteration  
       matrix = spread_history{iter};   % loop through spread_history var  
       subplot(1, iteration, iter); % create subplots for each iteration
       imagesc(matrix); % display the matrix as an image
       title(['Iteration: ' num2str(iter)]); % label
       axis square; % manipulate proportions
   end
end 