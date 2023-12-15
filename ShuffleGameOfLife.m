function randomVirusGame = ShuffleGameOfLife(matrix_x, matrix_y, start_row, start_col)
   % inputs: 
        % matrix_x: size of the matrix in x dir (int) 
        % matrix_y: size of the matrix in y dir (int)
        % start_row: where to start the infection in y, optional (int)
        % start_col: where to start the infection in x, optional (int) 
   matrix = zeros(matrix_x, matrix_y); % create matrix based on input
   matrix_len = max(matrix_x, matrix_y);  % get the length of the matrix 

   % find a starting spot on the matrix
   if nargin < 3 || isempty(start_row) || isempty(start_col)
       start_row = randi(matrix_x);   
       start_col = randi(matrix_y);
   end 
   % denote the infected individual with a 1 
   matrix(start_row, start_col) = 1;

   % store each iteration's matrix to visualize the spread
   spread_history = cell(1, matrix_x*matrix_y); % cell array based on size of input matrix

   % continue spreading until no new cells are infected in an iteration
   iteration = 1;
    while true
        disp(['Iteration: ' num2str(iteration)])
        disp(matrix) % alternate displau
        spread_history{iteration} = matrix; % store the current state for visualization
        
        new_infected = zeros(size(matrix)); % matrix to track newly infected cells
        % find current position of infected
        [row, col] = find(matrix == 1);
        
        % Check and update neighboring cells
        neighborhood = [-1, 0, 1]; % define neighbors, able to look/infect 1 cell away
        for i = 1:length(row)   % iterate through
            for dr = neighborhood
                for dc = neighborhood
                    if dr == 0 && dc == 0   % check if in current position
                        continue; % skip the current cell
                    end
                    % update position
                    new_row = row(i) + dr;
                    new_col = col(i) + dc;
                    
                    % check if the neighboring cell is within bounds
                    if new_row >= 1 && new_row <= matrix_len && new_col >= 1 && new_col <= matrix_len
                        % update the neighboring cell if it's 0
                        if matrix(new_row, new_col) == 0
                            matrix(new_row, new_col) = 1;
                            new_infected(new_row, new_col) = 1; % mark as newly infected
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
        % apply the shuffle to the matrix 
        [m, n] = size(matrix);  % return the x*y of matrix
        [~, idx] = sort(rand(m,n));
        matrix= matrix(sub2ind([m, n], idx, ones(m,1)*(:n)));  
    
        iteration = iteration + 1; % add to iteration count 
    end
    figure;
    colormap([1 1 1; 1 0 0; 0.5 0.5 0.5]); % define colormap: white, red, grey
    
    % static 
    for iter = 1:iteration  % 
        matrix = spread_history{iter};
        subplot(1, iteration, iter); % Create subplots for each iteration
        imagesc(matrix); % Display the matrix as an image
        title(['Iteration: ' num2str(iter)]);
        axis square;
        pause(0.1); % Pause to see each iteration, adjust as needed
    end
end 