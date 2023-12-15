function virusGame = GameOfLife(matrix_x,matrix_y, start_row, start_col)
   % inputs: 
        % matrix_x: size of the matrix in x dir (int) 
        % matrix_y: size of the matrix in y dir (int)
        % start_row: where to start the infection in y, optional (int)
        % start_col: where to start the infection in x, optional (int) 
   % see if you can add booleans or flags for different rules 

   matrix = zeros(matrix_x, matrix_y); % create matrix based on input
   %matrix_len = max(matrix_x, matrix_y);  % get the length of the matrix 

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
       disp(matrix) % different way to visualize, seeing actual matrix after each iteration
       spread_history{iteration} = matrix; % store the current state for visualization

       new_infected = zeros(size(matrix)); % matrix to track newly infected cells

       [row, col] = find(matrix == 1);

       % check and update neighboring cells
       neighborhood = [-1, 0, 1]; % allowed within a 1 cell limit 
       for i = 1:length(row) % loop through the row
           for dr = neighborhood  % examine neighboring cells 
               for dc = neighborhood
                   if dr == 0 && dc == 0 % checking if it's the same cell
                       continue; % skip the current cell
                   end

                   new_row = row(i) + dr;     % update new row index
                   new_col = col(i) + dc;     % update new column index 

                   % check if the neighboring cell is within bounds
                   if new_row >= 1 && new_row <= matrix_x && new_col >= 1 && new_col <= matrix_y
                       % Update the neighboring cell if it's '0'
                       if matrix(new_row, new_col) == 0 % if it is 0, update
                           matrix(new_row, new_col) = 1;
                           new_infected(new_row, new_col) = 1; % mark as newly infected
                       end
                   end
               end
           end
       end
       % mark previously infected cells as '2' for recovered 
       matrix(matrix == 1) = 2;

       % check if there are no new infections in this iteration
       if nnz(new_infected) == 0
           break; % Break the loop if no new infections occurred
       end

       % update the newly infected cells for the next iteration
       matrix(new_infected == 1) = 1;

       iteration = iteration + 1;
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