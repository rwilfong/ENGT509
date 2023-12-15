function combined_method5_method3(matrix_x, matrix_y, transmission_rate, vaccinated_num, start_col, start_row)
    % inputs 
    % matrix_x: size of the matrix in x dir (int)
    % matrix_y: size of the matrix in y dir (int)
    % transmission_rate: probability of transmission from I to S (float, 0-1)
    % vaccinated_num: number of individuals vaccinated (int)
    % start_col: where to start the infection in x, optional (int)
    % start_row: where to start the infection in y, optional (int)

    matrix = zeros(matrix_x, matrix_y); % create matrix based on input
    matrix_len = max(matrix_x, matrix_y);  % get the length of the matrix 
    
    % define starting position for infected
    if nargin < 6 || isempty(start_row) || isempty(start_col)
        start_row = randi(matrix_x);   
        start_col = randi(matrix_y);
    end 
    
    % mark the infected individual with a 1 
    matrix(start_row, start_col) = 1;

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
    spread_history = cell(1, matrix_x*matrix_y); % cell array based on size of input matrix

    iteration = 1;
    while true
        disp(['Iteration: ' num2str(iteration)])
        disp(matrix) % different way to visualize, seeing actual matrix after each iteration
        spread_history{iteration} = matrix; % store the current state for visualization

        new_infected = zeros(size(matrix)); % matrix to track newly infected cells

        [row, col] = find(matrix == 1);

        % check and update neighboring cells
        neighborhood = [-1, 0, 1]; % allowed within a 1 cell limit 
        for i = 1:length(row)
            for dr = neighborhood   % examine neighboring cells
                for dc = neighborhood
                    if dr == 0 && dc == 0   % check if it's the current cell
                        continue; % skip the current cell
                    end
                    % update position
                    new_row = row(i) + dr;
                    new_col = col(i) + dc;

                    % check if the neighboring cell is within bounds
                    if new_row >= 1 && new_row <= matrix_len && new_col >= 1 && new_col <= matrix_len
                        % update the neighboring cell probabilistically
                        if matrix(new_row, new_col) == 0 && rand <= transmission_rate && matrix(new_row, new_col) ~= -1
                            matrix(new_row, new_col) = 1;
                            new_infected(new_row, new_col) = 1; % Mark as newly infected
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
    end

    % create visualization 
    figure;
    colormap([0 1 0; 1 1 1; 1 0 0; 0.5 0.5 0.5]); % define colormap: green, white, red, grey

    for iter = 1:iteration  
        matrix = spread_history{iter};
        subplot(1, iteration, iter); % create subplots for each iteration
        imagesc(matrix); % display the matrix as an image
        title(['Iteration: ' num2str(iter)]);
        axis square;
    end
end
