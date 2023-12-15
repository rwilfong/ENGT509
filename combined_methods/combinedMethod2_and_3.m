function combinedMethod2_and_3(matrix_x, matrix_y, infected_num, vaccinated_num, start_col, start_row)
    matrix = zeros(matrix_x, matrix_y); % initialize matrix
    matrix_len = size(matrix, 1); % length of the matrix

    % starting spot for infection
    if nargin < 6 || isempty(start_col) || isempty(start_row)
        start_row = randi(matrix_len);
        start_col = randi(matrix_len);
    end

    % set infected individual
    matrix(start_col, start_row) = 1;

    % randomly vaccinate individuals
    vaccinated = 0;
    while vaccinated < vaccinated_num
        vac_row = randi(matrix_len);
        vac_col = randi(matrix_len);
        if matrix(vac_row, vac_col) == 0
            matrix(vac_row, vac_col) = -1; % Denote vaccinated with -1
            vaccinated = vaccinated + 1;
        end
    end

    % store iteration's matrix for visualization
    spread_history = cell(1, matrix_x * matrix_y);

    % continue spreading until no new cells are infected in an iteration
    iteration = 1;
    while true
        disp(['Iteration: ' num2str(iteration)])
        disp(matrix)
        spread_history{iteration} = matrix;

        new_infected = zeros(size(matrix)); % Track newly infected cells

        [row, col] = find(matrix == 1);

        % check and update neighboring cells
        neighborhood = [-1, 0, 1];
        for i = 1:length(row)
            infected_count = 0; % track infections caused by each infected cell
            for dr = neighborhood
                for dc = neighborhood
                    if dr == 0 && dc == 0
                        continue; % skip the current cell
                    end
                    % update positions
                    new_row = row(i) + dr;
                    new_col = col(i) + dc;

                    % check if neighboring cell is within bounds
                    if new_row >= 1 && new_row <= matrix_len && new_col >= 1 && new_col <= matrix_len
                        % update the neighboring cell based on conditions
                        if matrix(new_row, new_col) == 0 && infected_count < infected_num && matrix(new_row, new_col) ~= -1
                            matrix(new_row, new_col) = 1;
                            new_infected(new_row, new_col) = 1; % mark as newly infected
                            infected_count = infected_count + 1; % increment the infection count
                        end
                    end
                end
            end
        end

        % mark previously infected cells as '2'
        matrix(matrix == 1) = 2;

        % check for new infections in this iteration
        if nnz(new_infected) == 0
            break; % break the loop if no new infections occurred
        end

        % update only the newly infected cells for the next iteration
        matrix(new_infected == 1) = 1;

        iteration = iteration + 1;
    end

    % make visualization
    figure;
    colormap([0 1 0; 1 1 1; 1 0 0; 0.5 0.5 0.5]); % Define colormap: green, white, red, grey

    for iter = 1:iteration
        matrix = spread_history{iter};
        subplot(1, iteration, iter); % create subplots for each iteration
        imagesc(matrix); % display the matrix as an image
        title(['Iteration: ' num2str(iter)]);
        axis square;
    end
end
