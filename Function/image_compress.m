% image_compress.m
function image_compress(input_file, compressed_file)
    % Read the input image
    A = imread(input_file);
    A = double(A); % Convert A from uint8 to double representation

    % Select the gray values between the rows 51 and 306, and between the columns 51 and 306
    B = A(51:306, 51:306);

    % Compute the Martucci predictor and the prediction error matrix E
    E = zeros(size(B));
    for i = 1:size(B, 1)
        for j = 1:size(B, 2)
            predictors = [];
            if i > 1
                predictors = [predictors, B(i-1, j)];
            end
            if j > 1
                predictors = [predictors, B(i, j-1)];
            end
            if i > 1 && j > 1
                predictors = [predictors, B(i-1, j-1)];
            end
            if isempty(predictors)
                y_MAP = 128;
            else
                y_MAP = median(predictors);
            end
            E(i, j) = B(i, j) - round(y_MAP);
        end
    end

    % Concatenate all the columns of E into a prediction error vector evec
    evec = E(:);

    % Find the optimal block size L
    L_values = 50:50:2000;
    CL_values = zeros(size(L_values));
    for i = 1:length(L_values)
        L = L_values(i);
        CL_values(i) = Segment_L(L, evec);
    end
    [min_CL, min_idx] = min(CL_values);
    best_L = L_values(min_idx);

    % Encode the evec using the optimal block size L
    GR_encode(evec, best_L);

    % Write the compressed image to a file
    imwrite(uint8(E + 128), compressed_file);
end

