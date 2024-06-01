clear all; addpath("Function\")
%% 1.1.
ID = 151317891;
X = mod(ID, 5);

imageFileName = sprintf('Input/Image_%d.png', X);

% Read the image using imread
A = imread(imageFileName);
% Convert A from uint8 to double representation
A = double(A);

% Plot the histogram of the image
figure(1);
imshow(A/255);
title('Original Image')

figure(2);
histogram(A(:), 'BinMethod', 'integers');
title('Histogram of the Image');
xlabel('Pixel Intensityl');
ylabel('Frequency');

saveas(gcf, 'Output/image_histogram.png')
%% 1.2.
B = A(51:306, 51:306);
B_save = B/255;
% Write the new image to a file
newImageFileName = sprintf('Output/MyImage_%d.png', X);
imwrite(B_save, newImageFileName);

figure(3)
imshow(B_save)
title('MyImage_1.png');

% Find the size of the new image file on the disk
sd = dir(newImageFileName);
fileSize = sd.bytes;
fprintf('The size of the file %s on the disk is %d bytes.\n', newImageFileName, fileSize);


figure(4);
histogram(B(:), 'BinLimits', [0, 255], 'BinWidth', 1);
title('Histogram of the New Image');
xlabel('Pixel Intensity');
ylabel('Frequency');
saveas(gcf, 'Output/new_image_histogram.png')
% Compute the empirical entropy
p = histc(B(:), 0:255) / numel(B);
p = p(p>0); % Remove zero entries
entropy = -sum(p .* log2(p));
fprintf('The empirical entropy of the image is %f.\n', entropy);

%% 1.3.
% Initialize the prediction error matrix E
E = zeros(size(B));

% Compute the Martucci predictor for each element of B
for i = 1:size(B, 1)
    for j = 1:size(B, 2)
        predictors = [];
        % Check if the northern position exists
        if i > 1
            predictors = [predictors, B(i-1, j)];
        end
        % Check if the western position exists
        if j > 1
            predictors = [predictors, B(i, j-1)];
        end
        % Check if the northwestern position exists
        if i > 1 && j > 1
            predictors = [predictors, B(i-1, j-1)];
        end

        if isempty(predictors)
            y_MAP = 128;
        else
            y_MAP = median(predictors);
        end
        
        % Define the prediction error matrix E
        E(i, j) = B(i, j) - round(y_MAP);
    end
end

% Concatenate all the columns of E into a prediction error vector evec
evec = E(:);

% Display the prediction error matrix E
figure(5)
imagesc(E);
colorbar;
title('Prediction Error Matrix E');

% Save the figure to a file
saveas(gcf, 'Output/Prediction_Error_Matrix.png');
%% 1.4.
% Initialize the minimum codelength and the optimal parameter
minCL = inf;
p_opt = 0;

% Try different p values
for p = 0:15
    % Compute the codelength for the current p
    CL = GR_estimation(evec, p);
    
    % Update the minimum codelength and the optimal parameter
    if CL < minCL
        minCL = CL;
        p_opt = p;
    end
end

% Encode X using the optimal parameter
BITSTREAM = GR_encode(evec, p_opt);

%% 1.5.
X_Test = GR_decode(3);

%% 1.6.
% Call the function Segment_L for 40 block-length values
L_values = 50:50:2000;
CL_values = zeros(size(L_values));
for i = 1:length(L_values)
    L = L_values(i);
    CL_values(i) = Segment_L(L, evec);
end

% Find the block size that produces the smallest codelength
[min_CL, min_idx] = min(CL_values);
best_L = L_values(min_idx);

% Compute bpp_CL
bpp_CL = CL_values / numel(B);

% Plot bpp_CL versus L
figure(6);
plot(L_values, bpp_CL);
title('bpp_CL versus L');
xlabel('Block Length L');
ylabel('bpp_CL');
saveas(gcf, 'Output/bpp_CL_vs_L.png');

%% 1.7
image_compress(imageFileName, 'Output/imageCompressed.png')
Brec = reconstruct_image(E);

figure(7);
imshow(Brec/255)

%% 2.1.
% Compute the histogram of the prediction error vector evec
symbols = unique(evec);
counts = histc(evec, symbols);
prob = counts / sum(counts);

% Design the Huffman dictionary
dict = huffmandict(symbols, prob);

% Encode the prediction error vector evec using the Huffman dictionary
huffcode = huffmanenco(evec, dict);

% Verify the lossless encoding using huffmandeco
evec_decoded = huffmandeco(huffcode, dict);
islossless = isequal(evec, evec_decoded);

% Compute the code length of the Huffman encoded prediction error evec
codelengths = cellfun(@length, dict);
avg_codelength = sum(prob .* codelengths);

% Compute the required code length of the Huffman dictionary
dict_codelength = sum(cellfun(@(x) length(x), dict));

% Compute the total code length
total_codelength = length(huffcode) + dict_codelength;

%% 2.2.
newImageFileNameLossless = sprintf('Output/MyImage_2.jp2');
imwrite(B_save, newImageFileNameLossless, 'Mode', 'lossless');

info_png = dir(newImageFileName);
info_lossless_png = dir(newImageFileNameLossless);
fprintf('The size of the .png file is %d bytes.\n', info_png.bytes);
fprintf('The size of the .jp2 file is %d bytes.\n', info_lossless_png.bytes);

%% 2.3.

% Define the possible block sizes
k_values = [8, 16, 32, 64, 128, 256];

% Initialize the minimum codelength and the optimal block size
min_CL = inf;
best_k = 0;

% Iterate over each possible block size
for k = k_values
    % Initialize the total codelength for this block size
    total_CL = 0;

    % Segment the prediction error matrix E into k x k blocks
    for i = 1:k:size(E, 1)
        for j = 1:k:size(E, 2)
            block = E(i:min(i+k-1, size(E, 1)), j:min(j+k-1, size(E, 2)));

            % Vectorize the block
            evec = block(:);

            % Find the optimal p for this block
            min_CL_block = inf;
            for p = 0:8
                CL_block = GR_estimation(evec, p);
                if CL_block < min_CL_block
                    min_CL_block = CL_block;
                end
            end

            % Add the codelength of this block to the total codelength
            total_CL = total_CL + min_CL_block;
        end
    end

    % If the total codelength for this block size is smaller than the current minimum,
    % update the minimum codelength and the optimal block size
    if total_CL < min_CL
        min_CL = total_CL;
        best_k = k;
    end
end

fprintf('The optimal block size is %d x %d.\n', best_k, best_k);

%% 2.4.

% Define the possible block sizes
k_values = [256, 128, 64, 32, 16, 8];

% Initialize the quadtree
quadtree = cell(size(E));

% Iterate over each possible block size
for k = k_values
    % Segment the prediction error matrix E into k x k blocks
    for i = 1:k:size(E, 1)
        for j = 1:k:size(E, 2)
            block = E(i:min(i+k-1, size(E, 1)), j:min(j+k-1, size(E, 2)));

            % Vectorize the block
            evec = block(:);

            % Find the optimal p for this block
            min_CL_block = inf;
            for p = 0:8
                CL_block = GR_estimation(evec, p);
                if CL_block < min_CL_block
                    min_CL_block = CL_block;
                end
            end

            % If the block size is not the smallest possible size
            if k > 8
                % Segment the block into 4 smaller blocks
                smaller_blocks = mat2cell(block, [k/2, k/2], [k/2, k/2]);

                % Compute the total codelength of the 4 smaller blocks
                total_CL_smaller_blocks = 0;
                for m = 1:2
                    for n = 1:2
                        smaller_block = smaller_blocks{m, n};
                        evec_smaller_block = smaller_block(:);
                        min_CL_smaller_block = inf;
                        for p = 0:8
                            CL_smaller_block = GR_estimation(evec_smaller_block, p);
                            if CL_smaller_block < min_CL_smaller_block
                                min_CL_smaller_block = CL_smaller_block;
                            end
                        end
                        total_CL_smaller_blocks = total_CL_smaller_blocks + min_CL_smaller_block;
                    end
                end

                % If the total codelength of the 4 smaller blocks is smaller than the codelength of the block
                if total_CL_smaller_blocks < min_CL_block
                    % Keep the split into smaller blocks
                    quadtree{i, j} = smaller_blocks;
                else
                    % Keep the block
                    quadtree{i, j} = block;
                end
            else
                % Keep the block
                quadtree{i, j} = block;
            end
        end
    end
end

% Initialize the output image
Brec = zeros(size(E));

% Iterate over each cell in the quadtree
for i = 1:size(quadtree, 1)
    for j = 1:size(quadtree, 2)
        % If the cell contains a cell array (i.e., it was split into smaller blocks)
        if iscell(quadtree{i, j})
            smaller_blocks = quadtree{i, j};
            for m = 1:2
                for n = 1:2
                    smaller_block = smaller_blocks{m, n};
                    avg_value = mean(smaller_block(:));
                    Brec((i-1)*k/2+m, (j-1)*k/2+n) = avg_value;
                end
            end
        else
            % If the cell contains a matrix (i.e., it was not split)
            block = quadtree{i, j};
            avg_value = mean(block(:));
            Brec(i:i+k-1, j:j+k-1) = avg_value;
        end
    end
end

% Display the output image
figure(8)
imshow(Brec, []);
