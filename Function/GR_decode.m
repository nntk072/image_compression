% GR_decode.m
function X = GR_decode(p)
    % Read the BITSTREAM from the file
    fid = fopen('Output/BITSTREAM_FILE.bin', 'r');
    BITSTREAM = fread(fid, 'ubit1');
    fclose(fid);

    % Initialize the output array X
    X = zeros(size(BITSTREAM));

    i = 5; % Skip the first 4 bits which are for the parameter p
    j = 1;
    while i <= length(BITSTREAM)
        % Decode the unary code
        q = 0;
        while i <= length(BITSTREAM) && BITSTREAM(i) == 0
            q = q + 1;
            i = i + 1;
        end
        if i > length(BITSTREAM)
            break;
        end
        i = i + 1; % Skip the 1 bit

        % Decode the binary code
        r = 0;
        for k = 1:p
            if i > length(BITSTREAM)
                break;
            end
            r = r * 2 + BITSTREAM(i);
            i = i + 1;
        end

        % Compute the original value
        xi = q * 2^p + r;

        % Decode the sign bit if the original value is not zero
        if xi ~= 0 && i <= length(BITSTREAM)
            sign_bit = BITSTREAM(i);
            i = i + 1;
            if sign_bit == 1
                xi = -xi;
            end
        end

        % Store the original value in the output array X
        X(j) = xi;
        j = j + 1;
    end

    % Remove the trailing zeros in the output array X
    X = X(1:j-1);
end
