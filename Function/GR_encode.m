function BITSTREAM = GR_encode(X, p)
    BITSTREAM = bitget(p, 4:-1:1); % Initialize BITSTREAM with the bits for the parameter p
    for i = 1:length(X)
        xi = X(i);
        q = floor(abs(xi) / 2^p);
        r = mod(abs(xi), 2^p);
        unary_code = [zeros(1, q), 1];
        binary_code = bitget(r, p:-1:1);
        sign_bit = xi < 0;
        if xi ~= 0
            BITSTREAM = [BITSTREAM, unary_code, binary_code, sign_bit];
        else
            BITSTREAM = [BITSTREAM, unary_code, binary_code];
        end
    end
    fid = fopen('Output/BITSTREAM_FILE.bin', 'w');
    fwrite(fid, BITSTREAM, 'ubit1');
    fclose(fid);
end
