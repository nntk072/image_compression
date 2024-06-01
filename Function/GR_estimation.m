function CL = GR_estimation(X, p)
    % Compute the codelength
    CL = sum(p + floor(abs(X) / 2^p) + 1 + (X ~= 0));
    
    % Add 4 to account for the need to transmit the value of p
    CL = CL + 4;
end
