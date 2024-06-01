function CL = Segment_L(L, evec)
    n = length(evec);
    num_blocks = ceil(n / L);
    CL = 0;
    for i = 1:num_blocks
        block = evec((i-1)*L+1:min(i*L, n));
        min_CL = inf;
        for p = 0:8
            CL_block = GR_estimation(block, p);
            if CL_block < min_CL
                min_CL = CL_block;
            end
        end
        CL = CL + min_CL;
    end
end