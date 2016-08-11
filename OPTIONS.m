global optionsZ


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZIGZAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Measurement locations
l1 = 'location';    v1 = 1;
                    % 1 calculate location based on vertex utm
                    % 2 calculate location based on last probe point
l2 = 'z';           v2 = 2;
                    % 1 depth value (raw)
                    % 2 SWE value using mean SWE tube value at zigzag

optionsZ = struct(l1,v1,l2,v2);
    clear l* v*