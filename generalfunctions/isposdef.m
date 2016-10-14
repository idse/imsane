function s = isposdef(M)
    % Is matrix positive definite?
    %
    % s = isposdef(M)
    %
    % s:    boolean
    % M:    matrix

    [~,s] = chol(M);
    s = ~s;
end