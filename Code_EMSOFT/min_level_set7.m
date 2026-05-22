function minsize = min_level_set7(P,M,n)
    sdpvar x1 x2 x3 x4 x5 x6 x7;
    sdpvar alph;
    sdpindets = [x1 x2 x3 x4 x5 x6 x7];
    sdpvector = monolist([x1;x2;x3;x4;x5;x6;x7], [2]);
    % We want to leave off the constant term
    sdpvector = sdpvector(2:8);
    search_radius = n;

    V = sdpvector'*P*sdpvector;
    W = sdpvector'*M*sdpvector;
    sdpvar gam lambda;
    %F = sos(V-lambda*(W-search_radius)-gam);
    F = sos(-V+lambda*(W-search_radius)+gam);
    F = [F lambda>=0];
    [sol,v,Q] = solvesos(F,gam,sdpsettings('verbose',0));
    minsize = double(gam);