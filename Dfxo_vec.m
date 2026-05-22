function Dfxo=Dfxo_vec(Df,vars,xo)
    Dfxo = double(subs(Df,vars,xo));
    Dfxo = reshape(Dfxo,[],1);
end