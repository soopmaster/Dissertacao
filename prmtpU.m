function dalph = prmtpU(alph,Dfxo)
    % Variables for ODE: x = alph(:,k) (vectorised)
    %   U(xo,alph,y) = -alph*Df(xo)

    % Input variables:
    %   Dfxo_vec: vectorised jacobian in point xo

    %alph = reshape(x,dim,[]);
    %Dfxo = reshape(Dfxo_vec,dim,[]);

    dalph = alph*Dfxo;

    %dxdt = reshape(dalph,[],1);

end