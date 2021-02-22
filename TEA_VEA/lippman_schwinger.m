
%% eShanks given an A,b

function [dEz_for] = lippman_schwinger(Ez_for, L,U,P,Q,Mz0, deeps_str, N_order, k0, Z0, omega)
    
    
    N = size(deeps_str);
    TTopt = ones(N);
    %% ======================================================================%%
    % Ez_for if the forward simulation field

    % what's the different between Ez_proj and dEz_for
    dEz_for{1}=Ez_for; % dEz_for... the derivative of dEz...usualy computed by adjoint.



    %% it appears we are doing several forward simulations here to get the sources
    % at differenert iorders, but why?
    Mz_LS_cell{1} = Mz0;
    for iorder=2:N_order

        % even if I set TTopt to 1, the code below extracts the region where 
        % deeps_str > 0, i.e. it extracts only the region of the perturbation
        % to construct a source
        Mz = Mz0;
        % only setting source in the Topt region?
        Mz(TTopt==1)=1./(-1j*k0*Z0)*k0^2*deeps_str(TTopt==1).*dEz_for{iorder-1}(TTopt==1);
        Mz_LS_cell{iorder} = Mz;
        bprime = sparse(reshape(Mz, prod(N), 1));
        bprime = (1j*omega)*bprime;
        new_field = reshape(P.'*(U\(L\(Q.'* bprime))),N);
        dEz_for{iorder} = new_field;

    end

  

end