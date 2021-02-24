
%% eShanks given an A,b

function [eShanks] = vector_shanks_transformation(dEz_for, alpha, N_order)
    

    Shanks_cell{1} = dEz_for{1};


    %% for lippman-schwinger??
    for iorder=2:N_order


        %% initialization of the shanks transform of the field
        Shanks_cell{iorder} = Shanks_cell{iorder-1} + (alpha).^(iorder-1)*dEz_for{iorder};
    end
    y  = abs(10*randn(size(dEz_for{1})));

    eShanks{1} = Shanks_cell;
    eShanks{2} = {};
    for iorder=1:N_order-1

        delta1n = eShanks{1}{iorder+1}-eShanks{1}{iorder};
        denominator = conj(delta1n).*y;
        teaInverse = y./denominator;    
        eShanks{2}{iorder} = teaInverse;
    end

    for ip=3:N_order
        eShanks{ip} = {};
        for iorder=1:N_order-ip+1
            % this looks much like the shanks transformation
            %% for tea, we have to seperate out even and odd in ip or iorder
            if(mod(ip-1,2) == 0)%2k+2
                delta1n = eShanks{ip-1}{iorder+1}-eShanks{ip-1}{iorder};
                delta0n = eShanks{ip-2}{iorder+1}-eShanks{ip-2}{iorder};
                bilinear = conj(delta1n).*delta0n; %THIS ORDER IS IMPORTANT
                teaInverse = delta0n./bilinear;
            else %2k+1, odd cases
                delta1n = eShanks{ip-1}{iorder+1}-eShanks{ip-1}{iorder};
                denominator = conj(delta1n).*y;            
                teaInverse = y./denominator;
            end
            eShanks{ip}{iorder} = eShanks{ip-2}{iorder+1} + teaInverse;
        end

    end


end