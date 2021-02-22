
%% eShanks given an A,b

function [eShanks] = shanks_transformation(dEz_for, alpha, N_order)
    
    Shanks_cell{1} = dEz_for{1};


    %% for lippman-schwinger??
    for iorder=2:N_order

        %% initialization of the shanks transform of the field
        Shanks_cell{iorder} = Shanks_cell{iorder-1} + (alpha).^(iorder-1)*dEz_for{iorder};
    end

    eShanks{1} = Shanks_cell;
    eShanks{2} = {};
    for iorder=1:N_order-1

        eShanks{2}{iorder} = 1./(Shanks_cell{iorder+1} - Shanks_cell{iorder});

    end;

    for ip=3:N_order
        eShanks{ip} = {};
        for iorder=1:N_order-ip+1
            % this looks much like the shanks transformation
            eShanks{ip}{iorder} = eShanks{ip-2}{iorder+1} + 1./(eShanks{ip-1}{iorder+1} - eShanks{ip-1}{iorder});

        end

    end


end