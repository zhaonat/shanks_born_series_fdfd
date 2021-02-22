function [] = shanks_line_search(alpha_scan, N_order)

    best_new_fom = 0;
    best_alpha = 0;
    for alpha = alpha_scan
        eShanks =...
          shanks_transformation(A,b, Jz, deeps_str, alpha, N_order, k0, Z0, omega);
      best_field = eShanks{N_order}{1};
      fom = log10((abs(best_field(eta==1)))^2)-ref;
      if(fom > best_new_fom)
         best_new_fom = fom;
         best_alpha = alpha;
      end
      shanks_alpha_history = [shanks_alpha_history, best_alpha];
      
      
end