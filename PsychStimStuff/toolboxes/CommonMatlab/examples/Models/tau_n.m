function out = tau_n(V,theta_t_na,sigma_t_na,theta_t_nb,sigma_t_nb)
% function out = tau_n(V,theta_t_na,sigma_t_na,theta_t_nb,sigma_t_nb)
% tau_n=(0.087+11.4*gammaf(V,theta_tna,sigma_tna))*(0.087+11.4*gammaf(V,the
% ta_tnb,sigma_tnb)); 
out = (0.087+11.4*gammaf(V,theta_t_na,sigma_t_na))*(0.087+11.4*gammaf(V,theta_t_nb,sigma_t_nb));
