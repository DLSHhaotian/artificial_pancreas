function [noise_cc,noise_v] = sensorNoise(noise_cc_prepre,noise_cc_pre,noise_v_prepre,noise_v_pre,wcc,w)
% Syntax: [noise_cc,noise_v] = sensorNoise(noise_cc_prepre,noise_cc_pre,noise_v_prepre,noise_v_pre,wcc,w)
%         noise_cc_prepre: last of last noise cc
%         noise_cc_pre: last noise cc
%         noise_v_prepre: last of last noise v
%         noise_v_pre: last noise v
%         wcc: added term of cc
%         w: added term of v

noise_cc=1.23*noise_cc_pre-0.3995*noise_cc_prepre+wcc;
noise_v=1.013*noise_v_pre-0.2135*noise_v_prepre+w;
end