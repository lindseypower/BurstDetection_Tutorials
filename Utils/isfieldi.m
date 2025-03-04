% Utility function written by Lindsey Power (March 2025) 
% Required to run functions in bg_SWM.m

function [val, cfgOut] = isfieldi(cfg, lst)
   for i=1:length(lst)
       field = lst{i};
       val = isfield(cfg, field); 
       cfgOut = cfg; 
   end
end
