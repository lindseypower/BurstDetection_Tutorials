function [val, cfgOut] = isfieldi(cfg, lst)
   for i=1:length(lst)
       field = lst{i};
       val = isfield(cfg, field); 
       cfgOut = cfg; 
   end
end
