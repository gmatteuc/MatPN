function [theta] = qdiff(s1,s2,quant)
theta = quantile(s1,quant)-quantile(s2,quant);
end