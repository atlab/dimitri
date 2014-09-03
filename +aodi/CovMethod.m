%{
aodi.CovMethod (lookup) # covariance matrix estimation methods
method          : tinyint                # covariance estimation method
---
include_spont               : tinyint                       # 1=yes, 0=no
condition                   : tinyint                       # 0=all conditions, 1=first conditoin only, 2=second condition only
regularization              : enum('sample','diag','factor','glasso','lv-glasso') # 
hyperparam_space            : longblob                      # arrays of hyperparameter values
binwise=1                   : tinyint                       # 1= binwise, 0= trial-wise
%}

classdef CovMethod < dj.Relvar
    methods
        function fill(self)
            
            self.inserti({
                
            0   1 0  'sample'  {} 1
            5   1 0  'sample'  {} 0

            %             1   0 0  'sample'  {} 1
%             2   0 1  'sample'  {} 1
%             3   0 2  'sample'  {} 1
            
%            10  1 0  'diag'      {exp(-6:0.05:-1) exp(-4:0.05:0)}   1
%            30  1 0  'factor'    {exp(-6:0.05:-1) 0:70}   1
%            40  1 0  'factor'    {0:70 exp(-6:0.05:-1)}   1
%            80  1 0  'glasso'    {exp(-5:0.05:-1)}     1
%            90  1 0  'lv-glasso' {exp(-7:.05:-1) exp(-5:.05:-1)} 1
            
            100  1 0  'lv-glasso' {exp(-7:.05:-1.5) exp(-5:.05:-1.5)} 1
            105  1 0  'lv-glasso' {exp(-7:.05:-1.5) exp(-5:.05:-1.5)} 0
            })
        end
    end
end
