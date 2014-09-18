classdef pv
    methods(Static)
        function classify
            % classify pv cells based on the correlation matrix
            cellType = 'PV';
            keys = fetch(aodi.CovMatrix & 'filter_id=11' & 'method=100' & (aodi.Labels & struct('celltype',cellType)));
            
            auc = [];
            for key = keys'
                [C0,vars] = fetch1(aodi.CovMatrix & setfield(key,'method',0),'corr_matrix','variances');
                if cove.avgCorr(C0)<0.2
                    [cellTypes,C,S,L,selection] = fetch1(aodi.CovMatrix*aodi.TrialTraces*aodi.ActiveCells & key, ...
                        'celltypes','corr_matrix','sparse','lowrank','selection');
                    pv = strcmpi(cellTypes(selection),cellType);
                    vars = mean(reshape(vars,[],size(vars,4)));
                    P = L*L'-S;
                    
                    pp = sqrt(diag(diag(-P)));
                    q2 = pp\S/pp;
                    q3 = pp\P/pp;
                    q4 = pp\L*L'/pp;
                    m1 = vars;
                    m2 = mean(abs(L'));
                    m3 = mean(q2);
                    m4 = mc(q3,pv);
                    % m = sum(q4,2);
                    m5 = m1-m4'/2000 - mean(S)/100;
                    auc(end+1) = roc(m5,pv);
                end
            end
            hist(auc)
            mean(auc)
        end        
    end
end



function m = mc(C,pv)
[U,D] = eigs(inv(C),3);
m = sum(abs(U*D),2);
end



function auc  = roc(m,truth)
bins = linspace(min(m),max(m),1000);
h1 = hist(m(truth),bins);
h2 = hist(m(~truth),bins);
h1 = cumsum(h1);
h2 = cumsum(h2);
h1 = h1/h1(end);
h2 = h2/h2(end);
plot(bins',[h1;h2]')
auc = (h1-h2)*gradient(h2)';
end
