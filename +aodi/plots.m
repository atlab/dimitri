classdef plots
    
    properties(Constant)
        figPath = '~/dev/figures/interneurons'
    end
    
    methods(Static)
        
        
        function pvOri
            aodi.plots.pvOriPlot(true, false,'PV','noise correlations',          'ori-scor');
            aodi.plots.pvOriPlot(false,false,'PV','partial noise correlations',  'ori-pcor');
            aodi.plots.pvOriPlot(false,true, 'PV','avg conn strength (LV-GLASSO)',    'ori-sconn')
            aodi.plots.pvOriPlot(true, true, 'PV','avg conn strength (thresh. corrs)','ori-pconn')
        end        
        
        
        
        function pvOriPlot(doCorrs,doThresh,cellType,titl,fname)
            bins = [0 15 45];
            [A,B] = aodi.plots.pvOriHelper(doCorrs,doThresh,cellType,bins);
            bins2 = [bins(2:end) 90];
            
            style = {'v-' 'o-' 's-'};
            ti = {[cellType '-/-'],[cellType '-/+'] [cellType '+/+']};
            mx = max(reshape(A(:,1:2,:),[],1));
            if ~isempty(B)
                mn = min(0,min([reshape(A(:,1:2,:),[],1);reshape(B(:,1:2,:),[],1)]));
            else
                mn = min(0,min(reshape(A(:,1:2,:),[],1)));
            end
            s = 10^floor(log10(mx));
            yticks = (-10:10)*(floor(mx/s)*s);
            
            for i=1:3
                fig = Figure(1,'size',[60 50]);
                if i==3
                    A = nanmean(A,3);
                    B = nanmean(B,3);
                end
                if ~doThresh
                    plot(bins/2+bins2/2, squeeze(A(:,i,:)) ,style{i},...
                        'Color',[.4 .4 .4],'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor',[0 0 0],'MarkerSize',3);
                else
                    plot(bins/2+bins2/2, squeeze(A(:,i,:)) ,style{i},...
                        'Color',[.4 1 .4],'MarkerFaceColor',[.6 1 .6],'MarkerEdgeColor',[0 1 0],'MarkerSize',3);
                    hold on
                    plot(bins/2+bins2/2, squeeze(B(:,i,:)) ,style{i},...
                        'Color',[1 .4 .4],'MarkerFaceColor',[1 .6 .6],'MarkerEdgeColor',[1 0 0],'MarkerSize',3);
                    hold off
                end
                set(gca,'YTick',yticks,'YTickLabel',arrayfun(@(g) sprintf('%g',g), yticks, 'uni', false))
                set(gca,'XTick',union(bins,bins2))
                grid on
                ylabel(titl)
                xlabel '\Deltaori (degrees)'
                title(ti{i})
                ylim([mn mx]+[-(mn<0) 1]*0.05*mx)
                xlim([bins(1) bins2(end)])
                
                fig.cleanup
                set(gca,'Position',[0.22 0.20 0.70 0.70])
                fig.save(fullfile(aodi.plots.figPath,sprintf('%s-%u.eps', fname, i)))
            end
        end
        
        
        function [A,B] = pvOriHelper(doCorrs,doThresh,cellType,bins)
            A = [];
            B = [];
            minPerBin = 20;
            oriRel = pro(aodi.OriTuning,'stim_idx->s2','*');
            keys = fetch(aodi.CovMatrix & 'filter_id=11' & 'method=100' & ...
                (aodi.TrialTraces & 'evoked_bins>5') & (aodi.Labels & struct('celltype',cellType)) & oriRel);
            for key = keys'
                assert(1==count(aodi.CovMatrix & key), 'one matrix at a time please')
                assert(1==count(oriRel & key), 'must have one orientation map')
                assert(strcmp('lv-glasso', fetch1(aodi.CovMethod & key,'regularization')), ...
                    'lv-glasso regularizaton is required')
                [S,selection,cellTypes,L,pref,pval] = fetch1(...
                    aodi.CovMatrix*aodi.ActiveCells*aodi.TrialTraces*oriRel & key, ...
                    'sparse','selection','celltypes','lowrank','von_pref','von_p_value');
                p = sum(selection);
                pref = pref(selection)*180/pi;
                pval = pval(selection);
                C0 = fetch1(aodi.CovMatrix & setfield(key,'method',7),'corr_matrix');
                if cove.avgCorr(C0)>0.2
                    % limit to sites that have low average correlations
                    continue
                end
                disp(key)
                
                if doCorrs
                    % substitute S with thresholded correlations
                    if ~doThresh
                        S = C0;
                    else
                        [i,j] = meshgrid(1:p,1:p);
                        thresh = quantile(abs(C0(i(:)<j(:))),cove.sparsity(S));
                        S = C0.*(abs(C0)>thresh);
                    end
                end
                
                cellTypes = cellTypes(selection);
                pv = cellfun(@(c) strcmp(c,cellType), cellTypes);
                
                if any(pv)
                    fprintf('%d latent units (%d observed)\n', size(L'))
                    if doThresh
                        if ~doCorrs
                            d = (S-L*L')^.5;
                            S = -d\S/d;
                        end
                        S = S.*(1-eye(p));
                    else
                        % average correlations are stored in the excitatory
                        % connections
                        if ~doCorrs
                            S = corrcov(S-L*L');
                        end
                    end
                end
                if ~doCorrs && ~doThresh
                    % flip the sign back
                    S = -S;
                end
                
                [i,j] = ndgrid(1:size(S,1),1:size(S,2));
                ix = sub2ind(size(S),i(i<j),j(i<j));
                S = S(ix);
                D = ne7.rf.oriDiff(pref(i(ix)), pref(j(ix)));
                D(pval(i(ix))>0.05 | pval(j(ix))>0.05)=nan;
                pairTypes = 1+pv(i(ix))+pv(j(ix));
                
                A_  = [];
                B_  = [];
                c = inf;
                %                avgConn = mean(logical(S));
                for pairType = 1:3
                    ix = pairTypes==pairType;
                    b = sum(bsxfun(@ge,D(ix),bins),2);
                    if i<3
                        c = min(c,min(hist(b(b>0),1:length(bins))));
                    end
                    ix(ix) = b>0;
                    b = b(b>0);
                    if ~doThresh
                        A_ = cat(2,A_,accumarray(b, S(ix), [length(bins) 1], @mean));
                    else
                        ix = find(ix);
                        A_ = cat(2,A_,accumarray(b(S(ix)>0), S(ix(S(ix)>0)), [length(bins) 1], @mean, nan));
                        B_ = cat(2,B_,accumarray(b(S(ix)<0), S(ix(S(ix)<0)), [length(bins) 1], @mean, nan));
                    end
                end
                if c>=minPerBin
                    A = cat(3,A,A_);
                    B = cat(3,B,B_);
                end
            end
        end

        
        function pvDist
            aodi.plots.pvDistPlot(true, false,false,'PV','noise corrs',    'dist-scor');
            aodi.plots.pvDistPlot(false,false,false,'PV','partial corrs',  'dist-pcor');
            aodi.plots.pvDistPlot(true,true, false,'PV','% connectivity (sparse pcorrs)',    'dist-sconn')
            aodi.plots.pvDistPlot(false, true, false,'PV','% connectivity (thresh. corrs)','dist-pconn')
            aodi.plots.pvDistPlot(true,true, true, 'PV','neg/pos conn ratio (sparse pcorss)','dist-sratio')
            aodi.plots.pvDistPlot(false, true, true, 'PV','neg/pos conn ratio (thresh. corrs)','dist-pratio')
        end
        
        
        
        function pvDistPlot(doCorrs,doThresh,doRatio,cellType,titl,fname)
            bins = [0 30 75 150];
            [A,B] = aodi.plots.pvDistHelper(doCorrs,doThresh,cellType,bins);
            bins2 = [bins(2:end) 200];
            
            style = {'v-','o-','s-'};
            ti = {[cellType '-/-'],[cellType '-/+'],[cellType '+/+']};
            mx = max(reshape(A(:,1:2,:),[],1));
            mn = min(0,min(reshape(A(:,1:2,:),[],1)));
            s = 10^floor(log10(mx));
            yticks = [0 floor(mx/s)*s];
            
            for i=1:3
                if i==3
                    A = nanmean(A,3);
                    B = nanmean(B,3);
                end
                fig = Figure(1,'size',[60 50]);
                if ~doThresh
                    plot(bins/2+bins2/2, squeeze(A(:,i,:)) ,style{i},...
                        'Color',[.4 .4 .4],'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor',[0 0 0],'MarkerSize',3);
                else
                    if doRatio
                        plot(bins/2+bins2/2, (squeeze(A(:,i,:))+squeeze(B(:,i,:))).\squeeze(B(:,i,:)) ,style{i},...
                            'Color',[.4 .4 .4],'MarkerFaceColor',[.6 .6 .6],'MarkerEdgeColor',[0 0 0],'MarkerSize',3);
                        mx = 0.9;
                        mn = 0;
                        yticks = [0 .5];
                    else
                        plot(bins/2+bins2/2, squeeze(A(:,i,:)) ,style{i},...
                            'Color',[.4 1 .4],'MarkerFaceColor',[.6 1 .6],'MarkerEdgeColor',[0 1 0],'MarkerSize',3);
                        hold on
                        plot(bins/2+bins2/2, squeeze(B(:,i,:)) ,style{i},...
                            'Color',[1 .4 .4],'MarkerFaceColor',[1 .6 .6],'MarkerEdgeColor',[1 0 0],'MarkerSize',3);
                        hold off  
                    end
                end
                set(gca,'YTick',yticks,'YTickLabel',arrayfun(@(g) sprintf('%g',g), yticks, 'uni', false))
                set(gca,'XTick',bins)
                grid on
                ylabel(titl)
                xlabel 'Distance (\mum)'
                title(ti{i})
                ylim([mn mx]+[-(mn<0) 1]*0.05*mx)
                
                
                fig.cleanup
                set(gca,'Position',[0.22 0.20 0.70 0.70])
                fig.save(fullfile(aodi.plots.figPath,sprintf('%s-%u.eps', fname, i)))
            end
        end
        
        
        function [A,B] = pvDistHelper(doCorrs,doThresh,cellType,bins)
            A = [];
            B = [];
            minPerBin = 10;
            keys = fetch(aodi.CovMatrix & 'filter_id=11' & 'method=100' & ...
                (aodi.TrialTraces & 'evoked_bins>5') & (aodi.Labels & struct('celltype',cellType)));
            for key = keys'
                assert(1==count(aodi.CovMatrix & key), 'one matrix at a time please')
                assert(strcmp('lv-glasso', fetch1(aodi.CovMethod & key,'regularization')), ...
                    'lv-glasso regularizaton is required')
                
                [S,selection,cellTypes,L,xyz] = fetch1(...
                    aodi.CovMatrix*aodi.ActiveCells*aodi.TrialTraces & key, ...
                    'sparse','selection','celltypes','lowrank','cell_xyz');
                p = sum(selection);
                C0 = fetch1(aodi.CovMatrix & setfield(key,'method',7),'corr_matrix');
                if cove.avgCorr(C0)>0.08
                    continue
                end
                
                if doCorrs
                    % substitute S with thresholded correlations
                    if ~doThresh
                        S = C0;
                    else
                        [i,j] = meshgrid(1:p,1:p);
                        thresh = quantile(abs(C0(i(:)<j(:))),cove.sparsity(S));
                        S = C0.*(abs(C0)>thresh);
                    end
                end
                
                cellTypes = cellTypes(selection);
                pv = cellfun(@(c) strcmp(c,cellType), cellTypes);
                
                if any(pv)
                    fprintf('%d latent units (%d observed)\n', size(L'))
                    if doThresh
                        if ~doCorrs
                            S = -S;
                        end
                        S = S.*(1-eye(p));
                    else
                        % average correlations are stored in the excitatory
                        % connections
                        if ~doCorrs
                            S = corrcov(S-L*L');
                        end
                    end
                end
                if ~doCorrs && ~doThresh
                    % flip the sign back
                    S = -S;
                end
                
                xyz = xyz(selection,:);
                [i,j] = ndgrid(1:size(S,1),1:size(S,2));
                ix = sub2ind(size(S),i(i<j),j(i<j));
                dz = sqrt(sum((xyz(i(ix),3)-xyz(j(ix),3)).^2,2));  % limit to the same layer
                ix = ix(dz<200);
                S = S(ix);
                D = sqrt(sum((xyz(i(ix),1:3)-xyz(j(ix),1:3)).^2,2));
                pairTypes = 1+pv(i(ix))+pv(j(ix));
                
                A_  = [];
                B_  = [];
                c = inf;
%                avgConn = mean(logical(S));
                for pairType = 1:3
                    ix = pairTypes==pairType;
                    b = sum(bsxfun(@ge,D(ix),bins),2);
                    if i<3
                        c = min(c,min(hist(b,1:length(bins))));
                    end
                    if ~doThresh
                        A_ = cat(2,A_,accumarray(b, S(ix), [length(bins) 1], @mean));
                    else
                        A_ = cat(2,A_,accumarray(b, S(ix)>0, [length(bins) 1], @mean));
                        B_ = cat(2,B_,accumarray(b, S(ix)<0, [length(bins) 1], @mean));
                    end
                end
                if c>=minPerBin
                    A = cat(3,A,A_);
                    B = cat(3,B,B_);
                end
            end
        end
        
        
        
        function qual
            
            % site selection
            n1 = 15;   % of 31
            n2 = 4;   % of 24
            
            % getting data from
            key1 = fetch(pop.AodBinnedTraces, sprintf('ORDER BY nneurons DESC LIMIT 1 OFFSET %d',n1-1));
            key2 = fetch(aodi.TrialTraces & 'filter_id=2', sprintf('ORDER BY ntraces DESC LIMIT 1 OFFSET %d',n2-1));

            X1 = fetchn(aod.Traces & key1, 'trace');
            fprintf('Number of cells %d\n',length(X1))
            times1 = getTimes(aod.Traces & key1)/1000;
            times1 = times1-times1(1);
            dt1 = mean(diff(times1));
            
            [X2,Y2] = fetchn(vis2p.MaskTracesRaw & key2, 'calcium_trace','red_trace');
            times2 = fetch1(vis2p.MaskGroupRaw & key2, 'frame_timestamps')/1000;
            times2 = times2-times2(1);
            dt2 = mean(diff(times2));
            fprintf('Number of cells %d\n',length(X2))
            
            % plot random pair
            w = 0.1;  % time width of the smoothing kernel 
            
            for i=1:100
                h(1) = subplot(211);
                k = hamming(2*round(w/dt1)+1); k=k/sum(k);
                plot(times1, ne7.dsp.convmirr(double(X1{randi(end)}),k))

                h(2) = subplot(212);
                k = hamming(2*round(w/dt2)+1); k=k/sum(k);
                plot(times2, ne7.dsp.convmirr(double(X2{randi(end)}),k))
                
                linkaxes(h,'x')
                keyboard
            end
        end
        
        
        
        function sigCorr
            keys = fetch(aodi.CovMatrix & 'filter_id=11' & 'method=100' & ...
                (aodi.Labels & 'celltype="PV"'));
            label = 'PV';
            cor = zeros(0,3);
            for key = keys'
                sel = fetch1(aodi.ActiveCells & key, 'selection');
                C0 = fetch1(aodi.CovMatrix & setfield(key,'method',7), 'corr_matrix');
                if cove.avgCorr(C0) > 0.05
                    continue
                end
%                C0  = fetch1(aodi.CovMatrix & key, 'corr_matrix');
%                C0 = -corrcov(inv(C0));
                Sig = fetch1(aodi.SigCorrs & key, 'sig_corrs');
                Sig = Sig(sel,sel);
                
                cellTypes = fetch1(aodi.TrialTraces & key, 'celltypes');               
                
                ix = strcmp(cellTypes, label);
                ix = ix(sel);
                
                if any(ix)
                    q = zeros(size(C0));
                    q(ix,~ix) = 1;   % +/-
                    q(ix, ix) = 2;   % +/+

                    p = sum(sel);
                    [i,j] = meshgrid(1:p,1:p);
                    sig = Sig(i<j);
                    c = C0(i<j);
                    q = q(i<j);
                    
                    fig = Figure(1,'size',[80 80]);
                    scatter(sig(q==0),c(q==0),2,'r.')
                    hold on
                    scatter(sig(q==1),c(q==1),4,'gx')                    
                    scatter(sig(q==2),c(q==2),6,'ko','filled')
                     ix=q==0; r0 = regress(c(ix),[sig(ix) ones(sum(ix),1)]);
                     ix=q==1; r1 = regress(c(ix),[sig(ix) ones(sum(ix),1)]);
                     ix=q==2; r2 = regress(c(ix),[sig(ix) ones(sum(ix),1)]);
%                    ix=q==0; r0 = corr(c(ix),sig(ix));
%                    ix=q==1; r1 = corr(c(ix),sig(ix));
%                    ix=q==2; r2 = corr(c(ix),sig(ix));

                    cor(end+1,:) = [r0(1) r1(1) r2(1)]; %#ok<AGROW>
                    hold off
                    legend -/-  -/+  +/+
                    legend location southeast
                    h = refline;
                    set(h(3),'Color','r')
                    set(h(2),'Color','g')
                    set(h(1),'Color','k')
                    xlabel 'signal correlations'
                    ylabel 'noise corrs'
                    ylim([-.2 .2])
%                    drawnow
%                     fig.cleanup
%                     fname = fullfile(aodi.plots.figPath,...
%                         ['sigcorr' sprintf('%05d-%s--%02d.eps',key.mouse_id,key.exp_date,key.scan_idx)]);
%                     fig.save(fname)
                end
            end
            fig = Figure(1,'size',[70 70]);
            plot(1:3,cor,'^-','Color',[.4 .4 .7],'MarkerFaceColor',[.1 .1 .7],'MarkerEdgeColor',[0 0 .4])
            set(gca,'XTick',1:3,'XTickLabel',{'PV-/PV-','PV-/PV+','PV+/PV+'})
            ylabel('R of noise corr vs sig corrs')
            xlim([0.5 3.5])
            fig.cleanup
            fig.save('~/Desktop/csignoise.eps')
        end
        
        
        
        function pvOSI
            label = 'PV';
            keys = fetch(aodi.OriTuning & 'filter_id=4');
            for key = keys'
                [cellTypes,r2,p,pref,base,amp1,amp2,sharp] = fetch1(aodi.TrialTraces*aodi.OriTuning & key, ...
                    'celltypes', 'von_r2', 'von_p_value', 'von_pref', 'von_base', 'von_amp1', 'von_amp2', 'von_sharp');
                pref = base + amp1.*exp(sharp-1)+amp2.*exp(-sharp-1);
                ortho = base + (amp1+amp2).*exp(-1);
                osi = (pref-ortho)./(pref + ortho);
                pv = strcmp(cellTypes,label);                
            end
        end
        


        function oriDiff
            label = 'PV';
            keys = fetch(pro(aodi.OriTuning,'stim_idx->si')*(aodi.CovMatrix & 'method=100') & 'filter_id=11' & (aodi.Labels & struct('celltype',label)));
            pp = [];  pn = [];  nn = [];
            for key = keys'
                [cellTypes, pref, pval, r2, C0] = fetch1(...
                    pro(aodi.OriTuning,'stim_idx->si','*') * aodi.TrialTraces * aodi.CovMatrix & key, ...
                    'celltypes','von_pref','von_p_value','von_r2', 'corr_matrix');
                if cove.avgCorr(C0)<0.08
                    pv = strcmpi(cellTypes,label);
                    tuned = pval < 0.05;
                    fprintf('Tuned non-PV %2.1f%%  PV %2.1f%%\n', mean(tuned(~pv))*100, mean(tuned(pv))*100)
                    pv = pv(tuned);
                    r2 = r2(tuned);
                    pref = pref(tuned)*180/pi;

                    if sum(pv)
                        [i1,i2] = allPairs(find(pv),find(pv));
                        pp(end+1) = mean(ne7.rf.oriDiff(pref(i1),pref(i2)));

                        [i1,i2] = allPairs(find(pv),find(~pv));
                        pn(end+1) = mean(ne7.rf.oriDiff(pref(i1),pref(i2)));

                        [i1,i2] = allPairs(find(~pv),find(~pv));
                        nn(end+1) = mean(ne7.rf.oriDiff(pref(i1),pref(i2)));
                    end
                end
            end
        end
            
        
        function pvConn
            aodi.plots.pvConnPlot(true, false,false,'PV','noise correlations',          'PV-scorr')
            aodi.plots.pvConnPlot(false,false,false,'PV','partial corrs',               'PV-pcorr')
            aodi.plots.pvConnPlot(false,true, false,'PV','Connectivity (sparse partial)',    'PV-pconn')
            aodi.plots.pvConnPlot(false,true, true,'PV','Neg/pos ratio (sparse partial)',    'PV-pratio')
            aodi.plots.pvConnPlot(true, true, false,'PV','Connectivity (thresh. corrs)','PV-sconn')
            aodi.plots.pvConnPlot(true, true, true,'PV','Neg/pos ratio (thresh. corrs)','PV-sratio')

        end
        
        function pvConnPlot(doCorrs,doThresh,doRatio,cellType,titl,filepath)
            [enn,epn,epp,inn,ipn,ipp] = aodi.plots.pvConnHelper(doCorrs,doThresh,cellType);
            filepath = fullfile(aodi.plots.figPath, filepath);
            
            fig = Figure(1,'size',[60 80]);
            if doRatio
                plot(1:3,[inn;ipn;ipp]./([enn;epn;epp]+[inn;ipn;ipp]),'^-','Color',[.4 .4 .4],'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[0 0 0])
            else
                if isempty(inn)
                    plot(1:3,[enn;epn;epp],'^-','Color',[.4 .4 .4],'MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor',[0 0 0])
                else
                    plot(1:3,[enn;epn;epp],'^-','Color',[.4 .7 .4],'MarkerFaceColor',[.7 1 .7],'MarkerEdgeColor',[0 .4 0])
                    hold on
                    plot(1:3,[inn;ipn;ipp],'o-','Color',[.7 .4 .4],'MarkerFaceColor',[1 .7 .7],'MarkerEdgeColor',[.4 0 0])
                    hold off
                end
            end
            xlim([0.5 3.5])
            box off
            set(gca,'XTick',1:3, 'XTickLabel', {'-/-','-/+','+/+'},'YGrid','on')
            toStr = @(s) sprintf('%0.3g',s);
            set(gca,'YTickLabel', arrayfun(toStr,get(gca,'YTick'),'uni',false))
            title(titl)
            fig.cleanup
            fig.save(sprintf('%s.eps',filepath))
        end
        
        
        function [enn,epn,epp,inn,ipn,ipp] = pvConnHelper(doCorrs,doThresh,cellType)    
            keys = fetch(aodi.CovMatrix & 'filter_id=11' & 'method=100' & ...
                (aodi.TrialTraces & 'evoked_bins>5') & (aodi.Labels & struct('celltype',cellType)));
            epp = [];
            epn = [];
            enn = [];
            ipp = [];
            ipn = [];
            inn = [];
            for key = keys'
                assert(1==count(aodi.CovMatrix & key), 'one matrix at a time please')
                assert(strcmp('lv-glasso', fetch1(aodi.CovMethod & key,'regularization')), ...
                    'lv-glasso regularizaton is required')
                
                [S,selection,cellTypes,L] = fetch1(...
                    aodi.CovMatrix*aodi.ActiveCells*aodi.TrialTraces & key, ...
                    'sparse','selection','celltypes','lowrank');
                p = sum(selection);
                C0 = fetch1(aodi.CovMatrix & setfield(key,'method',7),'corr_matrix');
                if cove.avgCorr(C0)>0.05
                    continue
                end
                    
                if doCorrs
                    % substitute S with thresholded correlations
                    if ~doThresh
                        S = C0;
                    else
                        [i,j] = meshgrid(1:p,1:p);
                        thresh = quantile(abs(C0(i(:)<j(:))),cove.sparsity(S));
                        S = C0.*(abs(C0)>thresh);
                    end
                end
                
                cellTypes = cellTypes(selection);
                neurons = cellfun(@(c) strcmp(c,'neuron'), cellTypes);
                pv = cellfun(@(c) strcmp(c,cellType), cellTypes);
                
                if any(pv)
                    fprintf('%d latent units (%d observed)\n', size(L'))
                    if doThresh
                        if ~doCorrs
                            S = -S;
                        end
                        S = S.*(1-eye(p));
                        exc = S>0;
                        inh = S<0;
                        epp(end+1) = cove.avgCorr(exc(pv,pv)); %#ok<AGROW>
                        epn(end+1) = mean(mean(exc(neurons,pv)));%#ok<AGROW>
                        enn(end+1) = cove.avgCorr(exc(neurons,neurons));%#ok<AGROW>

                        ipp(end+1) = cove.avgCorr(inh(pv,pv));%#ok<AGROW>
                        ipn(end+1) = mean(mean(inh(neurons,pv)));%#ok<AGROW>
                        inn(end+1) = cove.avgCorr(inh(neurons,neurons));%#ok<AGROW>
                    else
                        % average correlations are stored in the excitatory
                        % connections
                        if ~doCorrs
                            S = corrcov(S-L*L');
                        end
                        epp(end+1) = cove.avgCorr(S(pv,pv));%#ok<AGROW>
                        epn(end+1) = mean(mean(S(neurons,pv)));%#ok<AGROW>
                        enn(end+1) = cove.avgCorr(S(neurons,neurons));%#ok<AGROW>
                    end
                end
            end
            if ~doCorrs && ~doThresh
                % flip the sign back
                epp = -epp;
                epn = -epn;
                enn = -enn;
            end
            

        end
        
        
        
        function pvNetwork
            rel = aodi.CovMatrix & 'method=100' & 'filter_id=11' ...
                & pro(aodi.OriTuning,'stim_idx->s2')...
                & (aodi.TrialTraces & 'evoked_bins>5')...
                & (aodi.Labels & 'celltype="pv"');
            usedKeys = [];
            for key = rel.fetch'
                C0 = fetch1(aodi.CovMatrix & setfield(key,'method',7),'corr_matrix');
                if cove.avgCorr(C0)<0.1
                    usedKeys = [usedKeys key];
                    aodi.plots.network(key,'pv')
                end
            end
        end
            
            
        
        function network(key,cellType)
            % keys = fetch(aodi.CovMatrix & 'method=100' & 'filter_id=11' & pro(aodi.OriTuning,'stim_idx->s2') & (aodi.TrialTraces & 'evoked_bins>5'))
            assert(1==count(aodi.CovMatrix & key), 'one matrix at a time please')
            assert(strcmp('lv-glasso', fetch1(aodi.CovMethod & key,'regularization')), ...
                'lv-glasso regularizaton is required')
            assert(1==count(aodi.OriTuning & rmfield(key,'stim_idx')))
            
            clf
            
            doInteractions = true;  % true=plot interactions, false=only cells
            
            % figure 4-G,H,I
            alpha = 0.05;  % tuning signficance levels
            xoffset = 0;
            yoffset = -10;
            zoffset = -10;
            
            zticks = 100:50:300;
            xticks = -200:50:200;
            yticks = -200:50:200;
            zm = .60;
            panx = 1.3;
            pany = 0;
            alphaMultiplier = 3.0;
            
            if nargin<2
                cellType = '';
            end
            fname = fullfile(aodi.plots.figPath,...
                ['network' sprintf('%05d-%s--%02d%s',key.mouse_id,key.exp_date,key.scan_idx,cellType)]);
            
            if doInteractions
                paperSize = [12 11];
                lineWidth = 0.5;
            end
            
            % get cell positions, tuning, and sparse interactions
            key = fetch(aodi.CovMatrix & key);
            assert(isscalar(key))
            [xyz,cellTypes] = fetch1(aodi.TrialTraces & key,'cell_xyz','celltypes');
            selection = fetch1(aodi.ActiveCells & key, 'selection');
            [ori,pval] = fetch1(aodi.OriTuning & rmfield(key,'stim_idx'), 'von_pref','von_p_value');
            [S,L] = fetch1(aodi.CovMatrix & key, 'sparse', 'lowrank');
            S = -corrcov(S);  % convert to partial correlations
            
%             % thresholded correlations
%             C0= fetch1(aodi.CovMatrix & setfield(key, 'method', 0), 'corr_matrix'); %#ok<SFLD>
%             sparsity = fetch1(aodi.CovMatrix & key, 'sparsity');
%             p = size(C0,1);
%             [i,j] = ndgrid(1:p,1:p);
%             C0 = C0.*(abs(C0)>quantile(abs(C0(i<j)), sparsity));
%             
%             % report everything
%             fprintf('Sparsity = %2.1f%%\n', 100*sparsity)
%             fprintf('Overlap = %2.1f%%\n',  100*sum(C0(i<j) & S(i<j))/sum(~~S(i<j)))
%             fprintf('Latent = %d\n', size(L,2));
%             fprintf('Negative interactions = %2.1f%%\n', 100*sum(S(i<j)<0)/sum(~~S(i<j)))
            
            zmin = min(xyz(:,3));
            zmax = max(xyz(:,3));
            if nargin>=2
                % limit to cell type 
                alphaMultiplier = alphaMultiplier*3;
                ixx = strcmpi(cellTypes,cellType);
                assert(any(ixx))
                S = S(ixx,ixx);
                selection = selection(ixx);
                xyz = xyz(ixx,:);
                ori = ori(ixx);
                pval = pval(ixx);                
            end
            
            x = xyz(:,1);
            y = xyz(:,2);
            z = xyz(:,3);
            
            % plot balls
            hue = mod(ori(:)/pi,1);
            sat = pval(:)<alpha;            
            val = (1-.2*(pval(:)>=alpha)).*selection(:);
            color = hsv2rgb([hue sat val]);
            clear hue sat val
            
            fragIdx = true(size(x)); 
          
            scatter3sph(x(fragIdx),y(fragIdx),z(fragIdx),'siz',3,'col',color(fragIdx,:))
            light('Position',[0.5 0.5 1],'Style','infinit','Color',[1 1 1])
            lighting gouraud
            axis vis3d
            axis equal
            set(gca,'ZDir','reverse')
            camproj perspective
            grid on
            
            % show interactions
            x = x(selection);
            y = y(selection);
            z = z(selection);
            
            if doInteractions
                [i,j] = ndgrid(1:size(S,1),1:size(S,2));
                % positive interactions
                ix = find(j(:)>i(:) & S(:)>0);
                for ix = ix'
                    lineAlpha = min(1,(alphaMultiplier*2.0*abs(S(ix))));
                    patch([x(i(ix)) x(j(ix))],[y(i(ix)) y(j(ix))], [z(i(ix)) z(j(ix))], [0 0], ...
                        'EdgeColor',[0 1 0]/2,'EdgeAlpha',lineAlpha,'LineWidth',lineWidth,'FaceColor','none')
                end
                
                % negative interactions
                ix = find(j(:)>i(:) & S(:)<0);
                for ix = ix'
                    lineAlpha = min(1,(alphaMultiplier*4.0*abs(S(ix))));
                    patch([x(i(ix)) x(j(ix))],[y(i(ix)) y(j(ix))], [z(i(ix)) z(j(ix))], [0 0], ...
                        'EdgeColor',[1 0 0]/2,'EdgeAlpha',lineAlpha,'LineWidth',lineWidth,'FaceColor','none')
                end
            end
            view(25-90, 65)
            zoom(zm)
            xlim([-100 100])
            ylim([-100 100])
            zlim([zmin zmax]+[-1 1]*0.05*(zmin-zmax))
            camPos = get(gca,'CameraPosition');
            set(gca,'CameraPosition',camPos*0.5)
            set(gca,'ZTick',zticks)
            set(gca,'XTick',xticks,'YTick',yticks)
            campan(panx,pany)
            set(gca,'fontsize',8,'linewidth',0.25,'TickLength',get(gca,'TickLength')*0.75)
            
            set(gcf,'PaperUnits','centimeters','PaperSize',paperSize,'PaperPosition',[0 0 paperSize])
            print('-dpdf','-r800',fname)
        end
        
    end
end


function [ix1, ix2] = allPairs(ix1,ix2)
[i,j] = meshgrid(ix1,ix2);
ix1 = i(i<j);
ix2 = j(i<j);
end