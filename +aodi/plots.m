classdef plots
    
    properties(Constant)
        figPath = '~/dev/figures/interneurons'
    end
    
    methods(Static)
        
        
        
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
            keys = fetch(aodi.CovMatrix & 'filter_id=4' & 'method=1' & ...
                (aodi.TrialTraces & 'evoked_bins>5' & 'min_trials>200'));
            label = 'PV';
            cor = zeros(0,3);
            for key = keys'
                sel = fetch1(aodi.ActiveCells & key, 'selection');
                C0 = fetch1(aodi.CovMatrix & key,'corr_matrix');
                S = fetch1(aodi.SigCorrs & key, 'sig_corrs');
                S = S(sel,sel);
                cellTypes = fetch1(aodi.TrialTraces & key, 'celltypes');               
                
                ix = strcmp(cellTypes, label);
                ix = ix(sel);
                
                if any(ix)
                    q = zeros(size(C0));
                    q(ix,~ix) = 1;   % +/-
                    q(ix, ix) = 2;   % +/+

                    p = sum(sel);
                    [i,j] = meshgrid(1:p,1:p);
                    s = S(i<j);
                    c = C0(i<j);
                    q = q(i<j);
                    
                    scatter(s(q==0),c(q==0),'r.')
                    hold on
                    scatter(s(q==1),c(q==1),'bx')                    
                    scatter(s(q==2),c(q==2),'ko','filled')
                    r0 = corr(s(q==0),c(q==0));
                    r1 = corr(s(q==1),c(q==1));
                    r2 = corr(s(q==2),c(q==2));
                    cor(end+1,:) = [r0 r1 r2];
                    hold off
                    h = refline;
                    set(h(3),'Color','r')
                    set(h(2),'Color','b')
                    set(h(1),'Color','k')
                    drawnow
                end
            end
            fig = Figure(1,'size',[60 60]);
            boxplot(cor)
            set(gca,'XTick',1:3,'XTickLabel',{'PV-/PV-','PV-/PV+','PV+/PV+'})
            ylabel('corr of sig vs noise corr')
            fig.cleanup
            fig.save('~/Desktop/signoise.eps')
        end
        
        
        
        function pvOSI()
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
            label = 'red';
            keys = fetch(aodi.OriTuning & 'filter_id=4');
            pp = [];  pn = [];  nn = [];
            for key = keys'
                [cellTypes, pref, pval, r2] = fetch1(...
                    aodi.OriTuning * aodi.TrialTraces & key, ...
                    'celltypes','von_pref','von_p_value','von_r2');
                pv = cellfun(@(c) strcmp(c,label), cellTypes);
                tuned = pval < 0.01 & r2 > 0.01;
                pv = pv(tuned);
                r2 = r2(tuned);
                pref = pref(tuned)*180/pi;
                disp(pref(pv))
                
                if sum(pv)>=2
                    [i1,i2] = allPairs(find(pv),find(pv));
                    pp(end+1) = mean(ne7.rf.oriDiff(pref(i1),pref(i2)));
                    
                    [i1,i2] = allPairs(find(pv),find(~pv));
                    pn(end+1) = mean(ne7.rf.oriDiff(pref(i1),pref(i2)));
                    
                    [i1,i2] = allPairs(find(~pv),find(~pv));
                    nn(end+1) = mean(ne7.rf.oriDiff(pref(i1),pref(i2)));
                end
            end
        end
            
        
        function pvConn
            aodi.plots.pvConnPlot(true, false,'PV','noise correlations',          'PV-scorr')
            aodi.plots.pvConnPlot(false,false,'PV','partial corrs',               'PV-pcorr')
            aodi.plots.pvConnPlot(false,true, 'PV','Connectivity (LV-GLASSO)',    'PV-pconn')
            aodi.plots.pvConnPlot(true, true, 'PV','Connectivity (thresh. corrs)','PV-sconn')
        end
        
        function pvConnPlot(doCorrs,doThresh,cellType,titl,filepath)
            [enn,epn,epp,inn,ipn,ipp] = aodi.plots.pvConnHelper(doCorrs,doThresh,cellType);
            filepath = fullfile(aodi.plots.figPath, filepath);
            
            for iPlot=1:3
                switch iPlot
                    case 1
                        xylabels ={'-/-' '-/+'};
                        x1 = enn;
                        x2 = inn;
                        y1 = epn;
                        y2 = ipn;
                        
                    case 2
                        xylabels ={'-/-' '+/+'};
                        x1 = enn;
                        x2 = inn;
                        y1 = epp;
                        y2 = ipp;
                        
                    case 3
                        xylabels ={'-/+' '+/+'};
                        x1 = epn;
                        x2 = ipn;
                        y1 = epp;
                        y2 = ipp;
                end
                
                fig = Figure(1,'size',[80 80]);
                scatter(x1,y1,'ko','filled')
                if ~isempty(x2)
                    hold on
                    scatter(x2,y2,'r^','filled')
                    legend positive negative
                    legend Location SouthEast
                    legend boxoff
                    for i=1:length(enn)
                        line([x1(i) x2(i)],[y1(i) y2(i)],'Color',[1 1 1]*0.7)
                    end
                    hold off
                end
                mx = max(max(x1), max(y1));
                ticks = 10^floor(log10(mx));
                ticks = 0:ticks:10*ticks;
                
                set(gca,'XTick',ticks,'YTick',ticks,'XTickLabel',ticks,'YTickLabel',ticks)
                axis(1.1*[0 1 0 1]*mx)
                set(refline(1),'Color',[1 1 1]*0.7,'LineWidth',.5,'LineStyle',':')
                xlabel(xylabels{1})
                ylabel(xylabels{2})
                grid on
                title(titl)
                fig.cleanup
                fig.save(sprintf('%s-scatter%d.eps',filepath,iPlot))
            end
            
            fig = Figure(1,'size',[83 35]);
            h = boxplot([enn' epn' epp'],'colors','k','labels',{'-/-' '-/+' '+/+'},'orientation','horizontal');
            set(h(1:2,:),'LineStyle','-','LineWidth',.25)
            set(h(7,:),'MarkerEdgeColor','k')
            xlim(min(xlim,xlim.*[0 1]))
            xlabel(titl)
            str = @(n) sprintf('%0.1g',n);
            set(gca,'XTickLabel', arrayfun(str, get(gca,'XTick'), 'uni', false))
            fig.cleanup
            fig.save(sprintf('%s-bars-pos.eps',filepath))
            
            if ~isempty(inn)
                fig = Figure(1,'size',[83 35]);
                h = boxplot([inn' ipn' ipp'],'colors','k','labels',{'-/-' '-/+' '+/+'},'orientation','horizontal');
                set(h(1:2,:),'LineStyle','-','LineWidth',.25)
                set(h(7,:),'MarkerEdgeColor','k')
                xlim(min(xlim,xlim.*[0 1]))
                xlabel([titl '(neg)'])
                fig.cleanup
                fig.save(sprintf('%s-bars-neg.eps',filepath))
            end

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
            counts = zeros(0,2);
            for key = keys'
                assert(1==count(aodi.CovMatrix & key), 'one matrix at a time please')
                assert(strcmp('lv-glasso', fetch1(aodi.CovMethod & key,'regularization')), ...
                    'lv-glasso regularizaton is required')
                
                [S,selection,cellTypes,L] = fetch1(...
                    aodi.CovMatrix*aodi.ActiveCells*aodi.TrialTraces & key, ...
                    'sparse','selection','celltypes','lowrank');
                p = sum(selection);
                C0 = fetch1(aodi.CovMatrix & setfield(key,'method',0),'corr_matrix');
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
                    counts(end+1,:) = [sum(pv), length(pv)];  %#ok<AGROW>
                    if doThresh
                        if ~doCorrs
                            S = -S;
                        end
                        S = S.*(1-eye(p));
                        exc = S>0;
                        inh = S<0;
                        epp(end+1) = mean(mean(exc(pv,pv))); %#ok<AGROW>
                        epn(end+1) = mean(mean(exc(neurons,pv)));%#ok<AGROW>
                        enn(end+1) = mean(mean(exc(neurons,neurons)));%#ok<AGROW>

                        ipp(end+1) = mean(mean(inh(pv,pv)));%#ok<AGROW>
                        ipn(end+1) = mean(mean(inh(neurons,pv)));%#ok<AGROW>
                        inn(end+1) = mean(mean(inh(neurons,neurons)));%#ok<AGROW>
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
        
        
        function network(key)
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
            
            fname = fullfile(aodi.plots.figPath,...
                ['network' sprintf('%05d-%s--%02d',key.mouse_id,key.exp_date,key.scan_idx)]);
            
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
            
            % thresholded correlations
            C0= fetch1(aodi.CovMatrix & setfield(key, 'method', 0), 'corr_matrix'); %#ok<SFLD>
            sparsity = fetch1(aodi.CovMatrix & key, 'sparsity');
            p = size(C0,1);
            [i,j] = ndgrid(1:p,1:p);
            C0 = C0.*(abs(C0)>quantile(abs(C0(i<j)), sparsity));
            
            % report everything
            fprintf('Sparsity = %2.1f%%\n', 100*sparsity)
            fprintf('Overlap = %2.1f%%\n',  100*sum(C0(i<j) & S(i<j))/sum(~~S(i<j)))
            fprintf('Latent = %d\n', size(L,2));
            fprintf('Negative interactions = %2.1f%%\n', 100*sum(S(i<j)<0)/sum(~~S(i<j)))
            
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