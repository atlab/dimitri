classdef plot
    methods(Static)
        
        
        function bowtie
            cm = ne7.vis.doppler(256)*0.9;
            alpha = 'von_p_value<0.01';
            for siteKey = fetch(pro(carfs.TraceSet) & carfs.Subfield & ...
                    carfs.Corr & (carfs.VonMises & alpha))'
                clf
                [R, masknums] = fetch1(carfs.Corr & siteKey, ...
                    'corr_matrix', 'masknums');
                traceKeys = repmat(siteKey, numel(masknums), 1);
                masknums = num2cell(masknums);
                [traceKeys.masknum] = deal(masknums{:});
                
                ixTuned = arrayfun(@(k) count(carfs.VonMises & k & alpha)>0, traceKeys);
                ixField = arrayfun(@(k) count(carfs.Subfield & k)>0, traceKeys);
                scale = 20;
                
                for refCell = find(ixTuned & ixField)'
                    for otherCell = setdiff(find(ixTuned & ixField),refCell)'
                        % plot tuning of refCell
                        von0 = fetch(carfs.VonMises & traceKeys(refCell), '*');
                        [~,j] = max([von0.von_amp1]);
                        von0 = von0(j);
                        von1 = fetch(carfs.VonMises & traceKeys(otherCell), '*');
                        [~,j] = max([von1.von_amp1]);
                        von1 = von1(j);
                        
                        color = cm(round(max(1, min(255, R(refCell,otherCell)/0.05*128+127.5))),:);
                        % ppolar([0 0], setfield(von0, 'von_pref', 0), color, scale)
                        hold on
                        
                        % translate and rotate coordinates of otherCell
                        [x0, y0] = fetch1(carfs.Subfield & traceKeys(refCell),'sub_x','sub_y');
                        [x1, y1] = fetch1(carfs.Subfield & traceKeys(otherCell),'sub_x','sub_y');
                        
                        x1 = x1 - x0;
                        y1 = y1 - y0;
                        
                        delta = von1.von_pref - von0.von_pref;
                        [x1, y1] = deal(cos(delta)*x1 + sin(delta)*y1, -sin(delta)*x1+cos(delta)*y1);
                        
                        % plot tuning of otherCell
                        von1.von_pref = delta;
                        ppolar([x1 y1], von1, color, scale)
                        drawnow
                    end
                end
                hold off
                axis image
                xlabel '\Deltax'' (degrees)'
                ylabel '\Deltay'' (degrees)'
                grid on                
                set(gcf, 'PaperSize', [20 20], 'PaperPosition', [0 0 20 20])
                print('-dpdf', sprintf('~/Desktop/bowtie-%04d-%2d.pdf', siteKey.mouse_id, siteKey.scan_idx))
            end
            
            function ppolar(xy, von, color, scale)
                rotation = pi/2;
                
                phi = linspace(0, 2*pi, 100);
                von = ne7.rf.VonMises2([von.von_base von.von_amp1 von.von_amp2 von.von_sharp von.von_pref]);
                
                y = von.compute(phi-rotation);
                
                %                 mx = max(1e-6,max(abs(y(:))));
                %                 plot(sin(phi)*mx,cos(phi)*mx,'k:','LineWidth',0.25)
                %                 hold on
                %                 plot(sin(phi)*mx/2,cos(phi)*mx/2,'k:','LineWidth',0.25)
                %                 plot([0 0],[-mx mx],'Color',[1 1 1]*0.7,'LineWidth',0.25)
                %                 plot([-mx mx],[0 0],'Color',[1 1 1]*0.7,'LineWidth',0.25)
                for i=1:size(y,1)
                    plot(xy(1)+scale*sin(phi).*y(i,:),xy(2)+scale*cos(phi).*y(i,:),'LineWidth',.5,'Color',color)
                end
                %                 axis equal
                %                 axis([-mx mx -mx mx])
                %                 axis off
                %                 hold off
            end
        end
        
        
        
        function sidebyside(r)
            
            cm = jet(5)*0.6;
            phi = linspace(0,2*pi,400);
            
            for siteKey = fetch(pro(carfs.TraceSet) & carfs.NoiseMapSet & carfs.VonMises & r)'
                clf
                cells = carfs.NoiseMap & (carfs.MatchedControl & 'mask_type="red"') & siteKey;
                nCells = cells.count;
                
                sub = 0;
                for key = cells.fetch'
                    sub = sub+1;
                    subplot(nCells,4,sub)
                    map = fetch1(carfs.NoiseMap & key, 'noise_map');
                    imagesc(-max(abs(map),[],3))
                    axis image
                    axis off
                    colormap gray
                    
                    title(sprintf('%d',key.masknum))
                    
                    sub = sub+1;
                    subplot(nCells,4,sub)
                    [pref,base,amp1,amp2,sharp,spatialFreq, pvalue] = fetchn(carfs.VonMises & key, ...
                        'von_pref','von_base','von_amp1','von_amp2','von_sharp','spatial_freq', 'von_p_value');
                    von = ne7.rf.VonMises2([base amp1 amp2 sharp pref]);
                    y = von.compute(phi);
                    colors = cm(log2(spatialFreq/0.01)+1,:);
                    lineWidth = 0.5 + 2.5*(pvalue < 0.05) + 2.0*(pvalue<0.01);
                    plotPolar(phi,y,colors,lineWidth)
                    
                    key = fetch(carfs.Trace & carfs.MatchedControl*pro(carfs.MatchedControl & key, 'mask_type->c') & 'mask_type="neuron"',  'masknum');
                    
                    sub = sub+1;
                    subplot(nCells,4,sub)
                    map = fetch1(carfs.NoiseMap & key, 'noise_map');
                    imagesc(-max(abs(map),[],3))
                    axis image
                    axis off
                    colormap gray
                    title(sprintf('%d',key.masknum))
                    
                    sub = sub+1;
                    subplot(nCells,4,sub)
                    [pref,base,amp1,amp2,sharp,spatialFreq, pvalue] = fetchn(carfs.VonMises & key, ...
                        'von_pref','von_base','von_amp1','von_amp2','von_sharp','spatial_freq', 'von_p_value');
                    von = ne7.rf.VonMises2([base amp1 amp2 sharp pref]);
                    y = von.compute(phi);
                    colors = cm(log2(spatialFreq/0.01)+1,:);
                    lineWidth = 0.5 + 2.5*(pvalue < 0.05) + 2.0*(pvalue<0.01);
                    plotPolar(phi,y,colors,lineWidth)
                    drawnow
                    
                    
                end
                fig = gcf;
                fig.PaperSize = [8 10];
                fig.PaperPosition = [0 0 fig.PaperSize];
                
                print('-dpdf','-r300',sprintf('~/Desktop/r%4d-%2d', key.mouse_id, key.scan_idx))
                
            end
            
        end
        
        
        function dspace(r)
            if nargin<1
                r = {};
            end
            
            f0 = fetchn(carfs.BestSpatialFrequency*carfs.Trace & 'mask_type="neuron"' & r, 'best_spatial_freq');
            f1 = fetchn(carfs.BestSpatialFrequency*carfs.Trace & 'mask_type="red"' & r,'best_spatial_freq') ;
            
            p = ranksum(f0,f1)
            bins = 0:5;
            a = hist(log2(f0/0.01)+0.5, bins);
            b = hist(log2(f1/0.01)+0.5, bins);
            a = a/sum(a);
            b = b/sum(b);
            
            bar(bins, [a;b]')
            set(gca, 'XTick', bins, 'XTickLabel', arrayfun(@(f) sprintf('%1.2g', f), 2.^bins*0.01, 'uni', false))
            legend controls sisters
            
        end
        
        
        function dor(r)
            
            cells = carfs.MatchedControl * carfs.VonMises & r & 'von_p_value<0.05';
            s1 = fetch(cells & 'mask_type="red"', 'von_pref', 'von_amp1');
            [d1, m1] = dj.struct.tabulate(s1,'von_pref','masknum');
            [a1, m1] = dj.struct.tabulate(s1,'von_amp1','masknum');
            [~,j] = nanmax(a1,[],2);
            d1 = arrayfun(@(i) d1(i,j(i)), 1:length(j));
            
            s0 = fetch(pro(carfs.MatchedControl & 'mask_type="red"','masknum')*pro(cells & 'mask_type="neuron"','mask_type->mt','masknum->m','von_pref','von_amp1'), 'masknum', 'von_pref', 'von_amp1');
            [d0, m0] = dj.struct.tabulate(s0,'von_pref','masknum');
            [a0, m0] = dj.struct.tabulate(s0,'von_amp1','masknum');
            [~,j] = nanmax(a0,[],2);
            d0 = arrayfun(@(i) d0(i,j(i)), 1:length(j));
            
            d0 = d0*180/pi;
            d1 = d1*180/pi;
            m = 180;
            nbins = 6;
            step = (m/2)/nbins;
            d0 = min(mod(d0,m),mod(-d0,m));
            d1 = min(mod(d1,m),mod(-d1,m));
            
            p = signrank(d0,d1)
            
            hist([d1; d0]', step/2:step:m/2)
            legend sister-sister control-control
            
        end
        
        
        
        
        function polar(r)
            if nargin<1
                r = {};
            end
            
            subplot(4,3,1:9)
            carfs.plot.dor(r)
            sub = 9;
            for siteKey = fetch(carfs.TraceSet & carfs.MatchedControl & r)'
                sub = sub + 1;
                subplot(4,3,sub)
                s0 = fetch(carfs.VonMises & siteKey & 'von_p_value<0.05' & (carfs.MatchedControl & 'mask_type="red"'), ...
                    'von_pref', 'von_amp1');
                [d0, m] = dj.struct.tabulate(s0, 'von_pref', 'masknum');
                [a0, m] = dj.struct.tabulate(s0, 'von_amp1', 'masknum');
                [~,j] = nanmax(a0,[],2);
                d0 = arrayfun(@(i) d0(i,j(i)), 1:length(j));
                compass(exp(1i*d0),'r')
                
                hold on
                
                s0 = fetch(carfs.VonMises & siteKey & 'von_p_value<0.05' & (carfs.MatchedControl & 'mask_type="neuron"'), ...
                    'von_pref', 'von_amp1');
                [d0, m] = dj.struct.tabulate(s0, 'von_pref', 'masknum');
                [a0, m] = dj.struct.tabulate(s0, 'von_amp1', 'masknum');
                [~,j] = nanmax(a0,[],2);
                d0 = arrayfun(@(i) d0(i,j(i)), 1:length(j));
                compass(exp(1i*d0),'k')
                
                hold off
            end
        end
        
    end
    
end


function plotPolar(phi,y,colors,lineWidth)
mx = max(1e-6,max(abs(y(:))));
plot(sin(phi)*mx,cos(phi)*mx,'k:','LineWidth',0.25)
hold on
plot(sin(phi)*mx/2,cos(phi)*mx/2,'k:','LineWidth',0.25)
plot([0 0],[-mx mx],'Color',[1 1 1]*0.7,'LineWidth',0.25)
plot([-mx mx],[0 0],'Color',[1 1 1]*0.7,'LineWidth',0.25)
for i=1:size(y,1)
    plot(sin(phi).*y(i,:),cos(phi).*y(i,:),'LineWidth',lineWidth(i),'Color',colors(i,:))
end
axis equal
axis([-mx mx -mx mx])
axis off
hold off
end