%{
carfs.VonMises (computed) # my newest table
-> carfs.GratingResponses
---
von_p_value                 : float                         # von mises shuffle p-values
von_r2                      : float                         # R-squared of response
von_pref                    : float                         # von mises preferred direction
von_base                    : float                         # von mises bases
von_amp1                    : float                         # von mises preferred amplitude
von_amp2                    : float                         # von mises anti-preferred amplitude
von_sharp                   : float                         # von mises sharpness
von_osi                     : float                         # orientation selectivity index (pref-ortho)/(pref+ortho)
von_dsi                     : float                         # direction selectivity index (pref-anti)/(pref+anti)
%}

classdef VonMises < dj.Relvar & dj.AutoPopulate
    
    properties
        popRel = carfs.GratingResponseSet
    end
    
    methods
        
        function plot(self)
            % plots the orientation tuning if it exists
            cm = jet(5)*0.6;
            phi = (0:400)/400*(2*pi);
            for key = fetch(carfs.Trace & self)'
                disp(key)
                cellType = fetch1(carfs.Trace & key, 'mask_type');
                [pref,base,amp1,amp2,sharp,spatialFreq, pvalue] = fetchn(carfs.VonMises & key, ...
                    'von_pref','von_base','von_amp1','von_amp2','von_sharp','spatial_freq', 'von_p_value');
                von = ne7.rf.VonMises2([base amp1 amp2 sharp pref]);
                y = von.compute(phi);
                colors = cm(log2(spatialFreq/0.01)+1,:);
                lineWidth = 0.5 + 2.5*(pvalue < 0.05) + 2.0*(pvalue<0.01);
                plotPolar(phi,y,colors,lineWidth)
                f = gcf;
                f.PaperSize = [1 1]*6.0;
                f.PaperPosition = [0 0 f.PaperSize];
                if strcmp(cellType, 'red')
                    hold on
                    plot(0,0,'r*','MarkerSize',40)
                    hold off
                end
                print('-dpng','-r150', ...
                    sprintf('~/Desktop/carfs/%d/%d-%d-ori.png',...
                    key.mouse_id, key.scan_idx, key.masknum))
            end
        end
        
    end
    
    methods(Access=protected)
        
        function makeTuples(self, key)
            
            [responses, traceKeys] = fetchn(carfs.GratingResponses & key, 'spike_responses');
            responses = cell2mat(cellfun(@(m) reshape(m,[1 size(m)]), responses, 'uni', false));
            nTraces = length(traceKeys);
            nShuffles = 10000;
            
            disp 'computing von Mises tuning...'
            [von, r2, p] = ne7.rf.VonMises2.computeSignificance(responses, nShuffles);
            for iTrace = 1:nTraces
                tuple = traceKeys(iTrace);
                tuple.von_r2 = r2(iTrace);
                tuple.von_base = von.w(iTrace,1);
                tuple.von_amp1 = von.w(iTrace,2);
                tuple.von_amp2 = von.w(iTrace,3);
                tuple.von_sharp= von.w(iTrace,4);
                tuple.von_pref = von.w(iTrace,5);
                tuple.von_p_value = p(iTrace);
                
                pref = tuple.von_base + tuple.von_amp1 + tuple.von_amp2*exp(-2*tuple.von_sharp);
                anti = tuple.von_base + tuple.von_amp2 + tuple.von_amp1*exp(-2*tuple.von_sharp);
                orth = tuple.von_base + (tuple.von_amp1 + tuple.von_amp2)*exp(-tuple.von_sharp);
                
                tuple.von_osi = (pref-orth)/(pref+orth);
                tuple.von_dsi = (pref-anti)/(pref+anti);
                self.insert(tuple)
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