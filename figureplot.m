function VDFSTCFP = figureplot

%    mediahead = load('mediahead.txt');
    mediaconc = load('mediaconc.txt');

%    conduithead = load('conduithead.txt');
    conduitconc = importdata('conduitconc.txt');
       
    ncol = 120;
    nlay = 21;
    
    x1 = 1:1:ncol;
    y1 = 1:1:nlay;
    
    colormap(flipud(colormap));

    for i = 249:1:249
  
        for j = 1:1:nlay
            for k = 1:1:ncol/10          
                for t = 1:1:10
                    
                    col = (k-1)*10 + t;
                    lay = (i-1)*nlay*ncol/10 + (j-1)*ncol/10 + k;
                    mconc(col, j) = mediaconc(lay, t);
                end
            end
        end    
        
%         for c = 1:1:18
%             for s = 1:1:10
%                 if c == 18 && s >= 5
%                 else
%                     node = (c-1)*10+s;
%                     n = (i-1)*18 + c;
%                     cconc(node) = conduitconc(n, s);
%                 end
%             end
%         end    
                
        mplot = mconc';
        
%         for i = 1:1:37
%             for j = 1:1:140
%                 mplot(38 - i, j) =  mplot(i, j);
%             end
%         end
        
        
%       [m, n] = size(mplot)
%       [m, n] = size(cconc)
                   
        [c, d] = contourf (x1, y1, mplot, 9);
        colorbar;

        set(gca,'ydir','reverse');
        set(gca, 'FontSize',18);
        hold on;
 
 plot ([0, 120], [10.5, 10.5] , '-r');
 hold on;

 plot ([0, 120], [11.5, 11.5] , '-r');
 hold on;
        
        
    xlabel('x (south - north)');
    ylabel('y (east - west)');
    title('Salinity Profile ');   
       
    end

end
