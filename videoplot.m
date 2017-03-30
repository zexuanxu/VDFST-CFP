function VDFSTCFP = videoplot

    vidObj = VideoWriter('movie.avi');
    open(vidObj);

%    mediahead = load('mediahead.txt');
    mediaconc = load('mediaconc.txt');

%    conduithead = load('conduithead.txt');
    conduitconc = importdata('conduitconc.txt');
       
    ncol = 120;
    nlay = 21;
    
    x1 = 1:1:ncol;
    y1 = 1:1:nlay;
    
    colormap(flipud(colormap));

    for i = 1:1:249
  
        for j = 1:1:nlay
            for k = 1:1:ncol/10          
                for t = 1:1:10
                    
                    col = (k-1)*10 + t;
                    lay = (i-1)*nlay*ncol/10 + (j-1)*ncol/10 + k;
                    mconc(col, j) = mediaconc(lay, t);
                end
            end
        end    
        
        for c = 1:1:12
            for s = 1:1:10
                node = (c-1)*10+s;
                n = (i-1)*12 + c;
                cconc(node) = conduitconc(n, s);
            end
        end    
                
        mplot = mconc';
        
%         for k = 1:1:nlay
%             for j = 1:1:ncol
%                 mplot(nlay + 1 - k, j) =  mplot(i, j);
%             end
%         end

% incorp conduit

%         for r = 1:1:ncol
%              mplot(11, r) = cconc(r);
%         end
        
%       [m, n] = size(mplot)
%       [m, n] = size(cconc)
                   
        [c, d] = contourf (x1, y1, mplot, 9);
        colorbar;

        set(gca,'ydir','reverse');
        set(gca, 'FontSize',18);
        hold on;

%         plot ([21, 21], [0, 29] , '-r');
%         hold on;
%         plot ([22, 22], [0, 29] , '-r');
%         hold on;
%         plot ([21, 139], [28, 28], '-r');
%         hold on;
%         plot ([21, 139], [29, 29], '-r');
%         hold on;   
%         plot ([138, 138], [0, 29] , '-r');
%         hold on;
%         plot ([139, 139], [0, 29] , '-r');
%         hold on;
        
        xlabel('Distance (*500 ft)');
        ylabel('Depth (*10 ft)');
        title('Salinity distribution along conduit (35 PSU, 0.0 ft at submarine spring)');     
        
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    
    close(vidObj);

end
