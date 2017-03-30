function VDFSTCFP = conduitconc

%    mediahead = load('mediahead.txt');
    mediaconc = load('mediaconc.txt');

%    conduithead = load('conduithead.txt');
    conduitconc = importdata('conduitconc.txt');
       
    ncol = 120;
    nlay = 21;
    
    x1 = 1:1:ncol;
    y1 = 1:1:nlay;
    
    for i = 24:25:249
  
        for j = 1:1:nlay
            for k = 1:1:ncol/10          
                for t = 1:1:10
                    
                    col = (k-1)*10 + t;
                    lay = (i-1)*nlay*ncol/10 + (j-1)*ncol/10 + k;
                    mconc(col, j) = mediaconc(lay, t);
                end
            end
        end    
        
        for s = 1:1:120
            if mconc(s, 11) > 10 && mconc(s+1, 11) < 10
                fprintf('%d\n', s);
                break;
            end    
        end
        
        % for conduit
        
%         for c = 1:1:12
%             for s = 1:1:10
%                 node = (c-1)*10+s;
%                 n = (i-1)*12 + c;
%                 cconc(node) = conduitconc(n, s);       
%             end
%         end    
%         
%         
%         for s = 1:1:120
%             if cconc(s) > 10 && cconc(s+1) < 10
%                 fprintf('%d\n', s);
%                 break;
%             end    
%         end
       
       
    end

end
