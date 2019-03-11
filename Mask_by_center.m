function [centroid] = Mask_by_center(img2,circleImage,circleImage2,oldCenter,threshold)
%check centroid
%binarize image
img4 = imbinarize(img2,threshold);

% mask out anything is is outside of the arena
img4(~circleImage) = 0;

%compute largest connected group of pixels to keep as fly
CC = bwconncomp(img4);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
if max(numPixels)>40
    img4(~circleImage2) = 0;
    CC = bwconncomp(img4);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
end

if isempty(idx)
    centroid =[0,0];
else
    temp1 = CC.PixelIdxList{idx};
    img3 = img2;
    img3(:,:) = 0;
    img3(temp1) = 1;
    
    %calculate centroid, major, minor, and orientation of fly
    stats = regionprops(img3,'Centroid',...
        'MajorAxisLength','MinorAxisLength','Orientation');
    
    %centroid (take mean)
    cent = cat(1, stats.Centroid);
    cent(isnan(cent(:,1)),:) = [];
    centroid = mean(cent,1);
    
    if ~isempty(oldCenter)
        while sqrt((oldCenter(1)-centroid(1)).^2+(oldCenter(2)-centroid(2)).^2)>50 && sum(numPixels)>0
            numPixels(idx) = 0;
            [~,idx] = max(numPixels);
            temp1 = CC.PixelIdxList{idx};
            img3 = img2;
            img3(:,:) = 0;
            img3(temp1) = 1;
            
            stats = regionprops(img3,'Centroid',...
                'MajorAxisLength','MinorAxisLength','Orientation');
            
            cent = cat(1, stats.Centroid);
            cent(isnan(cent(:,1)),:) = [];
            centroid = mean(cent,1);
            
        end
    end
    
end

if sum(numPixels)==0
    img4 = imbinarize(img2,threshold*2/3);
    % mask out anything is is outside of the arena
    img4(~circleImage) = 0;
    
    %compute largest connected group of pixels to keep as fly
    CC = bwconncomp(img4);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    if max(numPixels)>40
        img4(~circleImage2) = 0;
        CC = bwconncomp(img4);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
    end
    
    if isempty(idx)
        centroid =[0,0];
    else
        temp1 = CC.PixelIdxList{idx};
        img3 = img2;
        img3(:,:) = 0;
        img3(temp1) = 1;
        
        %calculate centroid, major, minor, and orientation of fly
        stats = regionprops(img3,'Centroid',...
            'MajorAxisLength','MinorAxisLength','Orientation');
        
        %centroid (take mean)
        cent = cat(1, stats.Centroid);
        cent(isnan(cent(:,1)),:) = [];
        centroid = mean(cent,1);
        
        if ~isempty(oldCenter)
            while sqrt((oldCenter(1)-centroid(1)).^2+(oldCenter(2)-centroid(2)).^2)>50 && sum(numPixels)>0
                numPixels(idx) = 0;
                [~,idx] = max(numPixels);
                temp1 = CC.PixelIdxList{idx};
                img3 = img2;
                img3(:,:) = 0;
                img3(temp1) = 1;
                
                stats = regionprops(img3,'Centroid',...
                    'MajorAxisLength','MinorAxisLength','Orientation');
                
                cent = cat(1, stats.Centroid);
                cent(isnan(cent(:,1)),:) = [];
                centroid = mean(cent,1);
                
            end
        end
        
        
    end
    
    
end
end