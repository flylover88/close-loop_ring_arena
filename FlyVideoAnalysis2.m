function [s,sAnalysis,sArena] = FlyVideoAnalysis2(folder,nClip)%,sArena)
close all
load([folder,nClip(1:end-4),'Data'])
v = VideoReader([folder,nClip,'.avi']);
f = v.FrameRate;
len = v.NumberofFrames;

visualization = 0;   
    
angVec = zeros(1,len);
%locations where the ratio of minor axis/major axis is less than 0.9
flag = zeros(1,len);
%locations where the thrust and the slip is beyond normal expectations
flag2 = zeros(1,len);
flag3 = zeros(1,len);
newVideo = cell(1,len);
Kinematics = struct('thrust',[],'slip',[],'yaw',[]);
Distances = struct('DistanceR',[],'xDist',[],'yDist',[]);

s = struct('Kinematics',Kinematics,...
    'Distances',Distances,...
    'MjrAxs',zeros(1,len),...
    'MinAxs',zeros(1,len),...
    'AngVec',zeros(1,len),...
    'Center',struct('x',zeros(1,len),'y',zeros(1,len)),...
    'Head',struct('x',zeros(1,len),'y',zeros(1,len)),...
    'Flags',struct('ratio',flag,'size',flag,'thrust',flag2,'slip',flag3),...
    'Status','good',...
    'LightOn',[]);

sAnalysis = struct('NumPixels',zeros(1,len),...
    'MaxIntensity',zeros(1,len),...
    'MedIntensity',zeros(1,len),...
    'SNRClamp',zeros(1,len),...
    'SNR',zeros(1,len));
    

s.LightOn = LEDon;

if len~=length(LEDon)
    s.status = 'bad';
end
s.time=0:1/f:(len-1)/f;
circleImage = sArena.Mask.Inside;
circleImage2 = ~sArena.bkgMask;

arenaCent=sArena.arenaCent;
rad=sArena.rad;
cF = sArena.cF;

badFrames = false(1,10800);
start = 1; i = start;
redo = 0;
%take first image


while i <= len

    if redo > 0
        i = i-1;
        minPx = minPx-5;
        maxPx = 300;
    else
        minPx = 30;
        maxPx = 300;
    end
        bA = vision.BlobAnalysis('AreaOutputPort',1,...
                        'BoundingBoxOutputPort',0,...
                        'OrientationOutputPort',1,...
                        'MaximumCount',1000,...
                        'MinimumBlobArea',minPx,...
                        'MaximumBlobArea',maxPx,...
                        'ExcludeBorderBlobs',1,...
                        'MinorAxisLengthOutputPort',1,...
                        'MajorAxisLengthOutputPort',1);
    
    %read in frames one by one
    img = v.read(i);
    
    %change image to 2D
    img2 = floor(mean(img,3));imgtemp = img2;
    
%     if i==start
%         centroid = Mask_by_center(img2,circleImage,circleImage2,[],max(max(img2))./2);
%     else
%         centroid = Mask_by_center(img2,circleImage,circleImage2,[s.Center.x(i-1),s.Center.y(i-1)],max(max(img2))./2);
%         if (centroid(1) == 0 || centroid(2) == 0 || badFrames(i-1)) && i<len
%             img = v.read(i+1);
%             %change image to 2D
%             img2 = floor(mean(img,3));
%             centroid = Mask_by_center(img2,circleImage,circleImage2,[],max(max(img2))./2);
%         end
%     end
    
   
   
    
    
     %binarize image
    if redo > 0 
        img3 = imbinarize(img2,6-redo*1);
    else
        img3 = imbinarize(img2,6);
    end
    
    % mask out anything is is outside of the arena
    img3(~circleImage) = 0;
    
    % increase threshold to remove salt and pepper noise
    kk = 1;
    while sum(sum(img3-medfilt2(img3)))>200 && max(max(img2))>8
        img3 = imbinarize(img2,6+kk);
        kk = kk+1;
    end
    img2 = img3;
    
    [area,centroid,mjrAx,minAx,ort] = step(bA,img2);
    
   
    if exist('area')
        [~,idx] = max(area);
    else
        idx=[];
    end
    %if no blob is found or found blob is too small either redo, copy last
    %frame or if it is the first frame set everthing to nan
    if isempty(idx) || max(area)<10
        if redo<4
            redo = redo+1;
        elseif i~=start
            badFrames(i) = true;
            s.MjrAxs(1,i) = s.MjrAxs(1,i-1);    %diameter
            s.MinAxs(1,i) = s.MinAxs(1,i-1);    %diameter
            s.Center.x(i) = s.Center.x(i-1);
            s.Center.y(i) = s.Center.y(i-1);
            s.Head.x(i) = s.Head.x(i-1);
            s.Head.y(i) = s.Head.y(i-1);
            s.Flags.ratio(1,i) = s.Flags.ratio(1,i-1);
            s.Flags.size(1,i) = s.Flags.size(1,i-1);
            redo = 0;
        else
            badFrames(i) = true;
            s.MjrAxs(1,i) = nan;    %diameter
            s.MinAxs(1,i) = nan;    %diameter
            s.Center.x(i) = nan;
            s.Center.y(i) = nan;
            s.Head.x(i) = nan;
            s.Head.y(i) = nan;
            s.Flags.ratio(1,i) = nan;
            s.Flags.size(1,i) = nan;
            redo = 0;
        end
    else
        temp=centroid;
        centroid=temp(idx,:);

        s.Center.x(i) = centroid(1).*cF;
        s.Center.y(i) = centroid(2).*cF;
        
        %Check to make sure change in centroid position is not more than 10pxls
        if i>start && s.Center.x(i-1)>0
            while sqrt(((s.Center.x(i) - s.Center.x(i-1))^2+...
                    (s.Center.y(i) - s.Center.y(i-1))^2))>100 
                img3(:,:) = 0;
                img3(img2) = 1;
                img3 = imfill(img3,'holes');
                
                [area,centroid,mjrAx,minAx,ort] = step(bA,img2);
                
                if isempty(centroid)
                    stats = [];
                    badFrames(i) = true;
                    mjrAx = s.MjrAxs(1,i-1);
                    minAx = s.MinAxs(1,i-1);
                    ort = angVec(1,i-1);
                    index = [];
                    centroid = [s.Center.x(i-1),s.Center.x(i-1)];
                    s.Center.x(i) = centroid(1).*cF;
                    s.Center.y(i) = centroid(2).*cF;
                    s.Flags.ratio(1,i) = s.Flags.ratio(1,i-1);
                else
                    index = find(isnan(area(:,1)));
                    centroid(isnan(centroid(:,1)),:) = [];
                    centroid = mean(centroid,1);
                    s.Center.x(i) = centroid(1).*cF;
                    s.Center.y(i) = centroid(2).*cF;
                end
            end
        end
        newVideo{i} = img3;
        

        s.MjrAxs(1,i) = mjrAx(idx);    %diameter

        s.MinAxs(1,i) = minAx(idx);    %diameter
        
        %orientation (take orientation associated with largest major axis value)
%         ort = cat(1, stats.Orientation);
%         idx = cat(1, stats.MajorAxisLength) == max(cat(1,stats.MajorAxisLength));
        ang(i) = -ort(idx)*180/pi;
        
        %Format where orientation is the angle is between the negative
        %horizontal axis and the major axis.
%         ang = 180-ang;
        if ang<0
            ang(i) = 360+ang(i);
        end
        angVec(1,i) = ang(i);
%         angVec(1,i) = ang+180;
        if angVec(1,i)>360
            angVec(1,i) = angVec(1,i)-360;
        end
        
        %Check to take the absolute smallest change in angle
        if i>start
            angChgVec = [angVec(1,i-1) - ang(i), angVec(1,i-1) - (ang(i)+180),...
                angVec(1,i-1) - (ang(i)-180)];
            [~, imin] = min(abs(angChgVec));
            angChg = angChgVec(imin);
            angVec(1,i) = angVec(1,i-1) - angChg;
        end
        if angVec(1,i)>360
            angVec(1,i) = angVec(1,i)-360;
        elseif angVec(1,i)<0
            angVec(1,i) = angVec(1,i)+360;
        end
        
        %Flag points where the general oval shape of fly is compromised
        if minAx/mjrAx > 0.9
            s.Flags.ratio(1,i) = 10;
        end
        if area<40
            s.Flags.size(1,i) = 10;
        end
        if (area(1)<40) && (redo < 5)
            redo = redo+1;
        else
            redo = 0;
        end
        
        mjrAL =mjrAx(idx);
        minAL = minAx(idx);
        %calculate major and minor axis endpoints
        xMajor1 = centroid(1) + mjrAL * cosd(angVec(1,i));
        xMajor2 = centroid(1) - mjrAL * cosd(angVec(1,i));
        yMajor1 = centroid(2) + mjrAL * sind(angVec(1,i));
        yMajor2 = centroid(2) - mjrAL * sind(angVec(1,i));
        
        xMinor1 = centroid(1) + minAL * cosd(angVec(1,i)+90);
        xMinor2 = centroid(1) - minAL * cosd(angVec(1,i)+90);
        yMinor1 = centroid(2) + minAL * sind(angVec(1,i)+90);
        yMinor2 = centroid(2) - minAL * sind(angVec(1,i)+90);
        
        %calculate head
        s.Head.x(i) = xMajor1;
        s.Head.y(i) = yMajor1;
        
        if visualization == 1 && redo == 0
            %plot image, centroid, ellipse, and the major/minor axis
            figure(10)
            %img3 = insertText(img3,[1 50],i,'AnchorPoint','LeftBottom');
            imagesc(img3)
            hold on
            plot(centroid(1),centroid(2), 'b*')
            if s.LightOn(i)
                plot(arenaCent(2), arenaCent(1),'ro','MarkerSize',70);
            end
            ellipsev2(rad, rad, 0, arenaCent(1), arenaCent(2),'r');
            ellipsev2(mjrAL, minAL, (ang(i))*pi/180, centroid(1), centroid(2),'r');
            plot([xMajor1,xMajor2], [yMajor1,yMajor2],'g')
            plot([xMinor1,xMinor2], [yMinor1,yMinor2],'g')
            plot(s.Head.x(i),s.Head.y(i),'bl.','MarkerSize',30)
            hold off
            title(num2str(i))
        end
    
        if mod(i,1000) == 0
            i
        end
    end
    i = i+1;
end
 s.AngVec = unwrap(angVec)-360;
[s] = postHocAnalysis(badFrames,s,sArena,len);
end
    
    