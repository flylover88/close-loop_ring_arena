clear;clc

currFolder = pwd;
% [status, message, messageid] = rmdir([currFolder '\clips']);
% if ~status
%     input('Clip directory still has files in it. Press enter to still delete?')
%     rmdir([currFolder '\clips'],'s');
% end
addpath([currFolder '\DataOld']);
addpath([currFolder '\DataNew']);
addpath([currFolder '\DataGen']);
addpath([currFolder '\Videos\ToDo']);
addpath([currFolder '\utils']);
addpath([currFolder '\Analysis']);

%% IF YOU HAVE BKGS CLIPS IN THE CLIP FOLDER ALONG WITH THEIR .MAT FILE YOU
%CAN SET justAnalysis to 1 and only do analysis
justAnalysis=0;


    
FPS = 30;
mF=10;
mD=10;
if ~justAnalysis
%find video files with mathcing data files
vid_Ring = dir([currFolder,'\Videos\ToDo\*.avi']);
%FicTrac Videos
nVid = length(vid_Ring);
Ctri = 1;
%checks to make sure there is a matching .mat file for the video
for z = 1:nVid
    tmp = vid_Ring(z).name;
    tmp2 = dir([currFolder,'\Videos\ToDo\',tmp(1:end-4),'.mat']);
    if size(tmp2,1) == 1
        trial{Ctri} = tmp(1:end-4);
        Ctri = Ctri +1;
    end
end
numTrial = Ctri - 1;
clear vid_Ring z Ctri nVid tmp tmp2



%%This first section will seperate the videos into clips, build the
%%background and then svae bkgs versions of the clips along with an index
%%of frames with the LED on(LEDon=1 for frames with LED one and 0 for
%%frames with LED off)
for z = 1:numTrial
    t=trial{z};
    tarFolder = [currFolder,'\Videos\ToDo\'];
    clear Date%added 11202018
    load([tarFolder,t,'.mat']);
    if ~exist('Date')%added 11202018
        Date = input('What was the date of this run?');
    end
    dates{z} = Date;%added 11202018
    
        
%     if exist('daMask')
%         Mask = daMask;
%     end
    
    
    
    %seperate runs into individual .avi clips
    try
        v = VideoReader([tarFolder,t,'.avi']);
    catch
        continue
    end
    cframe = 1:1:v.NumberofFrames;
    isFrameNdx = find(frame_timeStamp>0)';
    f = find(diff(isFrameNdx)>1);%find all point frame count jumps more than one-2clips meet
    
    lastFrames = unique([cframe(f) max(cframe)]);
    firstFrames = unique([min(cframe) cframe(f+1)]);
%     if max(lastFrames)~=max(cframe)
%         lastFrames = [lastFrames; max(cframe)];
%     end
%     if min(firstFrames)~=min(cframe)
%         firstFrames = [min(cframe); firstFrames];
%     end
    
    numFrames = lastFrames - firstFrames +1;
    
    
    
    mkdir([currFolder, '\clips\',t])
    addpath([currFolder '\clips\',t]);
    
    
    %break video in to clips
    c=1;
    ff=1;
    while c<length(numFrames)+1
        vo=VideoWriter([currFolder '\clips\',t,'\clip',int2str(ff),'.avi'],'Motion JPEG AVI');
        clipPath=[vo.Path];
        if numFrames(c)>FPS/2
            open(vo)
            for cc=firstFrames(c):lastFrames(c)
                if cc>max(cframe)
                    close(vo)
                    break
                end
                frame = read(v,cc);
                writeVideo(vo,frame)
            end
            close(vo)
            ff=ff+1;
        end
        c=c+1;
    end
    clear *rame* f c cc ff
    
    %Find arena parameters(this is why the .mat file is needed)
    sArena.arenaCent = [mean(Mask.Vertices(:,2)) mean(Mask.Vertices(:,1))];
    sArena.rad = (max(Mask.Vertices(:,1))-min(Mask.Vertices(:,1)) +...
        max(Mask.Vertices(:,2))-min(Mask.Vertices(:,2)))/4;
    sArena.cF = 40/sArena.rad;
   %At this point clips can be user deleted so that they wont be continued
   %down the pipe to be bkgs and analysised
    %% make background subtracted clips WARNING WILL DELETE SAVED BKGS VIDEOS IN CLIPS FOLDER
    clips = dir([clipPath,'\clip*.avi']); 
    numClips = size(clips,1);
    c=1;
    clear nclips
    for zz = 1:numClips
        v = VideoReader([clipPath,'\',clips(zz).name]);
        if v.Duration>.5
            if clips(zz).name(end-4) ~= 's'%make sure the clip isnt bkgs 
                nclips{c} = clips(zz).name;
                c=c+1;
            else
            clear v
            delete([clipPath,'\',clips(zz).name])
            end
        
        else
            clear v
            delete([clipPath,'\',clips(zz).name])
        end
    end
    
    numClips = c-1;
    clear c
    
    v = VideoReader([clipPath,'\',nclips{numClips}]);
    
    [X,Y] = meshgrid(1:v.Width,1:v.Height);
    bkgMask = (X - sArena.arenaCent(2)).^2 ...
        + (Y - sArena.arenaCent(1)).^2 < (sArena.rad/10)^2;
    
    bkgMask2 = (X - sArena.arenaCent(2)).^2 ...
        + (Y - sArena.arenaCent(1)).^2 < (sArena.rad/1.8)^2;
    % this value (1.8) might need to be changed if the recording window
    % changes
    
    
    
    %%
    %makes background from all clips for this video, didnt want to do over
    %all videos because the camera might have moved
    for c = 1:numClips
        clear v
        v = VideoReader([clipPath,'\',nclips{c}]);
        nframes = v.NumberofFrames;
        center = [v.Height/2 v.Width/2];
%         rad = mean(center)/4;
%         LEDon=zeros(1,nframes);
       %Read clips frame by frame, decides if LED is on/off
       %%
       b=1;
       bb=1;
        for cc = 1:nframes
            frame = read(v,cc);
            frame_sing = uint8(floor(mean(frame,3)));
            frame_blur = imgaussfilt(frame_sing);
            frame_bi = imbinarize(frame_blur,.15);
            frame_bi2 = imbinarize(frame_blur,.2);
            
            CC = bwconncomp(and(frame_bi,bkgMask));
            S = regionprops(CC,'Centroid','Area','BoundingBox');
            areas = cat(1,S.Area);
            flyMask = and(frame_bi2,bkgMask2);
            frame_blur = uint8(frame_blur.*uint8(~flyMask));%mask out fly and bring to normal value
            frame_blur(frame_blur==0)=30;
            if length(areas)>=1
                %determin if LED is on
                if max(areas)>500
                    LEDon = 1;
                else
                    LEDon = 0;
                end
            else
                LEDon=0;
            end
            %Make bkg for when LED is on
            if LEDon
%                 if exist('LEDon_bkg')
%                     LEDon_bkg = (LEDon_bkg + frame_blur)./2;
%                 else
                    LEDon_bkg1(:,:,b) = frame_blur;
%                     imshow(LEDon_bkg(:,:,b))
                    b=b+1;
%                 end
            else %Make bkg for when LED is off
%                 if exist('LEDoff_bkg')
%                     LEDoff_bkg = (LEDoff_bkg + frame_blur)./2;
%                 else
                    LEDoff_bkg1(:,:,bb) = frame_blur;
%                     imshow(LEDoff_bkg(:,:,bb))
                    bb=bb+1;
                    
%                 end
            end
        end
        if exist('LEDon_bkg1')
        LEDon_bkg2(:,:,c) = uint8(mean(LEDon_bkg1,3));
        end
        if exist('LEDoff_bkg1')
        LEDoff_bkg2(:,:,c) = uint8(mean(LEDoff_bkg1,3));
        end
    end

    if exist('LEDoff_bkg2')
    LEDoff_bkg = uint8(mean(LEDoff_bkg2,3));
    LEDoff_bkg(LEDoff_bkg<30) = 30;
    end
    if exist('LEDon_bkg2')
    LEDon_bkg = uint8(mean(LEDon_bkg2,3));
    LEDon_bkg(LEDon_bkg<30) = 30;
    end
    clear frame* S areas c cc LED*_bkg1 LED*_bkg2 b bb flyMask nframes 
    %% BKG has been made from all clips from an individual video by this point
    % Begin BKGS
    
    
    for c = 1:numClips
        clear v vo LEDon
        vo=VideoWriter([currFolder '\clips\',t,'\',nclips{c}(1:end-4),'_bkgs.avi'],'Motion JPEG AVI');
        v = VideoReader([clipPath,'\',nclips{c}]);
        nframes = v.NumberofFrames;
        center = [v.Height/2 v.Width/2];
        rad = mean(center)/4;
        
        open(vo)
        
       %again goes through each frame and decides if LED is on/off but this
       %time BKGS
        for cc = 1:nframes
            frame = read(v,cc);
            frame_sing = uint8(floor(mean(frame,3)));
            frame_blur = imgaussfilt(frame_sing);
            frame_bi = imbinarize(frame_blur,.14);
            CC = bwconncomp(and(frame_bi,bkgMask));
            S = regionprops(CC,'Centroid','Area');
            areas = cat(1,S.Area);
            if length(areas)>=1
                
                if max(areas)>400
                    LEDon(cc) = 1;
                else
                    LEDon(cc) = 0;
                end
            elseif cc>1
                LEDon(cc)=LEDon(cc-1);
            else
                LEDon(cc)=0;
            end
            
            if LEDon(cc)
                frame_bkgs = (frame_blur - LEDon_bkg).*uint8(Mask.Inside);
            else
                frame_bkgs = (frame_blur - LEDoff_bkg).*uint8(Mask.Inside);
            end
            %writes BKS video here
            writeVideo(vo,frame_bkgs)
            
        end
        close(vo)
%         LEDon2 = LEDon(c,:);
        sArena.f = FPS;
        sArena.Mask = Mask;
        sArena.bkgMask = bkgMask;
        save([currFolder '\clips\',t,'\',nclips{c}(1:end-4),'_Data'],'LEDon','sArena')
        
    end
    
    clear clips numClips nframes frame* CC c cc areas X Y v vo nclips LED* S zz rad 
    
end
clear z t trial
end
%% RUN ALL BELOW TO TURN BKGS CLIPS TO DATA AND FIGURES


%If any BKGS clips are not optimal they can be deleted here. The matching
%.mat can also be deleted but does not have to be. Any BKGS clips that do
%not have a corisponding Data.mat file will be deleted by the next section
%as well. At the end in the /clips folder you will have the raw clips, the
%BKGS clips and a Data.mat file that has assay info and LED index. The
%analysed data will go to /DataNew and low level figures will be in
%/Figures (saved in format Video*_clip*) Afterwards the user will need to
%move completed videos into an archive along with the contents of /clips,
%or just deleted them.




    %% Atthis point all clips that are longer than .5s have been bkgs and an index of frames with the LED on have been created
    % The pairs have been saved in /clips
    %find clips and LED data and make sure they are in pairs
    sVid = dir([currFolder, '\clips']);
    numVids = length(sVid);
    cc=1;
    for c = 1:numVids
        tmp=sVid(c).name;
        if length(tmp)>4
            if tmp(1:5)=='video'
                vids{cc}=sVid(c).name;
                cc=cc+1;
            end
        end
    end
    numVids = cc-1;
    clear sVid cc c tmp
    
 if ~exist('dates')   
    for z=1:numVids 
        dates{z} = input(['Enter the date for ' vids{z} ': ']);
     end
 end

    
    
    
 for z=1:numVids
 vidFolder=[currFolder,'\clips\',vids{z}];
 bkgsClips = dir([vidFolder,'\*bkgs.avi']); 
    
    cc=1;
    numClips = length(bkgsClips);
    for c = 1:numClips
        tmp = bkgsClips(c).name(1:end-4);
        LEDdata = dir([vidFolder, '\',tmp(1:end-4),'Data.mat']);
        if length(LEDdata)==1
            nClip{cc}=tmp;
            nLED{cc}=LEDdata.name(1:end-4);
            cc=cc+1;
        else
            delete([vidFolder, '\',tmp,'.avi'])
        end
    end
    numClips = cc-1;
    
    clear tmp c cc bkgsClips
%     sArena.f = FPS;
%     sArena.Mask = Mask;
%     sArena.bkgMask = bkgMask;

%Begin analysis and 
b=1;
clear fKin Data DataNew M
%%
    for c = 1:numClips
        
        [s,sAnalysis,sArena] = FlyVideoAnalysis2([vidFolder,'\'],nClip{c});%,sArena);
%         save([currFolder,'\DataNew\',nClip{c}], 's', 'sArena', 'sAnalysis')
%         movefile(

%Clip runs to in light exposure and 10 frames before (making multiple runs
%per clip if needed and ignoring clips without runs)
        clear L*
        Li = s.LightOn(1:length(s.Kinematics.thrust));
        r=s.Distances.DistanceR(1:length(s.Kinematics.thrust));
        ti = s.time(1:length(s.Kinematics.thrust));
        Li(r>12) = 0;
        Li(smooth(Li,5)>.5)=1;
        
        Lindx = find(diff(Li)==1)+1;%first frame in
        Londx = find(diff(Li)==-1);%last frame in
        
        if length(Lindx)>length(Londx)
            Londx = [Londx length(Li)];
        end
        
        tB = 10;%number of frames to use before fly goes in light(15 might be the max?)
        % NOTE first frame in it tB+1!!!
        tA = 300; %max number of frames after fly goes in light
        for cc = 1:length(Lindx)
            if Londx(cc)-Lindx(cc)>mF &&min(r(Lindx(cc):Londx(cc)))<mD
                pos = [s.Center.x(Lindx(cc)-tB:min(Londx(cc),Lindx(cc)+tA))', s.Center.y(Lindx(cc)-tB:min(Londx(cc),Lindx(cc)+tA))'];
                th = s.Kinematics.thrust(Lindx(cc)-tB:min(Londx(cc),Lindx(cc)+tA))';
                sl = s.Kinematics.slip(Lindx(cc)-tB:min(Londx(cc),Lindx(cc)+tA))';
                ya = s.Kinematics.yaw(Lindx(cc)-tB:min(Londx(cc),Lindx(cc)+tA))';
                speed = sqrt(diff(s.Center.x(Lindx(cc)-tB:min(Londx(cc),Lindx(cc)+tA)+1))'.^2 ...
                    + diff(s.Center.y(Lindx(cc)-tB:min(Londx(cc),Lindx(cc)+tA)+1))'.^2);
                [curv, theta] = CalcCurvature(s);
                curv = unwrap(rad2deg(curv((Lindx(cc)-tB:min(Londx(cc),Lindx(cc)+tA)))'));
                ti = [0:1:length(th)-1]'./sArena.f;
                M(:,:,b) = nan(tB+tA+1,8);
                M(1:length(pos),1:8,b) = [pos,th,sl,ya,curv,speed,ti];

                fKin{b,1} = M(:,:,b);
                fKin{b,2} = nClip{c};
                b=b+1;
                clear th sl ya speed curv theta pos ti
            end
        end
    end  
    clear c  b Li* Data Londx
    %%
    if ~exist('fKin')
        fKin = {};
    else
        
        
    %transform xy pos so that the fly is heading in the positive x
    %direction over the first 10 frames
    Data =fKin(:,1);
    DataNew = cell(size(Data)); 
    for i = 1:length(Data)
       currTrack = Data{i}(:,1:2)';
%        currTrack(1,:) = currTrack(1,:)-currTrack(1,11);                        % make everything relative to the 11th point of trajectory
%        currTrack(2,:) = (currTrack(2,:)-currTrack(2,11));
       xy  = currTrack(1,:)+1i*currTrack(2,:);                                 % change to polar coordinates for easier calculation
       vel_x=diff(currTrack(1,:)); vel_y=diff(currTrack(2,:));
       vel=vel_x+vel_y*1i;
       theta=angle(sum(vel(1:tB)));                                            % Use the first 10 steps for the heading direction
       xy2=xy.*exp(-theta*1i);%+(pi/2)*1i);                                       % shift by 90 degrees to point up instead of to the right
       newY = imag(xy2);
       newX=real(xy2);
       
       newY = newY-newY(tB);
       newX = newX-newX(tB);
       
       DataNew{i} = [newX',newY',Data{i}(:,3:8)];
    end
    clear xy* vel* newX newY Data i
    

    avgM=mean(M,3,'omitnan');%wants the nan padding
    
    %% BEGIN PLOTTING-----------------------------------------------------------
    timeEnd=3;
    tS=5; %vertical spacing runs in plots
    %plot tracks
    figure(21)
    subplot(4,1,1)
    hold on
    for c=1:length(DataNew)
        x=DataNew{c}(:,1);
        y=DataNew{c}(:,2) + c*tS;
        
        plot(x(1:tB),y(1:tB),'bl')%,'DisplayName','Light Off')
        plot(x(tB:end),y(tB:end),'r')%,'DisplayName','Light On')

        plot(x(5:5:end),y(5:5:end),'black.')
        plot(x(tB),y(tB),'o','MarkerSize',10)
        
         text(-30,c*tS,['Fly ' fKin{c,2}(5:end-5)])
        xlabel('(mm)')
        xlim([-3 7])
        ylim([-10 c*tS+15])
        title([vids{z} '-' dates{z}])%% added title on 11202018
        clear x y
    end
    hold off
%     print([vidFolder, '\Track'],'-dpdf','-fillpage')
    
    
    
    subplot(4,2,3)
    %plot thrust
    hold on
    ymin=0;ymax=0;
    for c=1:length(DataNew)
        clear th ti
        th=DataNew{c}(:,3);%  + c*20 -DataNew{c}(10,3);
        ti=DataNew{c}(:,8);
        
        plot(ti(1:end),th(1:end))
        
        ymax = max([ymax;th]);
        ymin = min([ymin;th]);

        ylabel('Thrust')
        
        
    end
    ylim([ymin ymax])
    xlim([0 timeEnd])
    plot([ti(tB+1) ti(tB+1)],[ymin ymax],'black:')
    plot(avgM(:,8),avgM(:,3),'black','LineWidth',1.5)
    hold off
%     print([vidFolder, '\Thrust'],'-dpdf')
    
    
    subplot(4,2,5)
    %plot slip
    hold on
    ymin=0;ymax=0;
    for c=1:length(DataNew)
        clear sl ti
        sl=DataNew{c}(:,4);%  + c*20 -DataNew{c}(10,4);
        ti=DataNew{c}(:,8);
        
        plot(ti(1:end),sl(1:end))
        
        ymax = max([ymax;sl]);
        ymin = min([ymin;sl]);
        
        ylabel('Slip')
        
        
    end
    ylim([ymin ymax])
    xlim([0 timeEnd])
    plot([ti(tB+1) ti(tB+1)],[ymin ymax],'black:')
    plot(avgM(:,8),avgM(:,4),'black','LineWidth',1.5)
    hold off
%     print([vidFolder, '\Slip'],'-dpdf')
    
    
    subplot(4,2,7)
    %plot yaw
    hold on
    ymin=0;ymax=0;
    for c=1:length(DataNew)
        clear ya ti
        ya=DataNew{c}(:,5);%  + c*3 -DataNew{c}(10,5);
        ti=DataNew{c}(:,8);
        
        plot(ti(1:end),ya(1:end))
        
        ymax = max([ymax;ya]);
        ymin = min([ymin;ya]);
        
        ylabel('Yaw')
        
    end
    ylim([ymin ymax])
    xlim([0 timeEnd])
    plot([ti(tB+1) ti(tB+1)],[ymin ymax],'black:')
    plot(avgM(:,8),avgM(:,5),'black','LineWidth',1.5)
    hold off
%     print([vidFolder, '\Yaw'],'-dpdf')
    
    
    subplot(4,2,4)
    %plot curve
    hold on
    ymin=0;ymax=0;
    for c=1:length(DataNew)
        clear cr ti
        cr=unwrap(DataNew{c}(:,6));%  +c*20 -unwrap(DataNew{c}(10,6));
        ti=DataNew{c}(:,8);
        
        plot(ti(1:end),cr(1:end))
        
        ymax = max([ymax;cr]);
        ymin = min([ymin;cr]);
        
        ylabel('Curve')
        
    end
    ylim([ymin ymax])
    xlim([0 timeEnd])
    plot([ti(tB+1) ti(tB+1)],[ymin ymax],'black:')
    plot(avgM(:,8),avgM(:,6),'black','LineWidth',1.5)
    hold off
%     print([vidFolder, '\Curve'],'-dpdf')
    
    
    
    subplot(4,2,6)
    %plot speed
    hold on
    ymin=mean(avgM(:,7));ymax=0;
    for c=1:length(DataNew)
        clear sp ti
        sp=DataNew{c}(:,7);%  + c*3 -DataNew{c}(10,7);
        ti=DataNew{c}(:,8);
         
        plot(ti(1:end),sp(1:end))
        
        ymax = max([ymax;sp]);
        ymin = min([ymin;sp]);
        
        ylabel('Speed')
        
    end
    ylim([ymin ymax])
    xlim([0 timeEnd])
    plot([ti(tB+1) ti(tB+1)],[ymin ymax],'black:')
    plot(avgM(:,8),avgM(:,7),'black','LineWidth',1.5)
    hold off
    
    
%     subplot(4,2,8)
%     hold on
%     P=10;
%     for c=1:length(DataNew)
%     plot(c*P,'DisplayName',fKin{c,2}(1:end-5))
%     end
%     hold off
% %     ylim([P c*P])
% %     xlim([.7 2.5])
%     lgd=legend('Location','Best');
% %     lgd.NumColumns=4;
%     ax= gca;
%     axis off
    
    
    print([vidFolder, '\Plot'],'-dpdf','-fillpage')
    end
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
%         Li0=Li; Li1=Li0;
%         Li1(diff(Li1)>0)=1;
%         Li0(diff(Li0)<0)=0;
%         
%         LiXY = s.LightOn(1:length(s.Distances.yDist));
%         rXY=s.Distances.DistanceR(1:length(s.Kinematics.thrust));
%         LiXY(rXY>12) = 0; 
%         
%         LiXY0=LiXY1;
%         LiXY1 = LiXY0;
%         LiXY1(diff(LiXY1)>0)=1;
%         LiXY0(diff(LiXY0)<0)=0;
%         
%         th0=s.Kinematics.thrust;  th1=s.Kinematics.thrust;
%         th0(Li0==1)=nan;
%         th1(Li1==0)=nan;
%         
%         sl0=s.Kinematics.slip;  sl1=s.Kinematics.slip;
%         sl0(Li0==1)=nan;
%         sl1(Li1==0)=nan;
%         
%         ya0=s.Kinematics.yaw;  ya1=s.Kinematics.yaw;
%         ya0(Li0==1)=nan;
%         ya1(Li1==0)=nan;
%         
%         y0=s.Distances.yDist; y1=s.Distances.yDist;
%         y0(LiXY0==1)=nan;
%         y1(LiXY1==0)=nan;
%         
%         x0=s.Distances.xDist; x1=s.Distances.xDist;
%         x0(LiXY0==1)=nan;
%         x1(LiXY1==0)=nan;
%         
%         
%         figure(1)
%         subplot(311)
%         hold on
%         plot(t,th0,'bl')
%         plot(t,th1,'r')
%         hold off
%         ylabel('thrust')
%         legend('Light out','Light in','Location','Best')
%         
%         subplot(312)
%         hold on
%         plot(t,sl0,'bl')
%         plot(t,sl1,'r')
%         hold off
%         ylabel('slip')
%         legend('Light out','Light in','Location','Best')
%         
%         subplot(313)
%         hold on
%         plot(t,ya0,'bl')
%         plot(t,ya1,'r')
%         hold off
%         ylabel('yaw')
%         legend('Light out','Light in','Location','Best')
%         xlabel('sec')
%         
% %         print([currFolder,'\Figures\',nClip{c},'TSY'],'-dpdf')
%         
%         figure(2)
%         hold on
% %         imagesc(
% %         plot(0,0,'ro','MarkerSize',60)
% % an = annotation('ellipse');
% % an.Color = 'red';
% % an.FaceColor = 'r';
% % an.Position = [.45 .45 .1 .1]; 
% % an.Units = 'points';
% 
% %         ellipsev2(12,12,0,0,0,'r')
%         plot(x0,y0,'bl','DisplayName','Light Off')
%         plot(x1,y1,'r','DisplayName','Light On')
%         hold off
%         legend('Location','Best')
%         title('Track')
%         xlim([-30 30])
%         ylim([-30 30])
% %         print([currFolder,'\Figures\',nClip{c},'Track'],'-dpdf')
%         
%     end
    
   
% rmdir([currFolder,'\clips'],'s')








%%

    
%%


