% clear all; close all; daqreset;
imaqreset;
close all
%% Overall, we want to acquire a 6 minute video. During the last three minutes we will
%stimulate with red light to activate ORNs. We have to acquire video and in
%sync stimulate.
cd('F:\CLoptogenetics')
currFolder = pwd;
addpath([currFolder '\utils'])
addpath([currFolder '\process'])
addpath([currFolder '\Analysis'])

global i LED_in t
%% User inputs
mask_make = 0; 
test = 0;
veiw = 0;
count_runsMax = 20;%number of runs wanted in this trial
d_passThru = 10;%minimum dist (mm) needed for a run to be "good"
pauseTime = 15*60;
numSecs=30*60;%%%%%%%%%%%%%%%%%%%%%%Change test duration

%% basic parameters
sampleRate=1000;
FPS=30;

frames=numSecs*FPS;

zone.fullRadiusMM=40;%mm
zone.lightRadiusMM=12.7;
zone.vidBorder = 7.5;

thresh.Buff = .05;
thresh_noiseMin=30;
thresh_noiseMax=300;

thresh.AppOverSpd = 1/sqrt(2);
thresh.MinorMajor = .9;
thresh.stop = 0.3 - thresh.Buff;%mm/s

thresh.walkSlowL = 0.3 + thresh.Buff;%??
thresh.walkSlowH = 5 - thresh.Buff;
thresh.walkFast = 5 + thresh.Buff;

thresh.bi = 30;
bkg_buildTime = 5;%frames
zone.blkFrameBuffer=10;

light_pulseWidth = 3;%seconds
Time_enterBuff = .3;
lightOnFirst_prob = .5;

% dilateSize_pix=3;
obsvWindow = frames;
%% Initialize Counters & states
zone.outCount=1;
count_runs=0;
r=0;

LED_old = -1;
% run = 1;
% firstFrameIn = 1;

is.Moving = 0;
is.Approaching = 0;
is.inVid = 0;
is.inLight =0;
%% Connect Video Source
if test
    v_in = 'C:\Users\tjd37\Desktop\RecVid\180809_1_20_video.avi';
    vid = VideoReader(v_in);
    
    Res = [vid.Width vid.Height];
    numBands=3;
    
else ~test;
    %initialize camera
%     vid = videoinput('avtmatlabadaptor64_r2009b', 2, 'F7M3_Mono8_640x512');
    vid = videoinput('gentl', 1, 'Mono8');
    src = getselectedsource(vid);
    %F7M3_YUV422_640x510'); %'F0M5_Mono8_640x480');
    % vid=videoinput('avtmatlabadaptor64_r2009b', 1, 'Mono8');

    Res = vid.VideoResolution;
    numBands = vid.NumberOfBands;
    %Set Camera Parameters
    vid.TriggerRepeat=frames+bkg_buildTime+mask_make-1; %handles.numframes;
    vid.FramesPerTrigger=1;
    vid.LoggingMode = 'Memory';
    triggerconfig(vid,'immediate')
    %triggerconfig(vid,'manual')%%%ADDED 11/29/2018
    %src.AcquisitionFrameRate = FPS;
disp('Camera Connected')
    % %preview video
    % preview(vid);
    % 
    % % waitforbuttonpress();
    % closepreview(vid);
end


%% Initialize Matrices
frame_timeStamp = NaN(frames,1);
LED_in = NaN(frames,2);
LEDon = zeros(frames,1);
frame_first = zeros(Res(2),Res(1),numBands,bkg_buildTime);
frame_noFly = zeros(Res(2),Res(1));

d = NaN(frames,1);
% Orientation = NaN(frames,1);
Centroid = NaN(frames,2);
disp('Variables & matrices mades')


%% Initialize save Location
%define where things will be saved
Date = datestr(now,'yymmdd');
time = datestr(now,'hhmmss');
if ~test
    DirLocation = strcat('F:\CLoptogeneticsT\');
else
%     DirLocation = strcat('C:\Users\tjd37\Desktop\test', Date);
end
if ~exist(DirLocation)
    mkdir(DirLocation);
end

%Find Trial Number in data.mat
FileLocation = strcat(DirLocation,'\TrialCount.mat');
if exist(FileLocation,'file')
    load(FileLocation,'Trials')
    Trials=Trials+1;
else
    Trials=1;
end

%create video file
video_num = strcat('\video',int2str(Trials));
VideoFile = strcat(DirLocation,video_num);
% diskLogger= VideoWriter([VideoFile,'.avi'], 'Grayscale AVI');
diskLogger= VideoWriter([VideoFile,'.avi'], 'Motion JPEG AVI');
v_out=diskLogger;
disp('Video File Initialized')

%% MASK MAKING & Determine Zones
if ~exist('Mask') && mask_make == 0
    disp('No mask found. mask_make changed to 1')
    mask_make = 1;
end
if ~test
    start(vid)
end
if mask_make
    if ~test
        
        frame_first(:,:,:,1)=getdata(vid);
    else
        frame_first(:,:,:,1)=read(vid,1);
    end
    message='Select rink perimeter. Double click when done.';
    Mask = maskUI(frame_first(:,:,:,1),message);
    
    zone.fullRadiusPix = mean(Mask.Pos(3:4));
    zone.pix2mm = 2*zone.fullRadiusMM/zone.fullRadiusPix;
    
    zone.vidRadius = zone.lightRadiusMM + zone.vidBorder;
    zone.center = [Mask.Pos(1)+Mask.Pos(3)/2,Mask.Pos(2)+Mask.Pos(4)/2]*zone.pix2mm;%x,y
end
if ~test
    stop(vid)
end
% pause(10)
%% BKG AND BLOB
disp(['Creating Background with first ',int2str(bkg_buildTime), ' frames'])
start(vid)
for i=mask_make+1:bkg_buildTime
    if ~test
%         trigger(vid);%ADDED 11/29/2018
        frame_first(:,:,:,i)=getdata(vid);
    else
        frame_first(:,:,:,i)=read(vid,i);
    end
end
stop(vid)
bkg = uint8(floor(mean(frame_first,4)));

blobAnalyser = vision.BlobAnalysis('AreaOutputPort',1,...
                        'BoundingBoxOutputPort',0,...
                        'OrientationOutputPort',1,...
                        'MaximumCount',1000,...
                        'MinimumBlobArea',thresh_noiseMin,...
                        'MaximumBlobArea',thresh_noiseMax*100,...
                        'ExcludeBorderBlobs',1,...
                        'MinorAxisLengthOutputPort',0,...
                        'MajorAxisLengthOutputPort',0);



%% initialize DAQ
disp('Initilizing DAQ')

s = daq.createSession('ni');
s.Rate=sampleRate;
s.IsContinuous=0;
s.TriggersPerRun = 1;

s.addAnalogOutputChannel('Dev1','ao0','voltage'); %drives the red light
s.addAnalogInputChannel( 'Dev1','ai0','voltage'); %receives the same command which drives the red light

% start the listener in the background
dataListener = s.addlistener('DataAvailable', @collectLED);
LED = 0;
%% Run Camera Loop
if ~test
    disp(['Pause for ',int2str(pauseTime),' Seconds'])
    preview(vid)
    pause(pauseTime)
end
start(vid)
open(v_out)
display('Running camera loop');
t=tic;

for i=bkg_buildTime+1:frames
    if ~test
        %trigger(vid)
        frame = getdata(vid);
        %
    else
        frame = read(vid,i);
    end
    time = toc(t);
[is,insta,zone,d,Centroid] = Run_realtime_opto(frame,thresh,is,bkg,zone,blobAnalyser,d,Centroid,obsvWindow,FPS,Mask);
%d is distance
d(i)
% is.inVid=1;
    if is.inVid
        
        writeVideo(v_out,frame(:,:,1));
        frame_timeStamp(i,1) = time;
        
       
            if zone.outCount > 0 %first frame in
                Time_in=toc(t);
                frame_in = i;
    %             count_runs = count_runs + 1;
            end
            if toc(t)-Time_in <= Time_enterBuff%LED sequence is off while the buffer time hasnt passed
                LED_new = -1;

            elseif toc(t)-Time_in > r*light_pulseWidth + Time_enterBuff
                if rand() < lightOnFirst_prob || r>1
                    LED_new = LED_old*-1;
                end
                r=r+1;
            end
            zone.outCount = 0;
        
    else
        if zone.outCount<zone.blkFrameBuffer
            %writeVideo(v_out,frame_noFly);
            %frame_timeStamp(i,1) = -time;
            zone.outCount=zone.outCount+1;
            if zone.outCount == 1
                if min(d(frame_in:i)) <= d_passThru
                    count_runs = count_runs+1;
                    disp(['Good run @',num2str(Time_in),' (',int2str(count_runs),' runs completed)'])
                end
            end
        end 
        r = 0;
        LED_new=-1;
    end
    
    if LED_new~=LED_old
        LED = (LED_new + 1)*2.5;%for on LED = 5; off =0
        LED_old = LED_new;
        s.queueOutputData(LED);
        s.startForeground();
    end
    if LED == 5
        LEDon(i,1) = 1;
    end
    
    if test && veiw
        imshow(frame)
    elseif ~test
        while toc(t)<i/FPS
            pause(.001)
        end
    end
    if count_runs >= count_runsMax
        disp(['Trial ended early, max runs (',int2str(count_runsMax),') reached'])
        toc(t)
        break
    end
end

s.queueOutputData(0)
s.startForeground()

disp(['Saving Data Files ',VideoFile,' in ',DirLocation])
close(v_out)
delete(dataListener)

if test
    disp(['Ran at ',int2str(frames/toc(t)), ' fps'])
end

save(FileLocation,'Trials')
save([VideoFile, '.mat'],'frame_timeStamp','LED_in','Mask','FPS','Date')

disp('done!')
disp([int2str(count_runs),' good run(s) collected'])
