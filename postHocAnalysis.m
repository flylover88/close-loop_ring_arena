function [s] = postHocAnalysis(badFrames,s,sArena,len)
%postHocAnalysis Makes corrections to the image processing done in
%FlyVideoAnalysis.m and generate a processed video of the fly kinematics 
%and movment and a structure detailing the results of the image processing.

%Thrust slip speed and yaw are in mm/s
%angs are in degrees

% 2017, Liangyu Tao

%This is an editted version of Liangyu's 
close all
angVec=s.AngVec;
f=sArena.f;
cF=sArena.cF;
badFrameNdx = find(badFrames);
for i = 1:length(badFrameNdx)
    j = badFrameNdx(i);
    if j == 1
        s.Center.x(j) = s.Center.x(j+1);
        s.Center.y(j) = s.Center.y(j+1);
        s.Head.x(j) = s.Head.x(j+1);
        s.Head.y(j) = s.Head.y(j+1);
        angVec(j) = angVec(j+1);
    elseif j == len
        s.Center.x(j) = s.Center.x(j-1);
        s.Center.y(j) = s.Center.y(j-1);
        s.Head.x(j) = s.Head.x(j-1);
        s.Head.y(j) = s.Head.y(j-1);
        angVec(j) = angVec(j-1);
    else
        s.Center.x(j) = mean([s.Center.x(j-1),s.Center.x(j+1)]);
        s.Center.y(j) = mean([s.Center.y(j-1),s.Center.y(j+1)]);
        s.Head.x(j) = mean([s.Head.x(j-1),s.Head.x(j+1)]);
        s.Head.y(j) = mean([s.Head.y(j-1),s.Head.y(j+1)]);
        angVec(j) = mean([angVec(j-1),angVec(j+1)]);
    end
end
Xchg=diff(s.Center.x);
Ychg=diff(s.Center.y);
Achg=diff(angVec);

XNdx=find(abs(Xchg)>5);
yNdx=find(abs(Ychg)>5);
ANdx=find(abs(Achg)>30);

Ndx=unique(sort([XNdx yNdx ANdx]));

if min(Ndx) == 1
    [~,ii] = min(Ndx)
    Ndx(ii) = [];
end

s.Center.x(Ndx+1) = s.Center.x(Ndx)+Xchg(Ndx-1);
s.Center.y(Ndx+1) = s.Center.y(Ndx)+Ychg(Ndx-1);
angVec(Ndx+1) = angVec(Ndx-1)+Achg(Ndx-1);
% time = 0:1/f:len/f-1/f;

%calculating displacement, thrust, slip
dspmt = sqrt(diff(s.Center.x(1:len)).^2 + diff(s.Center.y(1:len)).^2);
dspmt(dspmt>50) = nan;
speed = dspmt.*f;    %note that speed in this case is velocity
mvmtAngle = angVec(1:len-1) - unwrap(myatan(diff(s.Center.x(1:len)),...
    diff(s.Center.y(1:len)),'degrees',2));
thrustV = dspmt.*cosd(mvmtAngle)*f;
slipV = dspmt.*sind(mvmtAngle)*f;
mvmtAngle(mvmtAngle>180) = mvmtAngle(mvmtAngle>180)-360;
mvmtAngle(mvmtAngle<-180) = mvmtAngle(mvmtAngle<-180)+360;

%correct for appropiate direction
% critPoints = [];
idx = [find(s.Flags.ratio>0) find(s.Flags.size>0)];
idx = unique(sort(idx));

%find critical points where the general shape of fly is compromised
%(usually boundary points)
if ~isempty(idx)
    critPoints = [1,idx(diff(idx)>10),idx(end),len-1];
    all = 0;
else
    critPoints = [1,len-1];
    all = 1;
end
critPoints = [critPoints find(speed>50)];
critPoints(critPoints>len-1) = [];
critPoints = unique(critPoints);

%between these critical points, check to see that the general trend of
%velocity and speed matches
[angVecCor] = VelCheck(speed,thrustV,angVec,critPoints,all);

%correct for locations where the fly's yaw is consecutively above 90
%degrees.
dthe = diff(angVecCor);
dthe(dthe>180) = dthe(dthe>180) - 360;
dthe(dthe<-180) = dthe(dthe<-180) + 360;
idxTemp = find(abs(dthe)>90);
initLoc=find(diff(idxTemp)==1);
for j = 1:length(initLoc)
    angVecCor(idxTemp(initLoc(j)+1)) = angVecCor(idxTemp(initLoc(j)+1))+180;
end
angVecCor(angVecCor>360) = angVecCor(angVecCor>360)-360;

yaw = diff(angVecCor);
yaw(yaw>180) = yaw(yaw>180)-360;
yaw(yaw<-180) = yaw(yaw<-180)+360;
if abs(yaw(1)) > 100
    angVecCor(1) = angVecCor(1)+180;
end
angVecCor(angVecCor>360) = angVecCor(angVecCor>360)-360;

%correct for locations where large yall is found twice very close
% correct by thrust and change in yaw
yaw = diff(angVecCor);
yaw(yaw>180) = yaw(yaw>180)-360;
yaw(yaw<-180) = yaw(yaw<-180)+360;
ndx = find(abs(yaw(1:end-1)-yaw(2:end))>90);
ndx(diff(ndx)==1) = [];
ndx(ndx == len) = len-1; ndx(ndx == 1) = 2;
if ~isempty(ndx)
    for i = 1:length(ndx)
        case1 = [angVecCor(ndx(i))-angVecCor(ndx(i)-1),angVecCor(ndx(i))-angVecCor(ndx(i)+1)];
        case1(case1>180) = case1(case1>180)-360;
        case1(case1<-180) = case1(case1<-180)+360;
        angVecCorTemp = angVecCor(ndx(i))+180;
        angVecCorTemp(angVecCorTemp>360) = angVecCorTemp-360;
        case2 = [angVecCorTemp-angVecCor(ndx(i)-1),angVecCorTemp-angVecCor(ndx(i)+1)];
        case2(case2>180) = case2(case2>180)-360;
        case2(case2<-180) = case2(case2<-180)+360;
        if max(abs(case1))>max(abs(case2))
            angVecCor(ndx(i)) = angVecCorTemp;
        end
    end
end
angVecCor(angVecCor>360) = angVecCor(angVecCor>360)-360;

%correct for flags
yaw = diff(angVecCor);
yaw(yaw>180) = yaw(yaw>180)-360;
yaw(yaw<-180) = yaw(yaw<-180)+360;
ndx = find(abs(yaw)>90);
for i = 1:length(ndx)
    if s.Flags.ratio(ndx(i))>0
        case1 = [angVecCor(ndx(i)-1), angVecCor(ndx(i))+90, angVecCor(ndx(i)+1)];
        case1(case1>360) = case1(case1>360)-360;
        case1 = abs(diff(case1));
        case1(case1>180) = case1(case1>180)-360;
        case2 = [angVecCor(ndx(i)-1), angVecCor(ndx(i))-90, angVecCor(ndx(i)+1)];
        case2(case2<0) = case2(case2<0)+360;
        case2 = abs(diff(case2));
        case2(case2>180) = case2(case2>180)-360;
        
        case3 = [angVecCor(ndx(i)-1), angVecCor(ndx(i)), angVecCor(ndx(i)+1)];
        case3 = abs(diff(case3));
        case3(case3>180) = case3(case3>180)-360;
        
        case1 = max(abs(case1));
        case2 = max(abs(case2));
        case3 = max(abs(case3));
        caseAll = [case1,case2,case3];
        ndx1 = find(caseAll==min(caseAll));
        if ndx1==1
            angVecCor(ndx(i)) = angVecCor(ndx(i))+90;
        elseif ndx1==2
            angVecCor(ndx(i)) = angVecCor(ndx(i))-90;
        elseif ndx1==3
            angVecCor(ndx(i)) = angVecCor(ndx(i));
        end
        
        if max(abs(case1))>max(abs(case2)) && max(abs(case2))<max(abs(case3))
            angVecCor(ndx(i)) = angVecCor(ndx(i))-90;
        elseif max(abs(case1))<max(abs(case2)) && max(abs(case1))<max(abs(case3))
            angVecCor(ndx(i)) = angVecCor(ndx(i))+90;
        end
    end
end
angVecCor(angVecCor>360) = angVecCor(angVecCor>360)-360;
angVecCor(angVecCor<0) = angVecCor(angVecCor<0)+360;

tempAng = unwrap(deg2rad(angVecCor));
fineAng = smooth(tempAng,10,'loess')';
roughAng = smooth(tempAng,40,'loess')';
roughAng(abs(roughAng-tempAng)>0.15) = fineAng(abs(roughAng-tempAng)>0.15);
angVecCor = rad2deg(wrapTo2Pi(roughAng));

delX = (s.Center.x-sArena.arenaCent(2)*cF);
delY = (s.Center.y-sArena.arenaCent(1)*cF);
R = sqrt(delX.^2+delY.^2);

s.Distances.xDist = delX;
s.Distances.yDist = delY;
s.Distances.DistanceR = R;

xTemp = s.Center.x;s.Center.xUS = xTemp;
yTemp = s.Center.y;s.Center.yUS = yTemp;
fineX = smooth(xTemp,6,'loess')';
roughX = smooth(xTemp,30,'loess')';
fineY = smooth(yTemp,6,'loess')';
roughY = smooth(yTemp,30,'loess')';
temp3 = abs(roughX-xTemp);
roughX(temp3>1.5)=fineX(temp3>1.5);
s.Center.x = roughX;
temp3 = abs(roughY-yTemp);
roughY(temp3>1.5)=fineY(temp3>1.5);
s.Center.y = roughY;

s.Loss.x = sum((xTemp-s.Center.x).^2);
s.Loss.y = sum((yTemp-s.Center.y).^2);

%calculating thrust, slip, yaw
mvmtAngle = angVecCor(1:len-1) - myatan(diff(s.Center.x(1:len)),...
    diff(s.Center.y(1:len)),'degrees',2);
mvmtAngle(mvmtAngle>180) = mvmtAngle(mvmtAngle>180)-360;
mvmtAngle(mvmtAngle<-180) = mvmtAngle(mvmtAngle<-180)+360;

yaw = diff(angVecCor);
yaw(yaw>180) = yaw(yaw>180)-360;
yaw(yaw<-180) = yaw(yaw<-180)+360;
thrustV = dspmt.*cosd(unwrap(mvmtAngle))*f;
thrustV(abs(thrustV)>70) = NaN;
slipV = dspmt.*sind(mvmtAngle)*f;
slipV(abs(slipV)>50) =NaN;

s.Kinematics.thrust = thrustV;
s.Kinematics.slip = slipV;
s.Kinematics.yaw = yaw;
s.AngVec = angVecCor;

end