function [is,insta,zone,d,Centroid] = Run_realtime_opto(frame,thresh,is,bkg,zone,blobAnalyser,d,Centroid,obsvWindow,FPS,daMask)
    global i
    frame_bkgSub = uint8(frame - bkg);
    
    frame_singBand = ceil(mean(frame_bkgSub,3));
    
    frame_bi = imbinarize(frame_singBand,thresh.bi);
    
    [insta.Area,insta.Centroid,insta.Orientation] = step(blobAnalyser,or(frame_bi,daMask.Outside));

    if length(insta.Orientation)==1%only 1 blob seen
        is.Lost=0;
        Centroid(i,:) = insta.Centroid*zone.pix2mm - zone.center;
        Orientation(1) = insta.Orientation;
        
        d(i) = sqrt(Centroid(i,1)^2 + Centroid(i,2)^2);
        insta.speedAproach = (d(i-1) - d(i))*FPS;
    
        %speed/thrust/slip/yaw
        insta.velocity(:) = (Centroid(i,:) - Centroid(i-1,:))*FPS;
        insta.speed = sqrt(insta.velocity(1,1)^2 + insta.velocity(1,2)^2);

        %determine fly kinematic state
        if insta.speed < thresh.stop
            is.Moving = 0;
            is.Approaching = 0;
        end
        if insta.speed > thresh.walkSlowL
            is.Moving = 1;
            if insta.speedAproach/insta.speed >thresh.AppOverSpd
                is.Approaching = 1;
            else
                is.Approaching = 0;
            end
        end
    end
    
    %fly positional state
    if d(i)>zone.vidRadius%outside of zone
        is.inVid = 0;
    elseif d(i)<zone.vidRadius%inside zone
        is.inVid = 1;
        if d(i)<zone.lightRadiusMM
            is.inLight = 1;
        elseif d(i)>zone.lightRadiusMM
            is.inLight = 0;
        end
    end
    
    
    
    
