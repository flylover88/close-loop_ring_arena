function [angVecCor] = VelCheck(speed,thrustV,angVec,critPoints,all)
angVecCor = angVec;

%between these critical points, check to see that the general trend of
%velocity and speed matches
for j = 1:length(critPoints)-1
    case1 = length(find(thrustV(critPoints(j)+1:critPoints(j+1))>0))>...
        length(find(thrustV(critPoints(j)+1:critPoints(j+1))<0));
    case2 = length(find(speed(critPoints(j)+1:critPoints(j+1))>0))>...
        length(find(speed(critPoints(j)+1:critPoints(j+1))<0));

    if(case1 && ~case2 || ~case1 && case2)
        if all == 1
            angVecCor = angVecCor + 180;
        else
            angVecCorTemp = cell(2,23);
            cost = ones(23,23)*1000;
            track = cell(23,23);
            for init = -11:11
                for ending = -11:11
                    track{init+12,ending+12} = max(critPoints(j)+init,1):min(critPoints(j+1)+ending,length(angVec));
                    angVecCorTemp{init+12,ending+12} = angVecCor(track{init+12,ending+12});
                    angVecCorTemp{init+12,ending+12}(2:end-1) = ...
                        angVecCorTemp{init+12,ending+12}(2:end-1)+180;

                    yaw = diff(angVecCorTemp{init+12,ending+12});
                    yaw(yaw>180) = yaw(yaw>180)-360;
                    yaw(yaw<-180) = yaw(yaw<-180)+360;
                    if ~isempty(yaw)
                        cost(init+12,ending+12) = max(abs(yaw));
                    end
                end
            end
            [row,col] = find(abs(cost)==min(min(abs(cost))));
            angVecCor(track{row(1),col(1)})=angVecCorTemp{row(1),col(1)};
%           angVecCor(critPoints(j)+1:critPoints(j+1)) = angVec(critPoints(j)+1:critPoints(j+1))+180;
        end
    end
end            
end