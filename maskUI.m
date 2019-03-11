function data = maskUI(frame,message)
   figure(1);
   imshow(uint8(frame(:,:,1)));
   title(message)
   h=imellipse; %%%%%%HPARENT IS NOT DEFINED imellipse(hparent,position)
   wait(h);
%   vertices=wait(h);
   vertices=getVertices(h); %obtains vertices of the ellipse "h"
   Inside = h.createMask();
   Outside=~Inside;
%    Inside=+Inside;%convert each frame to a double and multiply.
   data.Inside=Inside;
%    Outside=+Outside;
   data.Outside=Outside;
%    data.vertices=vertices;
   data.Inside=Inside;
   data.Vertices= vertices;
   data.Pos = h.getPosition;
   close(1)
   