function [flipped_image] = Flip_Image(image, direction)
%UNTITLED2 Summary of this function goes here
%this function flips an image upside down if direction=0 and flips it right
%to left if direction = 1
%   Detailed explanation goes here
[x,y] = size(image);
flipped_image = image;

if (direction == 1)
   for n=1:1:x
       for m=1:1:y
           flipped_image(n,m) = image(n,y+1-m);
       end
   end
end

if (direction == 0)
   for n=1:1:x
       for m=1:1:y
           flipped_image(n,m) = image(x+1-n,m);
       end
   end
end


end

