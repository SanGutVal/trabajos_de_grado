function y = rejillaDifraccion(resolucion,X)
x1 = resolucion(1), y1 = resolucion(2);
%x = linspace(-0.01, 0.01, resolucion);
x = linspace(0, 1, x1);
%X = 30;
y = 255 + 255*cos(2*pi*(X)*x);
%y = y >= 0;
y = uint8(y);
y = repmat(y, y1, 1);
imshow(y)
%imwrite(y, 'rejilla.tif')

%figure(1), imagesc(rejillaDifraccion(480,30)), truesize;