x = linspace(-0.01, 0.01, 1200); %distancias en metros, 1024 pix in X
T = 0.8E-3 %Periodo de la rejilla de difraccion
y = cos(2*pi*(1/T)*x);
y = y >= 0;
y = repmat(y,1200,1);
imshow(y); 