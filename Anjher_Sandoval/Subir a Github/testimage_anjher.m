N = 300;
[x,y] = meshgrid(linspace(-1,1,N));
testimage = cos(2*pi*(10*y));
figure(1), imagesc(testimage), truesize;
TESTIMAGE = fftshift(fft2(testimage));
figure(2), imagesc(abs(TESTIMAGE)), truesize
