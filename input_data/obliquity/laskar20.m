
clear all
ts = [50,100,250,500,1000,2000,3000,4000];
mus = [34.244,34.068,34.432,35.082,35.950,36.779,37.267,37.620];
sigmas = [7.149,8.438,9.979,10.977,11.987,12.978,13.499,13.814];

x = 0:.1:80;
y = zeros(length(ts),length(x));
colorx = ["b","g","r","c","m","k","y","b"];

for i = 1:1:length(ts)
    y(i,:) = normpdf(x,mus(i),sigmas(i));
    plot(x,y(i,:),colorx(i));
    hold on
end

laskar_gaussian.x = x;
laskar_gaussian.pdfs = y;
laskar_gaussian.times = ts;
laskar_gaussian.mus = mus;
laskar_gaussian.sigmas = sigmas;