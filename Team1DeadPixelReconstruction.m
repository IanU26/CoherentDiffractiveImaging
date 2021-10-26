%This load the diffraction pattern with missing information at the center
%as well as the location of the dead cell. 
load('./diff_pat_dead_pixels.mat')
load('./dead_pixels.mat')
%load('./diff_pat_team_1.mat')
%dead_pixels is a 2D matrix that equals 1 where the diffraction pattern is
%missing. 

diff_pat_dead_pixels = diff_pat_dead_pixels;
%diff_pat = diff_pat * pi;

%The high and low are calculated here. If changing the support size, one
%must change it here and in the functions!
supportsize = 160;
idx_low = (size(diff_pat_dead_pixels,2))/2-supportsize+1; 
idx_high = (size(diff_pat_dead_pixels,2))/2+supportsize; 
idy_low = (size(diff_pat_dead_pixels,2))/2-supportsize+1; 
idy_high = (size(diff_pat_dead_pixels,2))/2+supportsize;


support = zeros(size(diff_pat_dead_pixels,2),size(diff_pat_dead_pixels,2));
support(idx_low:idx_high, idy_low:idy_high) = 1;




run_hio(diff_pat_dead_pixels,support,dead_pixels)

%Now we must include the location of the dead_pixel along with the
%diffraction pattern and the support. 
function [recon, x] = run_hio(diff_pat, support, dead_pixels)
supportsize = 160;
idx_low = (size(diff_pat,2))/2-supportsize+1; 
idx_high = (size(diff_pat,2))/2+supportsize; 
idy_low = (size(diff_pat,2))/2-supportsize+1; 
idy_high = (size(diff_pat,2))/2+supportsize;

rfact = [];
obj = rand(size(diff_pat,2),size(diff_pat,2))*0.1; %This will be our 0th iteration of the algorithm. Random array to start with. 
beta = 0.40;

for iter = 1:500
    dp = fftshift(fft2(ifftshift(obj))); %fourier trans of 0th iteration 
    R2 = abs(dp); %This will be used for R-Factor Calc
    
    %This is the key difference when there are dead_pixels in the
    %diffraction pattern. Essentially what we are doing is accepting there
    %is no possible way to obtain phase information from the pixels that
    %have no value so we ignore them. 
    for i = 1:size(dp,2)
        for j = 1:size(dp,2)
            if(dead_pixels(i,j) == 1)
                dp(i,j) = dp(i,j);
            else 
                dp(i,j) = diff_pat(i,j)*exp(1j*angle(dp(i,j))); %Now we set the dp to the diff_pat * phase of dp
            end 
        end
    end 
    
    
    
    obj_temp = fftshift(ifft2(ifftshift(dp))); %Send it back to real space

    s = logical((support.*(imag(obj_temp)>0))>0); %Apply the constraint

    for i = 1:size(s,2)
        for j = 1:size(s,2)

            if(s(i,j)==1)
                obj_new(i,j)= obj_temp(i,j);
            else 
                obj_new(i,j) = obj(i,j) - beta*obj_temp(i,j);
            end 
        end
    end 
    
    %So obj will be purely phase
    obj_new(idx_low:idx_high,idy_low:idy_high) = obj_new(idx_low:idx_high,idy_low:idy_high)./abs(obj_new(idx_low:idx_high,idy_low:idy_high));
    obj = obj_new;
    
    plot_results(diff_pat,obj);
    

end
recon = obj;
end


function plot_results(diff_pat, recon)
supportsize = 160;
idx_low = (size(diff_pat,2))/2-supportsize+1; 
idx_high = (size(diff_pat,2))/2+supportsize; 
idy_low = (size(diff_pat,2))/2-supportsize+1; 
idy_high = (size(diff_pat,2))/2+supportsize;

a = angle(recon(idx_low: idx_high,idy_low: idy_high));
a = -1*flip(a);
a = a - mean(a);

imagesc(a);
end

%{
%R factor plotting function 
function plot_rfactor(x)
y = [];
for i = 1:size(x,2)
    y(i) = i;
end
plot(y,x);
end
%}




%The code below is used to determine the support size
%{
load('./diff_pat_dead_pixels.mat')
supportsize = 224;
idx_low = (size(diff_pat_dead_pixels,2))/2-supportsize+1; 
idx_high = (size(diff_pat_dead_pixels,2))/2+supportsize; 
idy_low = (size(diff_pat_dead_pixels,2))/2-supportsize+1; 
idy_high = (size(diff_pat_dead_pixels,2))/2+supportsize;


support = zeros(size(diff_pat_dead_pixels,2),size(diff_pat_dead_pixels,2));
support(idx_low:idx_high, idy_low:idy_high) = 1;
fftshift(ifft2(ifftshift(support)));

imagesc(abs(fftshift(ifft2(ifftshift(support)))))

%}

