load('./yeast_diffpat.mat')

%The high and low are calculated here. If changing the support size, one
%must change it here and in the functions!
%The true support size will be the input*2 because it is added AND
%subtracted by the midpoint
xsupportsize = 98;
ysupportsize = 98;
idx_low = (size(diffpat1,2))/2-xsupportsize+1; 
idx_high = (size(diffpat1,2))/2+xsupportsize; 
idy_low = (size(diffpat1,2))/2-ysupportsize+1; 
idy_high = (size(diffpat1,2))/2+ysupportsize;

support = zeros(size(diffpat1,2),size(diffpat1,2));
support(idx_low:idx_high, idy_low:idy_high) = 1;

%Now we finally call the function run_hio with the only needed inputs being
%the diffraction pattern and a support
run_hio(diffpat1,support)


function [recon, x] = run_hio(diff_pat, support)

xsupportsize = 98;
ysupportsize = 98;
idx_low = (size(diff_pat,2))/2-xsupportsize+1; 
idx_high = (size(diff_pat,2))/2+xsupportsize; 
idy_low = (size(diff_pat,2))/2-ysupportsize+1; 
idy_high = (size(diff_pat,2))/2+ysupportsize;


obj = rand(size(diff_pat,2),size(diff_pat,2)); %This will be our 0th iteration of the algorithm. Random array to start with. 
beta = 0.10;
%ERF is used for ERROR 
erf = [];

for iter = 1:2000
   
    iter
    dp = fftshift(fft2(ifftshift(obj))); %fourier trans of 0th iteration 
    
    
    %This portion handles the beamstop. The beamstop is located at dp = -1
    for i = 1:size(dp,2)
        for j = 1:size(dp,2)
            if(diff_pat(i,j) == -1)
                dp(i,j) = dp(i,j);
            else 
                dp(i,j) = diff_pat(i,j)*exp(1j*angle(dp(i,j))); %Now we set the dp to the diff_pat * phase of dp
            end 
        end
    end 
    
    
    
    obj_temp = fftshift(ifft2(ifftshift(dp))); %Send it back to real space
    
    %s=1 represents the pixels inside the support with real-positive values
    s = logical((support.*(real(obj_temp)> 0))>0); %Apply the constraint

    %r = 1 represents pixels just inside the support. USED IN ERROR 
    r = logical((support)>0);
    
    
    %HIO positively updating the pixels inside the support 
    %HIO negatively updating the pixels outside the support or with
    %negative real values.
    for i = 1:size(s,2)
        for j = 1:size(s,2)

            if(s(i,j)==1)
                obj_new(i,j)= obj_temp(i,j);
            else 
                obj_new(i,j) = obj(i,j) - beta*obj_temp(i,j);
            end 
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The code below is for the ERROR 
    %This step is to calculate how well we are recreating the image. 
    %ERF is similar to the R factor we used in Team1DiffPatt
    for i = 1:size(r,2)
        for j = 1:size(r,2)

            if(r(i,j)==1)
                sig(i,j)= abs(obj_new(i,j));
            else 
                sig(i,j) = 0;
            end 
        end
    end
    
    
    H = fftshift(fft2(ifftshift(sig)));
    
    %erf = (abs(sum(diff_pat)-abs(sum(H)))/(abs(sum(diff_pat))));
    num = 0;
    denom = 0;
    
    for i = 1:size(dp,2)
        for j = 1:size(dp,2)
            if(diff_pat(i,j) == -1)
                continue
            else 
                num = num + abs( (abs( diff_pat(i,j) ) - abs( H(i,j) )) );
                denom = denom + (abs(diff_pat(i,j)));
            end 
        end
    end 
    erf(end+1) = num/denom;
    %This creates a 1D array of errors. Each input is an iteration
    %End of error stuff 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %So obj will be purely real
    obj_new(idx_low:idx_high,idy_low:idy_high) = real(obj_new(idx_low:idx_high,idy_low:idy_high));
    
    obj = obj_new;
    
    plot_results(diff_pat,obj);
    %plot_erf(erf);

end
recon = obj;
x = erf;
end


function plot_results(diff_pat, recon)
xsupportsize = 98;
ysupportsize = 98;
idx_low = (size(diff_pat,2))/2-xsupportsize+1; 
idx_high = (size(diff_pat,2))/2+xsupportsize; 
idy_low = (size(diff_pat,2))/2-ysupportsize+1; 
idy_high = (size(diff_pat,2))/2+ysupportsize;


a = real(recon(idx_low:idx_high, idy_low:idy_high));
%a = -1*flip(a);
a = a - mean(a);

imagesc(abs(a));
end


%This function will plot the erf 
function plot_erf(x)
y = [];
for i = 1:size(x,2)
    y(i) = i;
end
plot(y,x);
end
