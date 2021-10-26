%This section will load in the diffraction pattern of the image I am hoping
%to reconstruct as well as a support. The support is a matrix of zeros with
%a square of 1's at the center. 
load('./diff_pat_team_1.mat')
load('./support.mat')
diff_pat = diff_pat * pi;

%The high and low are calculated. These will be used to set the center of 
%the support to 1 based on the size of the diffraction pattern. 
idx_low = 1024-128+1; 
idx_high = 1024+128; 
idy_low = 1024-128+1; 
idy_high = 1024+128;

%Setting the center of the support to 1. 
support = zeros(size(diff_pat,2),size(diff_pat,2));
support(idx_low:idx_high, idy_low:idy_high) = 1;

%Running the HIO Algorithm with diff_pat and support above.
run_hio(diff_pat,support)

%Beginning of Hybrid Input Output Algorithm 
function [recon, x] = run_hio(diff_pat, support)
idx_low = 1024-128+1; 
idx_high = 1024+128; 
idy_low = 1024-128+1; 
idy_high = 1024+128;
rfact = [];
%This will be our 0th iteration of the algorithm. 
%We begin with a completely random matrix the size of the diffraction pattern.
obj = rand(size(diff_pat,2),size(diff_pat,2))*0.1;  
beta = 0.10;

%We will iterate 500 times to reconstruct the original image.
for iter = 1:500
    iter
    dp = fftshift(fft2(ifftshift(obj))); %fourier transformation of 0th iteration 
    R2 = abs(dp); %This will be used for R-Factor Calc
    
    dp = diff_pat.*exp(1j*angle(dp)); %Now we set the dp to the diff_pat * phase of dp. 
    
    obj_temp = fftshift(ifft2(ifftshift(dp))); %Send it back to real space with an inverse Fourier Transformation

    s = logical((support.*(imag(obj_temp)>0))>0); 
    
    %This nested loop is crucial because it updates our reconstruction IFF
    %the support==1 (at the center) AND IFF the imaginary part is positive.
    %If these constraints aren't met the reconstructed image is "taxed" by
    %the -beta*obj_temp(i,j). This effectively forces those values to zero
    %throughout the iterations. 
    for i = 1:size(s,2)
        for j = 1:size(s,2)

            if(s(i,j)==1)
                obj_new(i,j)= obj_temp(i,j);
            else 
                obj_new(i,j) = obj(i,j) - beta*obj_temp(i,j);
            end 
        end
    end
 
    
    %This forces obj to be a purely phase object.
    obj_new(idx_low:idx_high,idy_low:idy_high) = obj_new(idx_low:idx_high,idy_low:idy_high)./abs(obj_new(idx_low:idx_high,idy_low:idy_high));
    obj = obj_new;
    
    %R Factor 
    %Sums all the pixel values for the abs(original Diff pattern) 
    %Whichever plot_results / plot_rfactor is below will be what is output
    %R1 = abs(fftshift(fft2(ifftshift(obj))));
    %r = ( abs (sum(abs(R2),'all') - abs(sum(R1,'all'))) ) / (sum(abs(R2),'all'));
    %rfact(size(rfact,2)+1) = r;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    plot_results(diff_pat,obj);
    %plot_rfactor(rfact);    
    
    
end
recon = obj;
x = rfact;
end 

%This function will plot the image we reconstructed in the HIO algorithm.
function plot_results(diff_pat, recon)
idx_low = 1024-128+1; 
idx_high = 1024+128; 
idy_low = 1024-128+1; 
idy_high = 1024+128;

a = angle(recon(idx_low: idx_high,idy_low: idy_high));
a = -1*flip(a);
a = a - mean(a);
imagesc(abs(diff_pat));
imagesc(a);
end

%R factor plotting function. R factor is a measure of how well the image
%was reconstructed.
function plot_rfactor(x)
y = [];
for i = 1:size(x,2)
    y(i) = i;
end
plot(y,x);
end





