%This script is used to show what the diffraction pattern looks like and
%how it is completely unrecognizable without the HIO algorithm. 
load('./diff_pat_team_1.mat')
imagesc(abs(diff_pat))