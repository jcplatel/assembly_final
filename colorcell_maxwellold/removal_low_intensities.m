function [R,G,B,id_low_intensity] = removal_low_intensities (R,G,B,prct_R,prct_G,prct_B)
% 
bruit_R = prctile(R, prct_R);
bruit_G = prctile(G, prct_G);
bruit_B = prctile(B, prct_B);
is_valid_R = R > (bruit_R);
is_valid_G = G > (bruit_G);
is_valid_B = B > (bruit_B);

id_low_intensity = ~(is_valid_R | is_valid_G | is_valid_B);
% 
% percentile_th=10;
% I = sqrt(R.^2 + G.^2 + B.^2);
% % figure;histogram (R,20)
% % figure;histogram (G,20)
% % figure;histogram (B,20)
% thr = prctile(I, percentile_th);   % 5% les plus faibles
% id_low_intensity = (I <= thr);

nb_low_intensity = sum(id_low_intensity);
fprintf('Nombre de cellules low_intensity : %d\n', nb_low_intensity);
% On force ces cellules en blanc
%C_points(low_intensity, :) = repmat([1 1 1], sum(low_intensity), 1);


% On met les intensités de ces cellules à NaN pour qu'elles soient ignorées
R(id_low_intensity) = NaN;
G(id_low_intensity) = NaN;
B(id_low_intensity) = NaN;

end