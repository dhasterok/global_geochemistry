%----------------------------
% Generate an example HP plot
%----------------------------

% Explain why a HP vs age might be flat when we know things like depletion
% and cooling/lower melt fractions have occured

% Doesnt account for recycling, and just a concept - no definitive models
% output here



figure()

age_product = [0:1:4000];
age_source = age_product;

hp_product = 3.*ones(size(age_product));
hp_source = 1.*ones(size(age_product));

% Flat source, flat product - this assumed constantly replenished source,
% no temperature change, no recycling. Products created today are same as
% past, same rock, and only thing influencing it is decay.
plot(age_product,hp_product,'-b')
hold on
plot(age_source,hp_source,'-r')
hold off

% However we know the Earth has cooled - and can assume that average melt
% fraction has also decreased with time.


%y = 4, x = 4000
%y =  3, x = 0
% y = 3.25 x = 3500

hp_product_mf = hp_product .* exp(0.05.*age_product./800);
hold on
plot(age_product,hp_product_mf,'--b')
hold off





ylim([0 5])
xlim([0 4000])