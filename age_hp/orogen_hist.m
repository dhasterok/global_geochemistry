function [S,C,age] = orogen_hist(agemin,agemax);

age = [agemin:1:agemax+1];
list = readtable('Condie&Aster2013_orogens.xlsx');

s = zeros([height(list) length(age)]);
c = zeros([height(list) length(age)]);

for i = 1:height(list)
    inds = find(list.subduction_onset(i) - list.subduction_duration(i) <= age ...
        & age <= list.subduction_onset(i));
    
    indc = find(list.collision_onset(i) - list.collision_duration(i) <= age ...
        & age <= list.collision_onset(i));
    
    s(i,inds) = 1;
    c(i,indc) = 1;
end

S = sum(s,1);
C = sum(c,1);

t = [age(1:end-1); age(2:end)]-(age(2)-age(1))/2;
t = [t(:); flipud(t(:))];

n = [C(1:end-1); C(1:end-1)];
m = -fliplr([S(1:end-1); S(1:end-1)]);
N = [n(:); m(:)];

fill(t,N,[0.7 0.7 0.7]);
hold on;
xlim([agemin agemax]);
ylim([-25 25]);
set(gca,'XTick',[0:200:4000],'YTick',[]);
ax = get(gca);
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';

plot([agemin agemax],[0 0],'k-');
for i = 0:200:4000
    if mod(i,1000) == 0
        plot([i i],[-1 1],'k-');
    else
        plot([i i],[-0.75 0.75],'k-');
    end
end
th = text(agemin,2,'Collision');
th.Rotation = 90;
th.VerticalAlignment = 'bottom';
th.FontSize = 8;

th = text(agemin,-2,'Subduction');
th.Rotation = 90;
th.VerticalAlignment = 'bottom';
th.HorizontalAlignment = 'right';
th.FontSize = 8;

sc = [335 173; 750 1071; 1500 1800; 2450 2720];
txt = {'Pangea', 'Rodina', 'Nuna', 'Kenorland'};
plot(sc',20*ones(size(sc')),'k-','LineWidth',3);
for i = 1:4
    th = text(midpt(sc(i,:)),22,txt{i});
    th.HorizontalAlignment = 'center';
    th.VerticalAlignment = 'bottom';
end
hold off 

S = S(:);
C = C(:);
age = age(:);

return