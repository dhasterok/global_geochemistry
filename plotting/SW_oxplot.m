function oxplot(data,varargin)
% oxplot(data)

% igneous axes (default)
ax = [30 100; 0 5; 0 30; 0 20; 0 20; 0 20; 0 10; 0 10; 0 2.5];
if nargin == 2
    if varargin{1} == 1
        % sedimentary (axes)
        ax = [30 100; 0 3; 0 30; 0 20; 0 20; 0 80; 0 10; 0 10; 0 2.5];
    end
end

%percentage plots
oxlist = {'SiO2','Na2O','K2O'};

figure;
for i = 1:length(oxlist)
    if i ~= 1
        harker(data,oxlist{2},oxlist{3},ax(1,:),ax(i,:));
    end
end

return


function harker(data,el1,el2,xl,yl);

sio2 = data.sio2;
d = data{:,lower(el1)}+data{:,lower(el2)};
scatter(data.sio2(:),d(:),10,'*');
set(gca,'CLim',[-2 2]);
xlabel(['SiO2 (wt. %)']);
ylabel(['Na2O + K2O',' (wt. %)']);
xlim(xl);
ylim([0,15]);

return
