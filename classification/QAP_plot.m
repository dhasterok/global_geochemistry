function qap_name = QAP_plot(data,varargin)

% data = mineral_norm(data,QAP);
normtype = 'mode';
if nargin > 1
    if strcmpi(varargin{1},'norm')
        normtype = lower(varargin{2});
    end
end

volc = strcmpi(data.rock_origin,'volcanic');

switch normtype
    case 'cipw'
        cipw = cipwnorm(data);
        %cipw.Properties.VariableNames
        Q = cipw.Quartz;
        A = cipw.Orthoclase;
        P = cipw.Albite + cipw.Anorthite;
        F = cipw.Nepheline + cipw.Leucite;
        Ol = cipw.Olivine_Mg + cipw.Olivine_Fe;
        Cpx = cipw.Diopside_Mg + cipw.Diopside_Fe;
        Opx = cipw.Hypersthene_Mg + cipw.Hypersthene_Fe;
        An_index = cipw.An_index;
    case 'mode'
        Q = data.Quartz_norm;
        A = data.Orthoclase_norm;
        P = data.Albite_norm + data.Anorthite_norm;
        F = data.Leucite_norm + data.Nepheline_norm;
        Ol = data.Olivine_norm;
        Cpx = data.Diopside_norm + data.Hedenbergite_norm;
        Opx = data.Enstatite_norm + data.Ferrosillite_norm;
        An = (data.Anorthite_norm/278.21) ./ (data.Anorthite_norm/278.21 + data.Albite_norm/524.45);
    otherwise
        error('Unknown norm.');
end

qapgons = load_qapgons;
[nqap,~] = size(qapgons);

% Q = 'Quartz';
% A = 'Alkali Feldspar';
% P = 'Plagioclase';


ifoid = F > 0;

qap_name = cell(height(data),1);
qap_name(1:end) = {''};

X = zeros(size(Q));
Y = zeros(size(Q));
[X(~ifoid),Y(~ifoid)] = tern2xy(Q(~ifoid),A(~ifoid),P(~ifoid));
[X(ifoid),Y(ifoid)] = tern2xy(F(ifoid),A(ifoid),P(ifoid));
Y(ifoid) = -Y(ifoid);
for i = 1:nqap
    if sum(ifoid) == 0
        if any(qapgons{i,1}(:,1) > 0)
            [qg1,qg2] = tern2xy(qapgons{i,1}(:,1),qapgons{i,1}(:,2),qapgons{i,1}(:,3));
        else
            continue;
        end
        in = inpolygon(X,Y,qg1,qg2);
        qap_name(in) = {qapgons{i,2}};
    else
        if any(qapgons{i,1}(:,1) > 0)
            [qg1,qg2] = tern2xy(qapgons{i,1}(:,1),qapgons{i,1}(:,2),qapgons{i,1}(:,3));
        else
            [qg1,qg2] = tern2xy(qapgons{i,1}(:,4),qapgons{i,1}(:,2),qapgons{i,1}(:,3));
            qg2 = -qg2;
        end
        in = inpolygon(X,Y,qg1,qg2);
        qap_name(in & ~volc) = {qapgons{i,2}};
        
        qap_name(in & volc) = {qapgons{i,3}};
        if strcmp(qapgons{i,3},'tephrite')
            qap_name(in & volc & Ol > 10) = {'basanite'};
        elseif strcmp(qapgons{i,3},'phonolitic tephrite')
            qap_name(in & volc & Ol > 10) = {'phonolitic basanite'};
        elseif strcmp(qapgons{i,3},'andesite')
            qap_name(in & volc & data.sio2 <= 52) = {'basalt'};
        end
    end
end
mind = (Q + A + P + F) < 10;
qap_name(mind) = {''};


fprintf('Computing mafic names...\n\n');
%mafind = An_index > 50 & ~mind & Q./(Q + A + P) < 0.2 & (Q == 0 | Ol + Cpx + Opx > 40);
mafind = An_index > 50 & ~mind & Ol + Cpx + Opx > 40;
mgons = load_mgons;
[nmaf,~] = size(mgons);
[X,Y] = tern2xy(P,Opx,Cpx);
for i = 1:nmaf
    [qg1,qg2] = tern2xy(mgons{i,1}(:,1),mgons{i,1}(:,2),mgons{i,1}(:,3));
    in = inpolygon(X,Y,qg1,qg2);
    qap_name(in & mafind) = {mgons{i,2}};
end


fprintf('Computing M names...\n\n');
umgons = load_umgons;
[num,~] = size(umgons);
[X,Y] = tern2xy(Ol,Opx,Cpx);
for i = 1:num
    [qg1,qg2] = tern2xy(umgons{i,1}(:,1),umgons{i,1}(:,2),umgons{i,1}(:,3));
    in = inpolygon(X,Y,qg1,qg2);
    qap_name(in & mind) = {umgons{i,2}};
end

qap_name(data.sio2 < 20) = {'carbonatite'};
qap_name(20 < data.sio2 & data.sio2 < 33) = {'silicocarbonatite'};

%qap_name
%count and print out all the different rock_types
[type,~,ind] = unique(qap_name);
fprintf('No. Samples       QAP Name (calc. CIPW)\n');
fprintf('-----------   ----------------------------\n');
for i = 1:length(type);
    ns = length(find(ind == i));
    fprintf('%10i    %s\n',ns,type{i});
end
fprintf(' \n');

%plot the figure
figure;
subplot(121);
if sum(ifoid) == 0
    ternary('Q','A','P');
else
    ternary('Q','A','P','F');
end

% plot data
if sum(ifoid) == 0
    ternplot(Q(~mind & ~mafind),A(~mind & ~mafind),P(~mind & ~mafind),'.');
    ternplot(Q(mafind),A(mafind),P(mafind),'.');
else
    ternplot(Q(~mind & ~mafind),A(~mind & ~mafind),P(~mind & ~mafind),F(~mind & ~mafind),'*');
    ternplot(Q(mafind),A(mafind),P(mafind),F(mafind),'.');
end

% draw polygons
for i = 1:nqap    
    q = qapgons{i,1}(:,1);
    a = qapgons{i,1}(:,2);
    p = qapgons{i,1}(:,3);
    f = qapgons{i,1}(:,4);
    if any(q > 0)
        ternplot(q,a,p,'-k');
        [x,y] = tern2xy(q,a,p);
    elseif sum(ifoid) > 0
        ternplot(q,a,p,f,'-k');
        [x,y] = tern2xy(f,a,p);
        y = -y;
    else
        continue; 
    end
    t = text(mean(x),mean(y),qapgons{i,2});
    set(t,'HorizontalAlignment','center','VerticalAlignment','middle');
end


subplot(222);
ternary('An','Opx','Cpx');
ternplot(P(mafind),Opx(mafind),Cpx(mafind),'*');
for i = 1:nmaf  
    plag = mgons{i,1}(:,1);
    opx = mgons{i,1}(:,2);
    cpx = mgons{i,1}(:,3);
    
    ternplot(plag,opx,cpx,'-k');
    [x,y] = tern2xy(plag,opx,cpx);
    
    t = text(mean(x),mean(y),mgons{i,2});
    set(t,'HorizontalAlignment','center','VerticalAlignment','middle');
end


subplot(224);
ternary('Ol','Opx','Cpx');
ternplot(Ol(mind),Opx(mind),Cpx(mind),'*');
for i = 1:num  
    ol = umgons{i,1}(:,1);
    opx = umgons{i,1}(:,2);
    cpx = umgons{i,1}(:,3);
    
    ternplot(ol,opx,cpx,'-k');
    [x,y] = tern2xy(ol,opx,cpx);
    
    t = text(mean(x),mean(y),umgons{i,2});
    set(t,'HorizontalAlignment','center','VerticalAlignment','middle');
end

return