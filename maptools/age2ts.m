function varargout = age2ts(varargin)
% updated geologic timescale 7 Dec 2022 to International Commission on
% Stratigraphy (www.stratigraphy.org) International Chronostratigraphic
% Chart v.2022/10

if nargin > 1
    error('There may be at most one arguement to ts2age.');
end

%fid = fopen('geologic_time_scale.csv');
%ts = textscan(fid,'%s %d %d','Delimiter',',');
%fclose(fid);

ts = timescale;

if nargin == 0
    varargout = {ts};
    return
else
    for i = 1:length(ts(:,1))
        ts{i,1} = lower(ts{i,1});
    end
    agename = varargin{1};
    if length(agename) < 1
        ages = [NaN NaN NaN];
    end
end

if ischar(agename)
    ages = findage(agename,ts);
elseif iscell(agename)
    for i = 1:length(agename)
        ages = findage(agename{i},ts);
    end
end

varargout = {ages};

return

function age = findage(agename,ts)

i = strmatch(lower(agename),ts(:,1));
if isempty(i)
    age = [NaN NaN NaN];
else
    age = [ts{i,2} 0.5*(ts{i,2}+ts{i,3}) ts{i,3}];
end

return



function ts = timescale

ts = {'Meghalayan',0,0.0042;
'Northgrippian',0.0042,0.008276;
'Greenlandian',0.008276,0.0117;
'Cenozoic',0,66;
'Quaternary',0,2.58;
'Holocene',0,0.0117;
'Pleistocene',0.0117,2.58;
'Pleistocene Stage 4',0.0117,0.129;
'Chibanian',0.129,0.774;
'Calabrian',0.774,1.8;
'Gelasian',1.8,2.58;
'Neogene',2.58,23;
'Pliocene',2.58,5.3;
'Piacenzian',2.58,3.6;
'Zanclean',3.6,5.333;
'Miocene',5.333,23.03;
'Messinian',5.333,7.246;
'Tortonian',7.246,11.63;
'Serravallian',11.63,13.82;
'Langhian',13.82,15.97;
'Burdigalian',15.97,20.44;
'Aquitanian',20.44,23.03;
'Paleogene',23.03,66;
'Oligocene',23.03,33.9;
'Chattian',23.03,27.82;
'Rupelian',27.82,33.9;
'Eocene',33.9,56;
'Priabonian',33.9,37.71;
'Bartonian',37.71,41.2;
'Lutetian',41.2,47.8;
'Ypresian',47.8,56;
'Paleocene',56,66;
'Thanetian',56,59.2;
'Selandian',59.2,61.6;
'Danian',61.6,66;
'Mesozoic',66,252;
'Cretaceous',66,145;
'Late Cretaceous',66,100.5;
'Maastrichtian',66,72.1;
'Campanian',72.1,83.6;
'Santonian',83.6,86.3;
'Coniacian',86.3,89.8;
'Turonian',89.8,93.9;
'Cenomanian',93.9,100.5;
'Early Cretaceous',100.5,145;
'Albian',100.5,113;
'Aptian',113,121.4;
'Barremian',121.4,129.4;
'Hauterivian',129.4,132.6;
'Valanginian',132.6,139.8;
'Berriasian',139.8,145;
'Jurassic',145,201.3;
'Late Jurassic',145,161.5;
'Tithonian',145,149.2;
'Kimmeridgian',149.2,154.8;
'Oxfordian',154.8,161.5;
'Middle Jurassic',161.5,174.7;
'Callovian',161.5,165.3;
'Bathonian',165.3,168.2;
'Bajocian',168.2,170.9;
'Aalenian',170.9,174.7;
'Early Jurassic',174.7,201.4;
'Toarcian',174.7,184.2;
'Pliensbachian',184.2,192.9;
'Sinemurian',192.9,199.5;
'Hettangian',199.5,201.4;
'Triassic',201.4,251.902;
'Late Triassic',201.4,237;
'Rhaetian',201.4,208.5;
'Norian',208.5,227;
'Carnian',227,237;
'Middle Triassic',237,247.2;
'Ladinian',237,242;
'Anisian',242,247.2;
'Early Triassic',247.2,251.902;
'Olenekian',247.2,251.2;
'Induan',251.2,251.902;
'Paleozoic',251.902,538.8;
'Permian',251.902,298.9;
'Lopingian',251.902,259.51;
'Changhsingian',251.902,254.14;
'Wuchiapingian',254.14,259.51;
'Guadalupian',259.51,273.01;
'Capitanian',259.51,264.28;
'Wordian',264.28,266.9;
'Roadian',266.9,273.01;
'Cisuralian',273.01,298.9;
'Kungurian',273.01,283.5;
'Artinskian',283.5,290.1;
'Sakmarian',290.1,293.52;
'Asselian',293.52,298.9;
'Carboniferous',298.9,358.9;
'Pennsylvanian',298.9,323.2;
'Late Pennsylvanian',298.9,307;
'Gzhelian',298.9,303.7;
'Kasimovian',303.7,307;
'Middle Pennsylvanian',307,315.2;
'Moscovian',307,315.2;
'Early Pennsylvanian',315.2,323.2;
'Bashkirian',315.2,323.2;
'Mississippian',323.2,358.9;
'Late Mississippian',323.2,330.9;
'Serpukhovian',323.2,330.9;
'Middle Mississippian',330.9,346.7;
'Visean',330.9,346.7;
'Early Mississippian',346.7,358.9;
'Tournaisian',346.7,358.9;
'Devonian',358.9,419.2;
'Late Devonian',358.9,382.7;
'Famennian',358.9,372.2;
'Frasnian',372.2,382.7;
'Middle Devonian',382.7,393.3;
'Givetian',382.7,387.7;
'Eifelian',387.7,393.3;
'Early Devonian',393.3,419.2;
'Emsian',393.3,407.6;
'Pragian',407.6,410.8;
'Lochkovian',410.8,419.2;
'Silurian',419.2,443.8;
'Pridoli',419.2,423;
'Ludlow',423,427.4;
'Ludfordian',423,425.6;
'Gorstian',425.6,427.4;
'Wenlock',427.4,433.4;
'Homerian',427.4,430.5;
'Sheinwoodian',430.5,433.4;
'Llandovery',433.4,443.8;
'Telychian',433.4,438.5;
'Aeronian',438.5,440.8;
'Rhuddanian',440.8,443.8;
'Ordovician',443.8,485.4;
'Late Ordivician',443.8,458.4;
'Hirnantian',443.8,445.2;
'Katian',445.2,453;
'Sandbian',453,458.4;
'Middle Ordivician',458.4,470;
'Darriwilian',458.4,467.3;
'Dapingian',467.3,470;
'Early Ordivician',470,485.4;
'Floian',470,477.7;
'Tremadocian',477.7,485.4;
'Cambrian',485.4,538.8;
'Furongian',485.4,497;
'Cambrian Stage 10',485.4,489.5;
'Jiangshanian',489.5,494;
'Paibian',494,497;
'Miaolingian',497,509;
'Guzhangian',497,500.5;
'Drumian',500.5,504.5;
'Wulian',504.5,509;
'Cambrian Series 2',509,521;
'Cambrian Stage 4',509,514;
'Cambrian Stage 3',514,521;
'Terreneuvian',521,538.8;
'Cambrian Stage 2',521,529;
'Fortunian',529,538.8;
'Precambrian',538.8,4000;
'Proterozoic',538.8,2500;
'Neoproterozoic',538.8,1000;
'Ediacaran',538.8,635;
'Cryogenian',635,720;
'Tonian',720,1000;
'Mesoproterozoic',1000,1600;
'Stenian',1000,1200;
'Ectasian',1200,1400;
'Calymmian',1400,1600;
'Paleoproterozoic',1600,2500;
'Statherian',1600,1800;
'Orosirian',1800,2050;
'Rhyacian',2050,2300;
'Siderian',2300,2500;
'Archean',2500,4000;
'Neoarchean',2500,2800;
'Mesoarchean',2800,3200;
'Paleoarchean',3200,3600;
'Eoarchean',3600,4000;
'Hadean',4000,4650};

return
