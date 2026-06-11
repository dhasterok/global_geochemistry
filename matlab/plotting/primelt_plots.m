function primelt_plots()

c = 1;
for i = [0:0.01:1]
    for j = [0:0.01:1]
        if i + j > 1
            continue;
        end
        an(c) = i;
        ol(c) = j;
        c = c + 1;
    end
end
qz = 1 - an - ol;

% clinopyroxene
s_cpx = -0.074 + 0.1713./ol - 0.0135./ol.^2;

s_gtlhzhz = 1./(16.843 + 28.733*an - 14.183*exp(an));

domain = 2*ones(size(an));
domain(ol > 0.5 & qz < s_gtlhzhz) = 3;
domain(domain ~= 3 & qz > s_cpx) = 1;

% F1
ind = domain == 1;
Fproj = 6.2819 * an(ind).^2 - 14.7789 * an(ind).^3 + ...
    0.00825 * (1 ./ an(ind)).^2;

% F2
ind = domain == 2;
Fproj(ind) = ((-1.994 + 2.25 * qz(ind) + 0.041 ./ qz(ind)) + ...
    (-1.183 - 3.005 * qz(ind) + 13.774 * qz(ind).^2 - 12.615 * qz(ind).^3)) / 2 + ...
    exp(0.931 + 1.623.*qz(ind)) .* qz(ind).^0.245 .* ol(ind) + ...
    exp(0.769 - 7.514 * qz(ind)) .* qz(ind).^0.577 ./ ol(ind);
Fproj(ind & qz <= 0) = 0;

ind2 = ind & Fproj > 0 & qz < -0.1773 + 0.1404 ./ ol - 0.008434 ./ ol.^2;
Fproj(ind2) = -Fproj(ind2);

% F3
ind = domain == 3;
Fproj(ind) = -2.5345 + 5.329 * (qz(ind) + ol(ind)*0.348) + 0.3012 ./ ...
    (qz(ind) + 0.348 * ol(ind));

ind2 = ind & Fproj > 0 & qz < -0.0896 + 0.02002 ./ ol + 0.02989 ./ ol.^2;
Fproj(ind2) = -Fproj(ind2);

figure;
ternary('An','Ol','Qz');
hold on;

ternsurf(an,ol,qz,Fproj,0.02);
%ternplot(an,ol,qz,'o',{20,domain});

return