function [xTop, yTop, zTop, alpha] = getLocalOpt(ls, curvature, Lgb, ks, kb)

%%  get x and alpha;
load('Data_Ls_more.mat');
ls = ls/Lgb;
xTop = interp1(Ls, xTop, ls); 
alpha = interp1(Ls, Alpha, ls);
yTopRef = interp1(Ls,yTop, ls);

%% get y

load('stretch_para_more.mat');

Pi_1 = ls; % bar l_s (NOTE: susLength is already normalized)
Pi_2 = Lgb^2 * ks / kb ; % bar k_s
relevantGroup = Pi_1.^2 ./ Pi_2;

delta =  interp1(Ls, Delta, ls, 'linear','extrap');
z_inf = interp1(Ls, z_inf, ls, 'linear','extrap');

yTop = delta * relevantGroup + z_inf;
if abs(yTop - yTopRef)/yTopRef > 0.02
    a = 'tdz';
end
xTop = xTop;
zTop = yTop;

%% get ytop
load('Y_para_more.mat');
delta =  interp1(Ls, DeltaY, ls, 'linear','extrap');
y_inf = interp1(Ls, Y_inf, ls, 'linear','extrap');
p = delta * relevantGroup + y_inf;
yTop = p * Lgb * curvature;
yTop = yTop;

end