close all; clear; clc;

addpath('C:\Users\stefa\Desktop\Cremona\A2S1_MA\Homeworks\Homeworks\HW4');
% model = mphload('SprucePlate')
% mphcd(model);

load data.mat

if not(isfolder("plots"))
    mkdir("plots")
end
%% PLOTTING AND DATA EXTRACTION
close all;
% f = model.sol('sol2').getPVals();
numEval = 35;

plotting = 1; % plottare =1, non plottare=0
saving = 0; %salvare =1, non salvare=0

for ii=1:numEval

%     solnum = ii+44;
%     Zin(ii).imp = mphevalpoint(model, 'solid.F_Pz/solid.u_tZ', ...
%         'selection', 5, ...
%         'dataset','dset2', ...
%         'outersolnum', ii)';
%     Zin(ii).imp_dB = mag2db(abs(Zin(ii).imp));
%     sol = "sol"+num2str(solnum);
%     Zin(ii).names = model.sol(sol).getParamNames();
%     Zin(ii).vals = model.sol(sol).getParamVals() ;
%     mytitle(ii) = model.sol(sol).label();
    if (plotting)
        figure('Renderer', 'painters', 'Position', [100 100 1000 500]);
        plot(f, abs(Zin(ii).imp), LineWidth=2);
%         mytitle = model.sol(sol).label();
        title('Input Impedance',mytitle(ii), fontSize=20);
        xlabel('freq [Hz]', FontSize=16), ylabel('|Z_{in}| [Ns/m]', FontSize=16)
        grid minor
        if (saving)
            filename=strcat("inImp_",Zin(ii).names(1).toCharArray', num2str(Zin(ii).vals(1)), Zin(ii).names(2).toCharArray', num2str(Zin(ii).vals(2)));
            saveas(gcf, [".\plots\"+filename+".png"]);
        end
    end
end

% saveas(gcf, [".\plots\impedanceGrid.png"]);
% close all;
%% Serve a trovare la locazione e il valori d'impedenza della frequenza più vicina a quella desiderata
f_strings = [349.23, 440, 523.25, 659.25, 783.99];

numEval = 35;
for freq = f_strings
    Zmin(f_strings==freq).val = freq;
    for ii=1:numEval
%         minima = islocalmin(Zin(ii).imp_dB);
%         f_minima = f(minima);
%         err = abs(freq-f_minima);
        err_f = abs(f-freq);
        [min_err, min_loc] = min(err_f);
%         f_min = f_minima(min_loc);
        f_min = f(min_loc);
        imp_min = abs(Zin(ii).imp(f==f_min)); %imp_dB
        Zmin(f_strings==freq).f(ii) = f_min;
        Zmin(f_strings==freq).imp(ii) = imp_min;
        Zmin(f_strings==freq).xy(ii,1) = Zin(ii).vals(1);
        Zmin(f_strings==freq).xy(ii,2) = Zin(ii).vals(2);
        Zmin(f_strings==freq).minloc = min_loc;
    end
    [Zmin_valChosen, Zmin_locChosen] = min(Zmin(f_strings==freq).imp);
    Zmin(f_strings==freq).mincoord = Zin(Zmin_locChosen).vals;
end


%% OPTIMA POSITIONS
% poi in realtà questo blocco di codice serve anche a trovare i minimi
% della tavola, ma visto che il bridge è stato scelto a occhio l'algoritmo
% di ricerca dei minimi non serve.

x_coord = 0.05 : .1 : .95;
y_coord = 0.05 : .1 : 1.35;
sz = [5,7];
min_ind = 20;
close all

for kk=1:length(f_strings)
    Z_tmp = zeros(5, 7);
    ind = 0;
    % distribuisco le impedenze in una matrice
    for xi = 1:length(x_coord)/2
        for yi = 1:length(y_coord)/2
            ind = ind+1;
            Z_tmp(xi,yi) = Zmin(kk).imp(ind);
        end
    end
    % specchio e concateno per ottenere l'intera tavola
    Z_tmp = Z_tmp';
    Z_h = [Z_tmp; flipud(Z_tmp)];
    Z = [Z_h, fliplr(Z_h)];
    stringImp(kk).Z = Z;

%     [six_min, six_loc] = sort(Z_h(:));
    [six_min, six_loc] = sort(Z(:));
    six_min = six_min(1:min_ind);
    six_loc = six_loc(1:min_ind);
    opt_pos(kk).name = f_strings(kk);
    opt_pos(kk).imp = six_min;
    opt_pos(kk).ind = six_loc;
    opt_pos(kk).Zh = Z_h; 
    
    figure('Renderer', 'painters', 'Position', [800 0 500 600]);
    s = ones(length(x_coord)*length(y_coord), 1)*300;
    x_coord_long = repelem(x_coord, 14);
    y_coord_long = repmat(y_coord, 1, 10);
    scatter3(x_coord_long, y_coord_long, Z(:), s, db(Z(:)), 'filled');
    subtitle = ['f_0 = ',num2str(f_strings(kk)),' Hz'];
    title('Scatter plot of impedance', subtitle, FontSize=16)
    view(0,90)
    colorbar
    saveas(gcf, ["./plots/Scatter"+f_strings(kk)+".png"])
end
 
% BRIDGE POSITION

% x_coord_halflong = repmat(x_coord, 1, 7);
% sz = [7, 10];
sz = [14, 10];

% close all
% figure('Renderer', 'painters', 'Position', [10 10 600 700]);
% col = ["red" "green" "blue" "cyan" "magenta"];
% for kk = 1:length(f_strings)
%     [bridge_y_coord, bridge_x_coord] = ind2sub(sz, opt_pos(kk).ind);
% %     y_coord_flip = flip(y_coord);
%     x_coord_flip = flip(x_coord);
%     siz = 50000./db2mag(opt_pos(kk).imp(:));
%     scatter3(x_coord(bridge_x_coord), y_coord(bridge_y_coord), opt_pos(kk).imp(:), siz, col(kk), 'filled');
%     hold on
% end

xlim([0, 1]); ylim([0, 1.4])
view(0,90)


gs = [1,2,3,4];
gt = [2,3,4,5];
G = graph(gs, gt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET MANUALLY
% qui vengono impostate le coordinate del ponte
% xg = [0.15 0.25 0.65 0.75 0.85];
xg = [.65 .45 .35 .35 .35];
% yg = [0.85 0.95 1.15 1.15 1.15];
yg = [.55 .85 1.05 1.15 1.25];
% sg = [783.99, 440, 659.25, 349.23, 523.25];
% sg = [659.25, 440, 783.99, 523.25, 349.23];
sg = f_strings;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot vari del ponte 

zg = 7000*ones(1, length(xg));
Z_overlap = zeros(14,10);
for ii = 1:length(f_strings)
    Z_overlap = Z_overlap + stringImp(ii).Z;
end
figure('Renderer', 'painters', 'Position', [0 0 700 1000]);
s = ones(length(x_coord)*length(y_coord), 1)*600;
scatter3(x_coord_long, y_coord_long, Z_overlap(:), s, db(Z_overlap(:)), 'filled');
title('Scatter plot of total impedance', FontSize=20)
view(0,90)
colorbar

% figure('Renderer', 'painters', 'Position', [10 10 500 700]);
hold on
plot(G, 'XData',xg, 'YData',yg, 'ZData',zg, ...
    LineWidth=3, EdgeColor='k', NodeColor='none', EdgeAlpha=1, ...
    NodeLabelMode='auto', NodeLabel={})
xlim([0, 1]); ylim([0, 1.4])
saveas(gcf, "./plots/TotalImpedanceScatter.png")
% legend(["349", "440", "523", "659", "783", "bridge"], Location="bestoutside")
% 
figure('Renderer', 'painters', 'Position', [620 10 500 700]);
h = plot(G, 'XData',xg, 'YData',yg, ...
    LineWidth=5, EdgeColor='k', NodeColor='k',EdgeAlpha=1, ...
    NodeLabelMode='auto', NodeLabel={});
xlim([0, 1]); ylim([0, 1.4])
text(h.XData-.04, h.YData-0.05 ,h.NodeLabel, ...
    'VerticalAlignment','Bottom',...
    'HorizontalAlignment', 'left',...
    'FontSize', 7)
h.NodeLabel = {};
text(xg-0.1, yg, ['F4'; 'A4'; 'C5'; 'E5'; 'G5'], FontSize=20);
grid on
title('Final Shape of the Bridge', FontSize=20)
saveas(gcf, "./plots/BridgeFinalShape.png");



% % Z_overlap = Z_overlap/5;
% figure('Renderer', 'painters', 'Position', [0 0 300 400]);
% s = ones(length(x_coord)*length(y_coord), 1)*300;
% scatter3(x_coord_long, y_coord_long, Z_overlap(:), s, Z_overlap(:), 'filled');
% title('Scatter plot of impedance')
% view(0,90)
% colorbar

%% DECAY OF STRING
% close all
% clear string
% 
% for ii=1:length(sg)
%     string(ii).freq = sg(ii);
%     ii_x = xg(ii);
% %     string(ii).length = ii_x;
%     if(ii_x>0.5)
%         ii_x = 1-ii_x;
%     end
% 
%     ii_y = yg(ii);
% %     string(ii).length = ii_y;
%     if(ii_y>0.7)
%         ii_y = 1.4-ii_y;
%     end
% 
%     string(ii).pos = [ii_x; ii_y];
% end
% [so, sind] = sort(sg); 
% 
% save = 0;
% 
% for ii=1:length(sg)
% % ii=4;
% % ii=5;
% 
% zpos = [Zin.vals];
% zpos_x = zpos(1,:)';
% zpos_y = zpos(2,:)';
% pos_x = string(ii).pos(1);
% pos_y = string(ii).pos(2);
% indx = find(abs(zpos_x-pos_x)<0.001);
% indy = find(abs(zpos_y-pos_y)<0.001);
% ind = intersect(indx, indy);
% 
% 
% % L = string(ii).length;
% f0 = string(ii).freq;
% w0 = f0*2*pi;
% mu = 10.8e-3;
% 
% % T = (f0*L/2/pi)^2*mu;
% % T = (f0*L*2)^2*mu;
% T = 890;
% c = sqrt(T/mu);
% string(ii).length = c/f0/2;
% k = (w0)/c;
% Zc_string = T/c;
% Z_endString = 1j*Zc_string*cot(k*L);
% string(ii).tens = T;
% 
% 
% 
% Y_plate_f0 = 1/Zin(ind).imp(Zmin(find(sg(ii)==f_strings)).minloc);
% string(ii).adm = Y_plate_f0;
% string(ii).imp = (1/Y_plate_f0);
% % chi = j*Y_plate_f0*Z_endString/pi;
% chi = 1i*Y_plate_f0*Zc_string/pi;       
% 
% a1 = @(eps) chi + eps + sqrt(eps.^2+chi^2);
% a2 = @(eps) chi + eps - sqrt(eps.^2+chi^2);
% eps = linspace(-0.2, 0.2, 1000);
% 
% 
% % epsy = 0.002;
% epsy = abs(sqrt(imag(chi)^2 - real(chi)^2))/4;
% % epsy=0;
% 
% [V,D]= eig([chi, chi; chi, chi+2*epsy]);
% force = V(:,2);
% alpha1 = w0*a1(epsy);
% alpha2 = w0*a2(epsy);
% beta1 = w0+alpha1;
% beta2 = w0+alpha2;
% % disp(D)
% 
% % string(ii).eig_freq1 = real(beta1)/2/pi;
% % string(ii).eig_freq2 = real(beta2)/2/pi;
% string(ii).eig_freq1 = f0*(1+real(a1(epsy)));
% string(ii).eig_freq2 = f0*(1+real(a2(epsy)));
% string(ii).eps = epsy;
% string(ii).tau(1) = 1/(imag(alpha1));
% string(ii).tau(2) = 1/(imag(alpha2));
% 
% dur = 10;
% Fs = 48000;
% time = linspace(0,dur,dur*Fs);
% force1 = force(1)*exp(1i*beta1*time);
% force2 = force(2)*exp(1i*beta2*time);
% 
% 
% % v_b = abs(Y_plate_f0)*(force1+force2)*exp(1i*angle(Y_plate_f0));
% % v_b = (Y_plate_f0)*(force1+force2);
% % 
% mi = sqrt(epsy^2+chi^2);
% % % F0 = exp(1i*w0*time);
% % 
% F0=100;
% % force1 = F0/2/mi.*((chi-epsy+mi)*exp(1i*a1(epsy)*time)+(epsy-chi+mi)*exp(1i*a2(epsy)*time));
% % force2 = F0/2/mi/chi.*((epsy+mi)*(chi-epsy+mi)*exp(1i*a1(epsy)*time)+(epsy-mi)*(epsy-chi+mi)*exp(1i*a2(epsy)*time));
% v_b = F0*2*pi*chi/mi/Zc_string .* exp(1i*(epsy+chi)*w0*time).*(mi*cos(mi*w0*time)+1i*chi*sin(mi*w0*time)).*exp(1i*w0*time);
% % v_b = Y_plate_f0*(force1+force2).*exp(1i*w0*time);
% 
% % figure()
% figure('Renderer', 'painters', 'Position', [(300*(ii-1))+5 10 300 400]);
% subplot 211
% plot(time, real(v_b));
% % plot(time, real(force1));
% grid minor
% xlabel('time'); ylabel('v_B [m/s]')
% title('Velocity of plate', ["f_0="+num2str(sg(ii))+" Hz"]);
% % legend
% 
% string(ii).signal = real(v_b)./max(real(v_b));
% string(ii).env = mag2db(envelope(real(v_b), 500, 'rms'));
% 
% % figure()
% subplot 212
% plot(time, string(ii).env, LineWidth=1.2);
% % plot(time, real(force2))
% grid minor
% xlim([0,1]); ylim([-100, -20])
% title('Decaying of the plate velocity', ["f_0="+num2str(sg(ii))+" Hz"])
% xlabel('time [s]'); ylabel('v_{B} [dB]')
% 
%     if(save)
%         filename = ["ExD_vel_"+num2str(sg(ii))];
%         saveas(gcf, [".\plots\"+filename+".png"]);
%     end
% end

% f0 = sg(4);
% figure()
% plot(f0*eps, (1+real(a1(eps))-real(chi))*f0, f0*eps, (1+real(a2(eps))-real(chi))*f0, LineWidth=2)
% grid minor
% hold on
% % yline(f0, 'k--')
% plot(f0*eps, f0*(eps./eps), 'k--')
% hold on
% plot(f0*eps, f0*(1+2*eps), 'k--')
% xlabel('f0\cdot eps'); ylabel('f0\cdot \Re{(a_{\pm})}')
% legend('a_+', 'a_-')
% xlim([min(eps)*f0, max(eps)*f0])




%%

% for ii = 1:length(f_strings)
%     audiowrite(["./audio/string"+num2str(ii)+".wav"], string(ii).signal, Fs)
% end

%% DECAY SECONDO TENTATIVO
close all
clear string
save=1;

for ii=1:length(sg)
    string(ii).freq = sg(ii);
    ii_x = xg(ii);
    if(ii_x>0.5)
        ii_x = 1-ii_x;
    end
    ii_y = yg(ii);
    if(ii_y>0.7)
        ii_y = 1.4-ii_y;
    end
    string(ii).pos = [ii_x; ii_y];
end
[so, sind] = sort(sg); 

save = 0;

for ii=1:length(sg)
% ii=2;
% ii=5;

zpos = [Zin.vals];
zpos_x = zpos(1,:)';
zpos_y = zpos(2,:)';
pos_x = string(ii).pos(1);
pos_y = string(ii).pos(2);
indx = find(abs(zpos_x-pos_x)<0.001);
indy = find(abs(zpos_y-pos_y)<0.001);
ind = intersect(indx, indy);

f0 = sg(ii);
w0 = 2*pi*f0;
T = 890;
mu = 10.8e-3;
Z0 = sqrt(T*mu);
string(ii).imp = Zin(ind).imp(Zmin(sg(ii)==f_strings).minloc);
string(ii).adm = 1/string(ii).imp;
Yb = string(ii).adm;

chi = 1i*Yb*Z0/pi;
% chi = chi/20;
epsy = abs(sqrt(imag(chi)^2 - real(chi)^2))/5;
% epsy=0.0000002;
% epsy = 0.002;

a_p = @(eps) chi+eps+sqrt(eps.^2+chi^2);
a_m = @(eps) chi+eps-sqrt(eps.^2+chi^2);

string(ii).eig_freq1 = f0*(1+real(a_p(epsy)));
string(ii).eig_freq2 = f0*(1+real(a_m(epsy)));
string(ii).eps = epsy;
string(ii).tau(1) = 1/(imag(w0*a_p(epsy)));
string(ii).tau(2) = 1/(imag(w0*a_m(epsy)));
% string(ii).tau = 1/(w0*imag(chi));

mi = sqrt(epsy^2+chi^2);

dur=10;
Fs=48000;
time=linspace(0, dur, Fs*dur);
% time = 0:0.01:dur;

Vb = 2*pi*chi/(mi*Z0)*exp(1i*(epsy+chi)*w0*time).*(mi*cos(mi*w0*time)+1i*chi*sin(mi*w0*time)).*exp(1i*w0*time);
% string(ii).env = mag2db(envelope(real(Vb), 1000, 'rms'));
% string(ii).env = mag2db(abs(Vb));
Vb = Vb./max(abs(Vb));
string(ii).env = db(Vb);


figure('Renderer', 'painters', 'Position', [(300*(ii-1))+5 10 300 400]);
subplot 211
plot(time, real(Vb), LineWidth=1.2);
xlabel('time [s]'); ylabel('V_B [m/s]');
title('time evolution')
grid minor
xlim([0,2])

subplot 212
plot(time, string(ii).env, LineWidth=1.2);
xlabel('time [s]'); ylabel('V_B [dB]');
title('envelope decay')
grid minor
xlim([0,2])
sgtitle(["string with f_0 = "+num2str(sg(ii))+" Hz"]);


filename = ["ExD_vel_"+num2str(sg(ii))];
saveas(gcf, [".\plots\"+filename+".png"]);

%     find(abs(string(ii).env-(max(string(ii).env)-60))<0.002, 1)
string(ii).T60 = time(find(abs(string(ii).env-(max(string(ii).env)-60))<0.002, 1));

end

%% ULTERIORE TENTATIVO

close all
clear string
save=1;

for ii=1:length(sg)
    string(ii).freq = sg(ii);
    ii_x = xg(ii);
    if(ii_x>0.5)
        ii_x = 1-ii_x;
    end
    ii_y = yg(ii);
    if(ii_y>0.7)
        ii_y = 1.4-ii_y;
    end
    string(ii).pos = [ii_x; ii_y];
end
[so, sind] = sort(sg); 

save = 0;
time = linspace(0, 20, 48000*20);
notes = ["F4";"A4";"C5";"E5";"G5"];

for ii=1:length(sg)
    zpos = [Zin.vals];
    zpos_x = zpos(1,:)';
    zpos_y = zpos(2,:)';
    pos_x = string(ii).pos(1);
    pos_y = string(ii).pos(2);
    indx = find(abs(zpos_x-pos_x)<0.001);
    indy = find(abs(zpos_y-pos_y)<0.001);
    ind = intersect(indx, indy);
    
    f0 = sg(ii);
    w0 = 2*pi*f0;
    T = 700;
    mu = 10.8e-3;
    Z0 = sqrt(T*mu);
    string(ii).imp = Zin(ind).imp(Zmin(sg(ii)==f_strings).minloc);
    string(ii).adm = 1/string(ii).imp;
    Yb = string(ii).adm;
    
    X = 1i*w0*Yb*Z0/pi;
    eps = linspace(-0.02, 0.02, 1000);
    epsi = eps*w0;
    beta1 = @(eps) X + eps + sqrt(eps.^2+X^2) - real(X);
    beta2 = @(eps) X + eps - sqrt(eps.^2+X^2) - real(X);

    chi = 1i*imag(X);
    supp1 = @(eps) chi + eps + sqrt(eps.^2+chi^2);
    supp2 = @(eps) chi  + eps - sqrt(eps.^2+chi^2);
    
%     close all
%     figure()
%     figure('Renderer', 'painters', 'Position', [(300*(ii-1))+5 310 500 400]);
% %     plot(epsi, real(beta1(epsi)), epsi, real(beta2(epsi)), LineWidth=1.5)
%     plot(epsi, imag(beta1(epsi)), epsi, imag(beta2(epsi)), LineWidth=1.5)
%     hold on
% %     plot(epsi, zeros(1, length(epsi)), 'k--', LineWidth=1.5)
%     plot(epsi, imag(supp1(epsi)), 'k--', LineWidth=1.5)
%     hold on
% %     plot(epsi, 2*epsi, 'k--', LineWidth=1.5)
%     plot(epsi, imag(supp2(epsi)), 'k--', LineWidth=1.5)
%     grid minor
%     title(["Detuning Effect on Time Decay- "+notes(ii)], ["f_0="+num2str(sg(ii))+"Hz"], FontSize=15)
%     xlabel('\epsilon', FontSize=15); ylabel('\Im(\alpha)', FontSize=15)
%     xlim([min(epsi), max(epsi)])
% 
%     saveas(gcf, ["./plots/DetuningTime_"+notes(ii)+".png"]);
    
    epsy = 0.002*w0;
%     epsy = abs(sqrt(imag(chi)^2 - real(chi)^2));
    
    mi = sqrt(epsy^2+X^2);
    R = abs(((X+mi)*exp(1i*(X+mi)*time) - (X-mi)*exp(1i*(X-mi)*time))/2/mi).^2;
    
    figure('Renderer', 'painters', 'Position', [(300*(ii-1))+5 10 600 400]);
    plot(time, pow2db(R), LineWidth=1.5);
    xlim([0,2])
    title(["Decay curve - "+notes(ii)], ["f_0="+num2str(sg(ii))+"Hz"], FontSize=15)
    xlabel('time [s]', FontSize=15); ylabel('R(t)', FontSize=15);
    grid minor
    saveas(gcf, ["./plots/DecayTime_"+notes(ii)+".png"])
    
    string(ii).env = db(R);
    string(ii).eig1 = f0+real(beta1(epsy))/2/pi;
    string(ii).eig2 = f0+real(beta2(epsy))/2/pi;
    string(ii).tau1 = 1/(imag(beta1(epsy)));
    string(ii).tau2 = 1/(imag(beta2(epsy)));
    string(ii).T60_1 = -log(1e-6)*string(ii).tau1;
    string(ii).T60_2 = -log(1e-6)*string(ii).tau2;
%     string(ii).tau1 = 1/(imag(X+epsy));
%     string(ii).tau2 = 1/(imag(X-epsy));

    string(ii).T60 = time(find(abs(string(ii).env-(max(string(ii).env)-60))<0.002, 1));

%     figure()
%     subplot 311
%     plot(time, abs(exp(1i*(beta1(epsy)+w0)*time)));
%     grid minor
%     subplot 312
%     plot(time, abs(exp(1i*(beta2(epsy)+w0)*time)));
%     grid minor
%     subplot 313
%     plot(time, abs(exp(1i*(beta2(epsy)+w0)*time)+exp(1i*(beta1(epsy)+w0)*time)));
%     grid minor
    
end



