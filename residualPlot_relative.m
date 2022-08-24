clear all
clc

files = dir('*.mat');

maxF = 5000;        % Max frequency

% Initialise import indexing
powerCount = 1;

% Import data
for i = 1:length(files)
    load(files(i).name);
    if contains(convertCharsToStrings(files(i).name),'power')
        power(powerCount) = IntegratedFields;
        powerCount = powerCount + 1;
    end
end


% Initialise empty arrays
outPower4K = zeros(length(power(1).OutPower4K),length(power));
outPower77K = zeros(length(power(1).OutPower77K),length(power));
outPowerOVC = zeros(length(power(1).OutPowerOVC),length(power));
dispNorm4K = zeros(length(power(1).DisplacementNorm4K),length(power));
dispNorm77K = zeros(length(power(1).DisplacementNorm77K),length(power));
dispNormOVC = zeros(length(power(1).DisplacementNormOVC),length(power));

% Populate empty arrays
for i = 1:length(power)
    outPower4K(:,i) = power(i).OutPower4K;
    outPower77K(:,i) = power(i).OutPower77K;
    outPowerOVC(:,i) = power(i).OutPowerOVC;
    dispNorm4K(:,i) = power(i).DisplacementNorm4K;
    dispNorm77K(:,i) = power(i).DisplacementNorm77K;
    dispNormOVC(:,i) = (power(i).DisplacementNormOVC);
end

% End of data import

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kinetic Energy calculation

KE4K = dispNorm4K.^2;
KE77K = dispNorm77K.^2;
KEOVC = dispNormOVC.^2;

df = maxF / (length(dispNorm4K) + 1);
f = 0;

for i = 1:length(outPower4K)
    f = f + df;
    KE4K(i,:) = 0.5 * 4 * 2710 * ((2 * pi() * f)^2) * (KE4K(i,:));
    KE77K(i,:) = 0.5 * 4 * 22698 * ((2 * pi() * f)^2) * (KE77K(i,:));
    KEOVC(i,:) = 0.5 * 4 * 27900 * ((2 * pi() * f)^2) * (KEOVC(i,:));
end

% gamma = 0.000015;
% 
% outPower4K = 0.5 * 4 * gamma * (outPower4K.^2);
% outPower77K = 0.5 * 4 * gamma * (outPower77K.^2);
% outPowerOVC = 0.5 * 4 * gamma * (outPowerOVC.^2);

% Initialise empty arrays
outPower4K_resid = zeros(length(power(1).OutPower4K),length(power)-1);
outPower77K_resid = zeros(length(power(1).OutPower77K),length(power)-1);
outPowerOVC_resid = zeros(length(power(1).OutPowerOVC),length(power)-1);
KE4K_resid = zeros(length(power(1).DisplacementNorm4K),length(power)-1);
KE77K_resid = zeros(length(power(1).DisplacementNorm77K),length(power)-1);
KEOVC_resid = zeros(length(power(1).DisplacementNormOVC),length(power)-1);

% for ii = 1:length(power)-1
%     KE4K_resid(:,ii) = abs(KE4K(:,ii)-KE4K(:,length(power)));
%     KE77K_resid(:,ii) = abs(KE77K(:,ii)-KE77K(:,length(power)));
%     KEOVC_resid(:,ii) = abs(KEOVC(:,ii)-KEOVC(:,length(power)));
%     outPower4K_resid(:,ii) = abs(outPower4K(:,ii)-outPower4K(:,length(power)));
%     outPower77K_resid(:,ii) = abs(outPower77K(:,ii)-outPower77K(:,length(power)));
%     outPowerOVC_resid(:,ii) = abs(outPowerOVC(:,ii)-outPowerOVC(:,length(power)));
% end

for ii = 1:length(power)-1
    KE4K_resid(:,ii) = abs((KE4K(:,ii)-KE4K(:,length(power)))./KE4K(:,length(power)))*100;
    KE77K_resid(:,ii) = abs((KE77K(:,ii)-KE77K(:,length(power)))./KE4K(:,length(power)))*100;
    KEOVC_resid(:,ii) = abs((KEOVC(:,ii)-KEOVC(:,length(power)))./KE4K(:,length(power)))*100;
    outPower4K_resid(:,ii) = abs((outPower4K(:,ii)-outPower4K(:,length(power)))./KE4K(:,length(power)))*100;
    outPower77K_resid(:,ii) = abs((outPower77K(:,ii)-outPower77K(:,length(power)))./KE4K(:,length(power)))*100;
    outPowerOVC_resid(:,ii) = abs((outPowerOVC(:,ii)-outPowerOVC(:,length(power)))./KE4K(:,length(power)))*100;
end

%Options

dataTag = 'disp_resid';
axisFontSize = 18;
legendFontSize = 18;

dispLegend = {'Disp $ = 1e^{-03}$m', 'Disp $ = 1e^{-04}$m', 'Disp $ = 1e^{-05}$m', 'Disp $ = 1e^{-06}$m','FontSize', legendFontSize,'Interpreter','Latex'};

% Plotting and saving 

if not(isfolder('figures'))
    mkdir('figures')
end

f = figure('visible','off');

xPower = linspace(0,maxF,length(outPower4K));

semilogy(xPower,outPowerOVC_resid,'LineWidth',2.0)
xlabel('$Frequency \; (Hz)$','Interpreter','Latex','FontSize', axisFontSize)
ylabel('$P^0 (W)$','Interpreter','Latex','FontSize', axisFontSize)
ylim([10e-9 10e6])
grid on
legend(dispLegend{:}, 'location', 'nw');
curfilepath = append('./figures/PowerOVC_', dataTag, '.eps');
exportgraphics(gcf,curfilepath,'Resolution',300)

semilogy(xPower,KEOVC_resid,'LineWidth',2.0)
xlabel('$Frequency (Hz)$','Interpreter','Latex','FontSize', axisFontSize)
ylabel('Kinetic Energy $(J)$','Interpreter','Latex','FontSize', axisFontSize)
ylim([10e-12 10e10])
grid on
legend(dispLegend{:}, 'location', 'nw');
curfilepath = append('./figures/KEOVC_', dataTag, '.eps');
exportgraphics(gcf,curfilepath,'Resolution',300)

semilogy(xPower,outPower77K_resid,'LineWidth',2.0)
xlabel('$Frequency \; (Hz)$','Interpreter','Latex','FontSize', axisFontSize)
ylabel('$P^0 (W)$','Interpreter','Latex','FontSize', axisFontSize)
ylim([10e-10 10e10])
grid on
legend(dispLegend{:}, 'location', 'nw');
curfilepath = append('./figures/Power77K_', dataTag, '.eps');
exportgraphics(gcf,curfilepath,'Resolution',300)

semilogy(xPower,KE77K_resid,'LineWidth',2.0)
xlabel('$Frequency (Hz)$','Interpreter','Latex','FontSize', axisFontSize)
ylabel('Kinetic Energy $(J)$','Interpreter','Latex','FontSize', axisFontSize)
ylim([10e-10 10e8])
grid on
legend(dispLegend{:}, 'location', 'nw');
curfilepath = append('./figures/KE77K_', dataTag, '.eps');
exportgraphics(gcf,curfilepath,'Resolution',300)

semilogy(xPower,outPower4K_resid,'LineWidth',2.0)
xlabel('$Frequency \; (Hz)$','Interpreter','Latex','FontSize', axisFontSize)
ylabel('$P^0 (W)$','Interpreter','Latex','FontSize', axisFontSize)
ylim([10e-9 10e8])
grid on
legend(dispLegend{:}, 'location', 'nw');
curfilepath = append('./figures/Power4K_', dataTag, '.eps');
exportgraphics(gcf,curfilepath,'Resolution',300)

semilogy(xPower,KE4K_resid,'LineWidth',2.0)
xlabel('$Frequency (Hz)$','Interpreter','Latex','Interpreter','Latex','FontSize', axisFontSize)
ylabel('Kinetic Energy $(J)$','Interpreter','Latex','FontSize', axisFontSize)
ylim([10e-12 10e6])
grid on
legend(dispLegend{:}, 'location', 'nw');
curfilepath = append('./figures/KE4K_', dataTag, '.eps');
exportgraphics(gcf,curfilepath,'Resolution',300)
