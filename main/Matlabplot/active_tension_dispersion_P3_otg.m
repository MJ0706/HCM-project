function active_tension()
lamda = 1.5
%I = Invariant(lamda)
timevariantformulation()

end
function timevariantformulation()
%Constant values
lr = 1.85;
lo = 1.55 ;
CaoMax = 4.35;
Cao = 4.35;
B = 4.75;
to = 350; %800 ; 
ttrans = 420; %410; 
tau = 35; %40; 
BCL = 1180; %910; 

c = [1.05, 1.125, 1.15, 1.3, 1.45, 1.75]
T = 95E03; %400E03;  

% c  = [2.0,2.6,2.7, 3.0,7.0,18.0 ]
% T = 400E03;


%variables 
lamda = 1.0;
lambda = 1.0;
timepoint =linspace(0, BCL, BCL/2);
sigma_active = zeros(1, length(timepoint));
dlamda = lamda/to;
actmax = zeros(1, length(c));
for j = 1:length(c)
    Tmax = c(j)*T;

for i = 1:length(timepoint)
    ta = timepoint(i);
%     if ta<to
%         lambda = lambda + dlamda;
%     else
%         lambda = lambda - dlamda;
%     end
    lambda = lamda;
    
    lso = ActiveLength(lambda, lr, lo);
    
    deno = sqrt(exp((B*lso)-1));
    ECa50 = CaoMax/deno;
    
    if ta<ttrans
        Ct = 0.5*(1-cos(pi*ta/to));
    else
        Ct = 0.5*(1-cos(pi*ttrans/to))*exp(-((ta-ttrans)/tau));
    end


    CaTerm = Cao^2 /(Cao^2 + ECa50^2);
    Pact = Tmax*CaTerm*Ct;
    sigma_active(i) = Pact;
    
end
    
figure(1)
hold on
plot(timepoint, sigma_active/1000,'LineWidth',1.5)
hold on
actmax(j) = max(sigma_active)/1000.0

end

grid on
legend({'\kappa = 0','\kappa = 0.07','\kappa = 0.1','\kappa = 0.14','\kappa=0.18','\kappa=0.22' }, 'FontSize',16, 'Fontname','Times New Roman')%,'17','18','19')%,'change','pp') %,'0.35,0.05','0.35,0.09','0.35,0.5')
%legend({'Control','Non-Obstructive', 'Obstructive'}, 'FontSize',16)
title('Obstructive HCM','FontSize',18, 'Fontname','Times New Roman', 'Fontweight' , 'normal')
xlabel('Time (ms)','FontSize',18, 'Fontname','Times New Roman')
ylabel('Tension (kPa)','FontSize',18, 'Fontname','Times New Roman')
hold off
ax = gca;
%ax.FontWeight = 'bold';
% ax.XLim = [25 60];
%pbaspect([1 1 1])
%grid minor
%plot_func (list, sigma_total_xx)
saveas(gcf, 'activetension_p3_dispersion_otg', 'png')
actmax
end

function lso = ActiveLength(lamda, lr, lo)
ls = lamda*lr;
if ls<=lo
    lso = 0.002;
else 
    lso = ls-lo;

end
end