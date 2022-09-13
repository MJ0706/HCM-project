function active_tension()
lamda = 1.5
%I = Invariant(lamda)
timevariantformulation()

end
function timevariantformulation()
%Constant values
lr = 1.85;
lo1 = [1.55, 1.55, 1.55];
CaoMax = 4.35;
Cao = 4.35;
B = 4.75;
to1 = [350, 400, 350]; %450];
ttrans1 = [385, 430,  420];
tau1 = [20, 35, 35];
BCL1 = [1000, 910, 1180];
c = 1 %[1, 0.95, 1.8, 2.9]
T = [620E03, 400E03, 99.75E03]; %88.35E03];

%variables 
lamda = 1.0;
lambda = 1.0;
% timepoint =linspace(0, BCL, BCL/2);
% sigma_active = zeros(1, length(timepoint));
% dlamda = lamda/to;
actmax = zeros(1, length(T));
for j = 1:length(T)
    Tmax = T(j);
    lo = lo1(j);
    to = to1(j);
    ttrans = ttrans1(j);
    tau = tau1(j);
    BCL = BCL1(j);
    
    timepoint =linspace(0, BCL, BCL/2);
    timenorm = timepoint/BCL;
    sigma_active = zeros(1, length(timepoint));
    dlamda = lamda/to;
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
plot(timenorm, sigma_active/1000,'LineWidth',1.5)
hold on
actmax(j) = max(sigma_active)/1000.0
end

grid on
% legend({'k = 0','k = 0.14','k = 0.22','k = 0.265'}, 'FontSize',12)%,'17','18','19')%,'change','pp') %,'0.35,0.05','0.35,0.09','0.35,0.5')
legend({'Control','Non-Obstructive', 'Obstructive'}, 'FontSize',16, 'Fontname','Times New Roman')
%title('Isometric tension plot','FontSize',18)
xlabel('Normalized time','FontSize',22, 'fontWeight','bold', 'Fontname','Times New Roman')
ylabel('Tension (kPa)','FontSize',22, 'fontWeight','bold', 'Fontname','Times New Roman')
hold off
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
% ax.XLim = [25 60];
%pbaspect([1 1 1])
%grid minor
%plot_func (list, sigma_total_xx)
saveas(gcf, 'activetension_p1_3_wodispersion_otg', 'png')
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