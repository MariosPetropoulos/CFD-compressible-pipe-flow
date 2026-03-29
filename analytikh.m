clear; clc;

%% Δεδομένα ροής
gamma = 1.4;
Rgas  = 287;      % αν το χρειαστείς μετά για T,ρ,κτλ
Tt_in = 288.0;
pt_in = 1.0332e5;

Mis_out = 0.59;   % ίδιο με αυτό που χρησιμοποίησες για p_out στη C++

%% Φόρτωμα γεωμετρίας από C++
% geometry.txt: 2 στήλες -> x, S(x) (area-like, π.χ. r^2)
geom = load("geometry_d.txt");
data=load('mach.txt');
data2=load("mach_b.txt");
data3=load("mach_c.txt");
data4=load("mach_d.txt");
data5 = load("mach_wave.txt");

Machh=data(:,2);
Machb = data2(:,2);
Machc = data3(:,2);
Machd = data4(:,2);
Mach_wave = data5(:,2);

x = geom(:,1);
A = geom(:,2);    % "S" από τον solver = διατομή (ή r^2)

N = numel(x);

%% Ορισμός της F(M) της σχέσης Area–Mach
F = @(M) (1./M) .* ...
    ((2./(gamma+1)) .* (1 + 0.5*(gamma-1).*M.^2)).^((gamma+1)/(2*(gamma-1)));

%% Βρες A* από το Mis_out
% Θεωρούμε ότι η φυσική έξοδος είναι στο x ≈ 1.0
[~, idx_exit] = min(abs(x - 1.0));
Ae = A(idx_exit);

F_out  = F(Mis_out);
A_star = Ae / F_out;   % Ae/A* = F(Mis_out)  =>  A* = Ae / F(Mis_out)

fprintf('Exit area Ae = %g, A* = %g\n', Ae, A_star);

%% Υπολογισμός αναλυτικού Mach(x) (ισεντροπική λύση, subsonic branch)

M_analytic = zeros(N,1);

for i = 1:N
    Aratio = A(i) / A_star;   % A(x)/A*
    
    % Λύνουμε F(M) = Aratio για M in (0,1) (subsonic λύση)
    % Με απλή bisection: F(M) μειώνεται από +∞ (M→0) σε 1 (M→1-).
    f = @(M) F(M) - Aratio;
    
    M_lo = 1e-4;
    M_hi = 0.999;
    
    for it = 1:50
        M_mid = 0.5*(M_lo + M_hi);
        if f(M_mid) > 0
            % F(M_mid) > Aratio -> ρίζα πιο δεξιά (σε μεγαλύτερο M)
            M_lo = M_mid;
        else
            % F(M_mid) < Aratio -> ρίζα πιο αριστερά (σε μικρότερο M)
            M_hi = M_mid;
        end
    end
    
    M_analytic(i) = 0.5*(M_lo + M_hi);
end

%% (Προαιρετικό) Υπολογισμός T(x), p(x), ρ(x), u(x) από M(x)

T = Tt_in ./ (1 + 0.5*(gamma-1).*M_analytic.^2);
p = pt_in .* (T./Tt_in).^(gamma/(gamma-1));
rho = p ./ (Rgas*T);
a = sqrt(gamma*Rgas*T);
u = M_analytic .* a;

%% (Προαιρετικό) Plot για σύγκριση με CFD (αν έχεις ήδη M_num)
% figure
% hold on
% 
% % Αναλυτική λύση
% plot(x, M_analytic, 'k-', 'LineWidth', 2)
% 
% % FVS
% plot(x, Machh, '--', ...
%     'Color', [0.7 0.0 0.0], ...   % κόκκινο
%     'LineWidth', 1)
% 
% plot(x, Machb, '-', ...
%     'Color', [0.7 0.0 0.0], ...
%     'LineWidth', 1)
% 
% % Roe
% plot(x, Machc, '--', ...
%     'Color', [0.0 0.1 0.5], ...   % σκούρο μπλε
%     'LineWidth', 1)
% 
% plot(x, Machd, '-', ...
%     'Color', [0.0 0.1 0.5], ...
%     'LineWidth', 1)
% 
% legend('M_{analytic}', ...
%        'FVS 1ης τάξης', 'FVS 2ης τάξης', ...
%        'Roe 1ης τάξης', 'Roe 2ης τάξης', ...
%        'Location', 'NorthEast')
% 
% xlabel('x')
% ylabel('Mach number')
% title('ΚΑΤΑΝΟΜΗ MACH')
% grid on
% box on
% set(gca,'FontSize',10)

% figure
% % plot(x, M_analytic, 'k')
% % hold on
% plot(x,Machd,'r')

figure

subplot(2,1,1)
plot(x, Machd, 'b', 'LineWidth', 2)

xlabel('x')
ylabel('M')
grid on
subplot(2,1,2)
plot(x, Mach_wave, 'r', 'LineWidth', 2)
xlabel('x')
ylabel('M')
sgtitle('Ηχητικός Λαιμός')

grid on