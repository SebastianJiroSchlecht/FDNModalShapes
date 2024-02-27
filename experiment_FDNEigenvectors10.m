% Sebastian J. Schlecht, Thursday, 15. February 2024
clear; clc; close all;

rng(1);

%% Define FDN
N = 3;
m = [13 19 23];
g = 0.98;
A = randomOrthogonal(N) * diag(g.^m);
b = ones(N,1);
c = ones(1,N);
c = [0 0 1];
d = zeros(1,1);

%% Modal decomposition of FDN
[residues, poles, direct, isConjugatePolePair] = dss2pr(m,A,b,c,d);
% sort according to pole frequency
residues = sortby(residues, angle(poles));
isConjugatePolePair = sortby(isConjugatePolePair, angle(poles));
poles = sortby(poles, angle(poles));

%% get eigenvectors compactly
numPoles = numel(poles);
rv = zeros(N,numPoles); % right eigenvector
lv = zeros(N,numPoles); % left eigenvector
for itN = 1:numPoles
    pole = poles(itN);

    % compute the adjugate
    E = diag(pole.^(m));
    P = E - A;
    adjP = adjugate(P);
    
    % find rank 1 decomposition
    [V,S,W] = svds(adjP,1);

    denominator = W' * (V .* m.' .* pole.^(m'-1) );

    rv(:,itN) = V / sqrt(denominator);
    lv(:,itN) = W / conj(sqrt(denominator));

end


%% Test residue match
res_compact = (lv' * b ) .* (c * rv).';
max(abs([residues - res_compact])) % residues match


%% For plotting only
% check impulse response
irLen = 1000;
ir_impz = dss2impz(irLen,m,A,b,c,d);
ir_pr = pr2impz(res_compact, poles, direct, isConjugatePolePair,irLen);

% expand eigenvalues across delay lines by multiplying the eigenvalue
for it = 1:N
    RV{it} = rv(it,:) .* poles.'.^((0:m(it)-1)).';
    LV{it} = lv(it,:) .* conj(poles.').^((m(it)-1):-1:0).';
end

% get different excitation position by delaying b

for it = 1:1:23
    inputDelay = poles.'.^[0, 0, it];
    inputDelay2 = 0*inputDelay;
    ind = 26;
    inputDelay2(ind,:) = inputDelay(ind,:);
    res_compact = (inputDelay .* lv' * b ) .* (c * rv).';
    ir_pr_pos(:,it) = pr2impz(res_compact, poles, direct, isConjugatePolePair,irLen);
end


%% Plot
figure; hold on; grid on;
plot(ir_impz)
plot(ir_pr+1)
legend('Impulse response (time-domain)','IR pole residue')
xlabel('Time (samples)')
ylabel('Impulse response value');

figure; hold on; grid on;
plot(real(poles),imag(poles),'x');
legend('Poles')
xlabel('Real axis')
ylabel('Imaginary axis')

figure; hold on; grid on;
plot(angle(poles),abs(residues),'x');
plot(angle(poles),abs(res_compact),'sq');
legend('Residues (Modal)', 'Residues (Compact Eigenvector)');
xlabel('Pole angle')
ylabel('Residue Magnitude')

figure; hold on; grid on;
plot(angle(poles),angle(residues),'x');
plot(angle(poles),angle(res_compact),'sq');
legend('Residues (Modal)', 'Residues (Compact Eigenvector)');
xlabel('Pole angle')
ylabel('Residue Angle')

% plot slow Eigenvectors
figure; hold on; grid on; 
for it = 1:N
    plot(real(RV{it}(:,7))+it);
    set(gca,'ColorOrderIndex',1);
end
% legend('Mode 1','Mode 2','Mode 3')
ylabel('Delay line number');
xlabel('Delay line samples')

for it = 1:N
    figure; hold on
    plot(real(RV{it}(:,11)),'-','LineWidth',8,'Marker','.','MarkerSize',19,'Color','black');
    plot(real(LV{it}(:,11)),'--','LineWidth',8,'Marker','.','MarkerSize',19,'Color',[1 1 1]*0.5);
    set(gca,'ColorOrderIndex',1);
    axis off
    xlim([0 max(m)])
    set(gcf,'Position',[0 0 1000 150])
    saveas(gcf,sprintf('./results/mode%d.png',it))
end

% plot all Eigenvectors
figure; hold on;
RVs = vertcat(RV{:});
plotMatrix(real(RVs))
plot3([1,numPoles+1],[1 1]*0+1,[100,100],'-k','LineWidth',3);
plot3([1,numPoles+1],[1 1]*m(1)+1,[100,100],'-k','LineWidth',3);
plot3([1,numPoles+1],[1 1]*sum(m([1 2]))+1,[100,100],'-k','LineWidth',3);
% colorbar
xlabel('Eigenvalue index $i$');
ylabel('State space index');
xlim([1 numPoles+1])
ylim([1 numPoles*2])

set(gcf,'Units', 'inches', 'Position', [0 0 3.5 4.9]);
exportgraphics(gcf,'./results/EigenvectorFDN.pdf')


figure; hold on;
plotMatrix(real(rv))
% plot separating lines
plot3([1,numPoles+1],[1 1]*1,[100,100],'-k','LineWidth',3);
plot3([1,numPoles+1],[1 1]*2,[100,100],'-k','LineWidth',3);
plot3([1,numPoles+1],[1 1]*3,[100,100],'-k','LineWidth',3);
colorbar('horiz')
% xlabel('Eigenvalue index $i$');
ylabel('Delay index');
xlim([1 numPoles+1])
ylim([1 N+1])

set(gcf,'Units', 'inches', 'Position', [0 0 3.5 1.1]);
exportgraphics(gcf,'./results/EigenvectorFDN_compact.pdf')


% moving intap
figure; hold on; grid on;

fs = 8000;
% soundsc(ir_pr_pos(:),8000);
IR_pos = fft(ir_pr_pos);
ff = linspace(0,fs,1000);
pp = (1:1:23) / 23;

colororder([parula(44)])
waterfall(pp, ff, db(IR_pos));
xlabel('Relative Intap Position')
ylabel('Frequency (Hz)')
zlabel('Magnitude (dB)')
ylim([100 fs/2/2])
zlim([-20 20])
view([67 57])

set(gcf,'Units', 'inches', 'Position', [0 0 3.5 2.5]);
exportgraphics(gcf,'./results/varyingIntap.pdf')

