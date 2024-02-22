% Sebastian J. Schlecht, Thursday, 15. February 2024
clear; clc; close all;

rng(1);

%% Define FDN
N = 4;
% m = [13 19 23];
% m = [3 9 13];
m = [142 231 303 343];
g = 0.99; 1; 0.99; 1; 
A = randomOrthogonal(N) * diag(g.^m);
b = randn(N,1);
c = randn(1,N);
d = randn(1,1);

%% Modal decomposition of FDN
[residues, poles, direct, isConjugatePolePair] = dss2pr(m,A,b,c,d);

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
    V = V * S;

    denominator = W' * (V .* m.' .* pole.^(m'-1) );

    rv(:,itN) = V / denominator;
    lv(:,itN) = W;

    

    ok = 1;
end


%% Test residue match

res_compact = (lv' * b) .* (c * rv).';

max(abs([residues - res_compact])) % residues match


%% For plotting only
% check impulse response
irLen = 1000;
ir_impz = dss2impz(irLen,m,A,b,c,d);
ir_pr = pr2impz(residues, poles, direct, isConjugatePolePair,irLen);

% expand eigenvalues across delay lines by multiplying the eigenvalue
for it = 1:N
    RV{it} = rv(it,:) .* poles.'.^((0:m(it)-1)).';
    LV{it} = lv(it,:) .* conj(poles.').^((m(it)-1):-1:0).';
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

figure; hold on; grid on; % plot Eigenvectors
% find low frequency modes for nicer visualization
[~,slowInd] = mink(angle(poles), 3);
for it = 1:N
    plot(real(RV{it}(:,slowInd))+it);
    set(gca,'ColorOrderIndex',1);
end
legend('Mode 1','Mode 2','Mode 3')
ylabel('Delay line number');
xlabel('Delay line samples')




