function B = controle_recalage_composites(fichier_table)
%% =========================================================================
%  CONTROLE_RECALAGE_COMPOSITES (17 nov., demandée par l'utilisateur)
%  Pour chaque session de la table : l'image composite utilisée
%  (compositecalcium_V2.tif, sinon compositecalcium.tif ; pages rouge, vert,
%  bleu, GCaMP) est-elle bien superposée au GCaMP et aux masques suite2p ?
%  Rien n'est modifié.
%
%  Trois mesures :
%   1. GCaMP de la composite contre meanImg de suite2p (même image, mesure
%      fiable) : décalage par corrélation de phase. Doit être [0 0], sinon
%      TOUS les canaux sont décalés par rapport aux masques.
%   2. Vert contre GCaMP (le vert contient la fuite du GCaMP) : corrélation
%      de phase, avec la netteté du pic.
%   3. Chaque canal contre les masques : décalage entier (±6 px) qui
%      maximise le contraste « dans les masques / hors masques », et le gain
%      de contraste (comme ajuster_canaux_masques).
%  Doute si : 1 non nul ; 2 >= 2 px avec pic net (> 0,05) ; 3 >= 2 px et
%  gain > 3 %.
%  =========================================================================

T = readtable(fichier_table);
n = height(T);
B = table(strcat(T.id_animal, '_', T.id_session, '_', T.id_champ), repmat({''}, n, 1), ...
    nan(n, 2), nan(n, 2), nan(n, 1), nan(n, 6), nan(n, 3), repmat({''}, n, 1), ...
    'VariableNames', {'session', 'composite', 'gcamp_vs_meanImg_px', 'vert_vs_gcamp_px', 'vert_pic', ...
    'masques_RVB_px', 'masques_gain_pct', 'doute'});
for i = 1:n
    d = T.dossier_suite2p{i};
    f = fullfile(d, 'compositecalcium_V2.tif');
    if ~isfile(f), f = fullfile(d, 'compositecalcium.tif'); end
    if ~isfile(f), B.doute{i} = 'pas de composite'; continue; end
    [~, B.composite{i}] = fileparts(f);
    I = cell(1, 4);
    for k = 1:4, I{k} = double(imread(f, k)); end
    F = load(fullfile(d, 'Fall.mat'), 'ops', 'stat', 'iscell');
    txt = {};
    % 1
    mi = double(F.ops.meanImg);
    if isequal(size(mi), size(I{4}))
        [s, ~] = decalage_phase(I{4}, mi);
        B.gcamp_vs_meanImg_px(i, :) = s;
        if any(s ~= 0), txt{end+1} = sprintf('GCaMP/meanImg [%d %d]', s); end %#ok<AGROW>
    else
        txt{end+1} = sprintf('tailles différentes (composite %dx%d, meanImg %dx%d)', size(I{4}), size(mi)); %#ok<AGROW>
    end
    % 2
    [s, pic] = decalage_phase(I{2}, I{4});
    B.vert_vs_gcamp_px(i, :) = s; B.vert_pic(i) = pic;
    if norm(s) >= 2 && pic > 0.05, txt{end+1} = sprintf('vert/GCaMP [%d %d]', s); end %#ok<AGROW>
    % 3
    masque = false(size(I{4}));
    for m = find(F.iscell(:, 1) == 1)'
        if iscell(F.stat), st = F.stat{m}; else, st = F.stat(m); end
        x = double(st.xpix(:)) + 1; y = double(st.ypix(:)) + 1;
        ok = x >= 1 & x <= size(masque, 2) & y >= 1 & y <= size(masque, 1);
        masque(sub2ind(size(masque), y(ok), x(ok))) = true;
    end
    noms = {'rouge', 'vert', 'bleu'};
    for c = 1:3
        [dec, gain] = meilleur_decalage(I{c}, masque);
        B.masques_RVB_px(i, 2 * c - [1 0]) = dec; B.masques_gain_pct(i, c) = 100 * gain;
        if norm(dec) >= 2 && gain > 0.03
            txt{end+1} = sprintf('%s/masques [%d %d] +%.0f%%', noms{c}, dec, 100 * gain); %#ok<AGROW>
        end
    end
    B.doute{i} = strjoin(txt, ' ; ');
    fprintf('%2d/%d %-28s %s\n', i, n, B.session{i}, B.doute{i});
end
end


function [dec, gain] = meilleur_decalage(I, masque)
contraste = @(J) mean(J(masque)) / mean(J(~masque));
c0 = contraste(I); best = c0; dec = [0 0];
for dx = -6:6
    for dy = -6:6
        v = contraste(decaler(I, dx, dy));
        if v > best, best = v; dec = [dx dy]; end
    end
end
gain = best / c0 - 1;
end


function J = decaler(I, dx, dy)
J = zeros(size(I));
[h, w] = size(I);
xs = max(1, 1 - dx):min(w, w - dx); ys = max(1, 1 - dy):min(h, h - dy);
J(ys + dy, xs + dx) = I(ys, xs);
end


function [s, force] = decalage_phase(a, b)
a = a - mean(a(:)); b = b - mean(b(:));
X = fft2(a) .* conj(fft2(b)); X = X ./ max(abs(X), eps);
C = real(ifft2(X));
[force, idx] = max(C(:));
[dy, dx] = ind2sub(size(C), idx);
n = size(C);
s = [mod(dx - 1 + floor(n(2) / 2), n(2)) - floor(n(2) / 2), mod(dy - 1 + floor(n(1) / 2), n(1)) - floor(n(1) / 2)];
end
