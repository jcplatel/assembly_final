function [ar, ag, ab, info] = ajuster_canaux_masques(ar, ag, ab, dossier)
%% =========================================================================
%  AJUSTER_CANAUX_MASQUES (17 nov., demandée par l'utilisateur)
%  Après image_alignment_multimodal : vérifie que chaque canal (rouge, vert,
%  bleu) tombe bien sur les cellules. Pour chaque canal, on cherche le
%  décalage entier (jusqu'à ±6 px) qui maximise le contraste « moyenne dans
%  les masques des cellules (suite2p, iscell = 1) / moyenne hors masques ».
%  On ne l'applique que s'il fait au moins 2 px ET améliore le contraste de
%  plus de 3 % ; sinon le canal ne bouge pas. (Avec 1 %, le bleu de la
%  session GT 444119 220919 plane0 était décalé et les couleurs moins
%  justes : 89,6 % au lieu de 91,1 %.)
%  17 nov., contrôle des 57 sessions : corrigés à 3 % = 444113 220930 plane1
%  (vert, +11 %) et 444175 221121 plane1 (les 3 canaux, +4 à 8 %).
%  =========================================================================

d = [regexprep(char(dossier), '[\\/]+$', '') filesep];
[F.stat, F.iscell] = charger_masques(d);   % 18 nov. : masques corrigés s'ils existent
taille = size(ar);
masque = false(taille);
roi = find(F.iscell(:, 1) == 1);
for m = roi(:)'
    if iscell(F.stat), st = F.stat{m}; else, st = F.stat(m); end
    x = double(st.xpix(:)) + 1; y = double(st.ypix(:)) + 1;
    ok = x >= 1 & x <= taille(2) & y >= 1 & y <= taille(1);
    masque(sub2ind(taille, y(ok), x(ok))) = true;
end

noms = {'rouge', 'vert', 'bleu'};
I = {ar, ag, ab};
info = struct('canal', {noms}, 'decalage_px', zeros(3, 2), 'contraste_avant', zeros(1, 3), ...
    'contraste_meilleur', zeros(1, 3), 'applique', false(1, 3));
for c = 1:3
    contraste = @(J) mean(J(masque)) / mean(J(~masque));
    c0 = contraste(I{c}); meilleur = c0; dec = [0 0];
    for dx = -6:6
        for dy = -6:6
            v = contraste(decaler(I{c}, dx, dy));
            if v > meilleur, meilleur = v; dec = [dx dy]; end
        end
    end
    info.contraste_avant(c) = c0; info.contraste_meilleur(c) = meilleur; info.decalage_px(c, :) = dec;
    if norm(dec) >= 2 && meilleur / c0 - 1 > 0.03
        I{c} = decaler(I{c}, dec(1), dec(2));
        info.applique(c) = true;
    end
end
ar = I{1}; ag = I{2}; ab = I{3};
end


function J = decaler(I, dx, dy)
% translation entière, bords remplis de 0
J = zeros(size(I));
[h, w] = size(I);
xs = max(1, 1 - dx):min(w, w - dx); ys = max(1, 1 - dy):min(h, h - dy);
J(ys + dy, xs + dx) = I(ys, xs);
end
