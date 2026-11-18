function cols = couleurs_assemblees(n, fond)
%COULEURS_ASSEMBLEES  Couleurs des assemblées, partout dans le pipeline
%   (interface et figures PNG).
%
%   cols = couleurs_assemblees(n)            % fond noir (défaut) : couleurs vives
%   cols = couleurs_assemblees(n, 'blanc')   % fond blanc : les mêmes, plus foncées
%
%   Teintes choisies ENTRE les couleurs brainbow (rouge, jaune = RG, vert,
%   cyan = GB, bleu, magenta = RB) pour ne jamais les confondre : azur,
%   vert tilleul, violet, vert printemps, rose ; puis les mêmes en clair.
%   Au-delà de 10 assemblées : d'autres teintes intermédiaires, toujours les
%   mêmes d'un affichage à l'autre.

if nargin < 2 || isempty(fond), fond = 'noir'; end
vives = [0.18 0.71 1.00;    % azur           #2EB5FF
         0.66 0.94 0.24;    % vert tilleul   #A8F03D
         0.63 0.44 1.00;    % violet         #A170FF
         0.18 0.96 0.69;    % vert printemps #2EF5B0
         1.00 0.37 0.68;    % rose           #FF5EAD
         0.62 0.85 1.00;    % azur clair     #9ED9FF
         0.85 1.00 0.54;    % tilleul clair  #D9FF8A
         0.80 0.69 1.00;    % lavande        #CCB0FF
         0.61 1.00 0.85;    % menthe         #9CFFD9
         1.00 0.68 0.84];   % rose clair     #FFADD6
foncees = [0.04 0.52 0.84;  % azur
           0.37 0.68 0.00;  % vert tilleul
           0.48 0.25 0.88;  % violet
           0.00 0.66 0.48;  % vert printemps
           0.88 0.19 0.49;  % rose
           0.31 0.66 0.91;  % azur clair
           0.55 0.78 0.25;  % tilleul clair
           0.62 0.50 0.91;  % lavande
           0.24 0.79 0.63;  % menthe
           0.95 0.48 0.71]; % rose clair
switch fond
    case 'noir',  base = vives;
    case 'blanc', base = foncees;
    otherwise, error('couleurs_assemblees : fond ''noir'' ou ''blanc''.');
end

if n > size(base, 1)
    % teintes supplémentaires, écartées d'au moins 20° des teintes brainbow
    % (multiples de 60°)
    k = (1:n - size(base, 1))';
    teinte = mod(205 + 137.508 * k, 360);
    ecart = mod(teinte, 60);
    teinte(ecart < 20) = teinte(ecart < 20) - ecart(ecart < 20) + 20;
    teinte(ecart > 40) = teinte(ecart > 40) - ecart(ecart > 40) + 40;
    if strcmp(fond, 'noir'), v = 1; s = 0.65; else, v = 0.80; s = 0.85; end
    base = [base; hsv2rgb([teinte / 360, repmat(s, numel(k), 1), repmat(v, numel(k), 1)])];
end
cols = base(1:max(n, 0), :);
end
