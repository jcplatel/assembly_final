function spike_prob = cascade_spikes(DFF, opts, fichier_cache)
%% =========================================================================
%  CASCADE_SPIKES (18 nov. 2026, demandée par l'utilisateur)
%  Déconvolution par CASCADE (Rupprecht et al. 2021, modèles GCaMP8) des
%  traces dF/F (cellules x images, NON lissées). Rend spike_prob : nombre de
%  potentiels d'action estimé par image (pas une probabilité ; peut dépasser
%  1) ; NaN sur les 32 premières et dernières images (fenêtre de 64).
%  CASCADE tourne en Python (cascade_lancer.py, même dossier).
%  opts (champs facultatifs) :
%    cascade_modele  : 'GC8_EXC_10Hz_smoothing150ms_high_noise' (données
%                      GCaMP8 à 10 Hz ; nos sessions sont à 10,4 Hz)
%    cascade_python  : python de l'environnement CASCADE
%    cascade_dossier : dossier du code CASCADE (GitHub HelmchenLabSoftware)
%  fichier_cache (facultatif) : si le fichier existe et a été calculé avec le
%  même modèle et la même taille, il est relu ; sinon il est écrit.
%  =========================================================================
if nargin < 2 || isempty(opts), opts = struct(); end
if nargin < 3, fichier_cache = ''; end
modele = lire(opts, 'cascade_modele', 'GC8_EXC_10Hz_smoothing150ms_high_noise');
python = lire(opts, 'cascade_python', 'C:\Users\jcplatel\Desktop\matlab\Cascade\venv_cascade\Scripts\python.exe');
dossier = lire(opts, 'cascade_dossier', 'C:\Users\jcplatel\Desktop\matlab\Cascade');

controle = [sum(double(DFF(:)), 'omitnan'), sum(abs(double(DFF(:))), 'omitnan')];   % même DFF ?
if ~isempty(fichier_cache) && isfile(fichier_cache)
    C = load(fichier_cache, 'spike_prob', 'modele', 'controle');
    if strcmp(C.modele, modele) && isequal(size(C.spike_prob), size(DFF)) && ...
            isfield(C, 'controle') && all(abs(C.controle - controle) <= 1e-6 * max(1, abs(controle)))
        spike_prob = double(C.spike_prob);
        return
    end
end

script = fullfile(fileparts(mfilename('fullpath')), 'cascade_lancer.py');
base = tempname;
entree = [base '_entree.mat']; sortie = [base '_sortie.mat'];
dF_traces = double(DFF); %#ok<NASGU>
save(entree, 'dF_traces', 'modele', '-v7');
% bibliothèques SSL / DLL de l'environnement conda d'origine du venv
env_conda = 'C:\Users\jcplatel\.conda\envs\suite2p';
ancien_path = getenv('PATH');
setenv('PATH', [env_conda '\Library\bin;' env_conda '\DLLs;' ancien_path]);
nettoyer = onCleanup(@() setenv('PATH', ancien_path));
cmd = sprintf('"%s" "%s" "%s" "%s" "%s"', python, script, entree, sortie, dossier);
[st, txt] = system(cmd);
delete(entree);
if st ~= 0 || ~isfile(sortie)
    error('cascade_spikes : CASCADE a échoué :\n%s', txt);
end
S = load(sortie, 'spike_prob'); delete(sortie);
spike_prob = double(S.spike_prob);
if ~isempty(fichier_cache)
    sauver_cache(fichier_cache, single(spike_prob), modele, controle);
end
end

function sauver_cache(f, spike_prob, modele, controle) %#ok<INUSD>
save(f, 'spike_prob', 'modele', 'controle', '-v7.3');
end

function v = lire(s, champ, defaut)
if isfield(s, champ) && ~isempty(s.(champ)), v = s.(champ); else, v = defaut; end
end
