# -*- coding: utf-8 -*-
"""
CASCADE_LANCER (18 nov. 2026, demandé par l'utilisateur)
Appelé par cascade_spikes.m. Lit un .mat (dF_traces : cellules x images, en
dF/F, non lissé ; modele : nom du modèle CASCADE), écrit un .mat avec
spike_prob (nombre de potentiels d'action estimé par image, NaN aux bords).
Usage : python cascade_lancer.py entree.mat sortie.mat dossier_cascade
"""
import sys
import os
import numpy as np
import scipy.io as sio

entree, sortie, dossier_cascade = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, dossier_cascade)
os.chdir(dossier_cascade)   # CASCADE cherche Pretrained_models/ dans le dossier courant
from cascade2p import cascade  # noqa: E402

D = sio.loadmat(entree)
traces = np.asarray(D['dF_traces'], dtype=np.float64)
modele = str(np.squeeze(D['modele']))
dossier_modeles = os.path.join(dossier_cascade, 'Pretrained_models')
if not os.path.isdir(os.path.join(dossier_modeles, modele)):
    cascade.download_model(modele, verbose=1)
spike_prob = cascade.predict(modele, traces, model_folder=dossier_modeles, verbosity=1)
sio.savemat(sortie, {'spike_prob': spike_prob.astype(np.float32), 'modele': modele})
print('CASCADE termine :', spike_prob.shape)
