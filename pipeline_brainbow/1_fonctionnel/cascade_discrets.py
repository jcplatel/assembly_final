# -*- coding: utf-8 -*-
"""
CASCADE_DISCRETS (19 nov. 2026, demandé par l'utilisateur)
Potentiels d'action discrets à partir de spike_prob (sortie de CASCADE),
avec cascade2p.utils_discrete_spikes.infer_discrete_spikes (réglages du
modèle : noyau de lissage 'auto'). Le placement des PA utilise des tirages
au hasard : graine fixée, une par paquet de cellules (résultat reproductible
pour un même nombre de processus).

Usage : python cascade_discrets.py entree.mat sortie.mat dossier_cascade [n_processus]
  entree.mat : spike_prob (cellules x images, NaN aux bords), modele, graine
  sortie.mat : n_pa (cellules x images, nombre de PA par image), modele
"""
import os
import sys
import numpy as np
import scipy.io as sio
from multiprocessing import Pool


def paquet(args):
    sp, modele, dossier, graine = args
    sys.path.insert(0, dossier)
    os.chdir(dossier)
    from cascade2p.utils_discrete_spikes import infer_discrete_spikes
    np.random.seed(graine)
    _, spikes = infer_discrete_spikes(sp, modele, model_folder=os.path.join(dossier, 'Pretrained_models'),
                                      verbosity=0)
    n = np.zeros(sp.shape, dtype=np.float32)
    for i, s in enumerate(spikes):
        for t in np.asarray(s, dtype=int):
            if 0 <= t < sp.shape[1]:
                n[i, t] += 1
    return n


if __name__ == '__main__':
    entree, sortie, dossier = sys.argv[1], sys.argv[2], sys.argv[3]
    n_proc = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    d = sio.loadmat(entree, squeeze_me=True)
    sp = np.atleast_2d(np.asarray(d['spike_prob'], dtype=float))
    modele = str(d['modele'])
    graine = int(d['graine']) if 'graine' in d else 1
    n_proc = max(1, min(n_proc, sp.shape[0]))
    morceaux = np.array_split(np.arange(sp.shape[0]), n_proc)
    taches = [(sp[m, :], modele, dossier, graine + k) for k, m in enumerate(morceaux)]
    if n_proc == 1:
        res = [paquet(taches[0])]
    else:
        with Pool(n_proc) as pool:
            res = pool.map(paquet, taches)
    n_pa = np.vstack(res)
    sio.savemat(sortie, {'n_pa': n_pa, 'modele': modele})
    print('PA discrets : %d cellules, %d PA, %d processus' % (sp.shape[0], int(n_pa.sum()), n_proc))
