# -*- coding: utf-8 -*-
"""
RECALER_AVEC_DECALAGES (18 nov. 2026, demandé par l'utilisateur)
Refait le film recalé d'un plan suite2p à partir des tifs bruts Bruker, en
réappliquant EXACTEMENT les décalages trouvés par suite2p à l'époque
(ops.yoff, xoff : rigide ; ops.yoff1, xoff1 : non rigide, dans Fall.mat).
Pas de nouvelle estimation : le film retombe sur les masques et les traces
existants. Sert quand le film recalé (data.bin) a été perdu et que les
exports (MotCorr.ome.tif, NWB) sont incomplets.

Les tifs bruts : un fichier par cycle (..._CycleNNNNN_Ch2_000001.ome.tif),
avec une page par plan (nplanes pages). L'ordre des cycles est celui de
ops.filelist ; le plan p = page p de chaque fichier (suite2p lit les pages
à la suite et les répartit entre les plans).
Comme suite2p, les uint16 sont divisés par 2 et passés en int16.

Usage :
  python recaler_avec_decalages.py Fall.mat dossier_TSeries sortie.bin
Écrit sortie.bin (int16, images Ly x Lx à la suite, comme data.bin de suite2p).
"""
import os
import sys
import numpy as np
import scipy.io as sio
import tifffile
from suite2p.registration import register, nonrigid
from suite2p.default_ops import default_ops

f_fall, dossier, sortie = sys.argv[1], sys.argv[2], sys.argv[3]
n_test = int(sys.argv[4]) if len(sys.argv) > 4 else 0   # test : seulement les n_test premières images
o = sio.loadmat(f_fall, variable_names=['ops'], squeeze_me=True, struct_as_record=False)['ops']
ops = default_ops()
for k in ['nplanes', 'iplane', 'Ly', 'Lx', 'nonrigid', 'block_size', 'bidiphase', 'nframes',
          'maxregshiftNR', 'smooth_sigma', 'snr_thresh']:
    if hasattr(o, k):
        ops[k] = getattr(o, k)
ops['bidi_corrected'] = True          # do_bidiphase = 0 à l'époque
ops['bidiphase'] = 0
nplanes = int(o.nplanes)
Ly, Lx = int(o.Ly), int(o.Lx)
ops['block_size'] = [int(b) for b in np.atleast_1d(o.block_size)]
yoff = np.asarray(o.yoff).astype(int)
xoff = np.asarray(o.xoff).astype(int)
yoff1 = np.asarray(o.yoff1, dtype=np.float32)
xoff1 = np.asarray(o.xoff1, dtype=np.float32)
nframes = int(o.nframes)

# plan de ce Fall.mat : d'après le chemin de reg_file (.../plane1/data.bin)
reg = str(o.reg_file).replace('\\', '/')
plan = int(reg.split('/plane')[-1].split('/')[0])

# blocs non rigides : ceux de l'époque (ops.yblock, xblock, nblocks), vérifiés
yblock, xblock, nblocks, block_size, _ = nonrigid.make_blocks(Ly=Ly, Lx=Lx, block_size=ops['block_size'])
yb0 = np.asarray(o.yblock); xb0 = np.asarray(o.xblock)
if not (np.array_equal(np.asarray(yblock), yb0) and np.array_equal(np.asarray(xblock), xb0)):
    print('blocs recalculés différents de ceux de ops : on garde ceux de ops', flush=True)
    yblock = [yb0[i] for i in range(yb0.shape[0])]
    xblock = [xb0[i] for i in range(xb0.shape[0])]
    nblocks = [int(v) for v in np.atleast_1d(o.nblocks)]
blocks = (yblock, xblock, nblocks)

fichiers = [os.path.join(dossier, os.path.basename(str(f).replace('\\', '/'))) for f in np.atleast_1d(o.filelist)]
if n_test:
    fichiers = fichiers[:n_test]
    nframes = n_test
manquants = [f for f in fichiers if not os.path.isfile(f)]
if manquants:
    sys.exit('Fichiers bruts manquants : %d (premier : %s)' % (len(manquants), manquants[0]))
print('plan %d, %d fichiers, %d images attendues' % (plan, len(fichiers), nframes), flush=True)

# images du plan, dans l'ordre (toutes les pages à la suite, réparties entre les plans)
def images_du_plan():
    n = 0
    for f in fichiers:
        with tifffile.TiffFile(f) as t:
            for page in t.pages:
                if n % nplanes == plan:
                    yield page.asarray()
                n += 1

lot = 500
tampon = []
i0 = 0
with open(sortie, 'wb') as fo:
    def ecrire(tampon, i0):
        fr = np.stack(tampon)
        if fr.dtype == np.uint16:
            fr = (fr // 2).astype(np.int16)
        else:
            fr = fr.astype(np.int16)
        k = slice(i0, i0 + fr.shape[0])
        fr = register.shift_frames(fr, yoff[k], xoff[k], yoff1[k], xoff1[k], blocks, ops)
        fo.write(np.asarray(fr, dtype=np.int16).tobytes())
        return i0 + fr.shape[0]
    for im in images_du_plan():
        tampon.append(im)
        if len(tampon) == lot:
            i0 = ecrire(tampon, i0); tampon = []
            print('  %d / %d' % (i0, nframes), flush=True)
    if tampon:
        i0 = ecrire(tampon, i0)
if i0 != nframes:
    sys.exit('ATTENTION : %d images écrites, %d attendues' % (i0, nframes))
print('Fini : %d images -> %s' % (i0, sortie))
