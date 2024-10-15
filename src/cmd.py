# vim: set expandtab shiftwidth=4 softtabstop=4:

import requests
from chimerax.core.commands import CmdDesc, StringArg, BoolArg, FloatArg


def predict_class(session, file_path, rescale_cryoem=False, resolution=1.0):
    url = 'https://ligands.cs.put.poznan.pl/api/predict'
    params = {
        'rescale_cryoem': rescale_cryoem,
        'resolution': resolution,
    }
    print(f'loading files {file_path}')
    with open(file_path, 'rb') as ligand:
        session.logger.info('yep')
        r = requests.post(url, files={'file': ligand}, data=params)
    results = r.json()
    print(results)

# CmdDesc contains the command description.  For the
# "hello" command, we expect no arguments.
blob_desc = CmdDesc(required=[("file_path", StringArg)],
                    optional=[("rescale_cryoem", BoolArg),
                              ('resolution', FloatArg),
                              ], )
