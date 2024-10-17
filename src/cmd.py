import tempfile

import numpy as np
import requests
from chimerax.core.commands import CmdDesc, StringArg, BoolArg, FloatArg, run

from .cryoem_utils import cut_ligand


def validate_class(session, ligand_id):

    models = session.models.list()  # get all models in the session
    cif_model, map_model, residue = None, None, None
    # Loop through models to find the PDB structure and extract the relevant information
    for i, model in enumerate(models):
        if (
            model.opened_data_format and model.opened_data_format.name == "mmCIF"
        ):  # Check if the model is a PDB structure, hopefully
            cif_model = model
            try:
                residue_command = f"select {ligand_id}"
                residue = run(session, residue_command)
            except Exception:
                residue = None
        elif (
            model.opened_data_format
            and model.opened_data_format.name == "CCP4 density map"
        ):  # Check if the model is a density map, hopefully.
            map_model = model

        print("=====================================")
    if map_model is None:
        print("Could not find density map")
    elif cif_model is None:
        print("Could not find cif file")
    elif residue is None:
        print(f"Incorrect ligand id: {ligand_id}")
    else:
        blob = cut_ligand(map_model, cif_model, residue)
        if blob is None:
            print("Could not cut ligand")
        else:
            url = "https://ligands.cs.put.poznan.pl/api/predict"
            params = {"rescale_cryoem": False}
            with tempfile.NamedTemporaryFile(suffix=".npz", delete=True) as tmp_file:
                np.savez_compressed(tmp_file, blob)
                tmp_file.seek(0)
                r = requests.post(url, files={"file": tmp_file}, data=params)
            results = r.json()
            print(results)
    # TODO: handle the case where there are multiple PDB structures in the session??


# CmdDesc contains the command description.  For the
# "hello" command, we expect no arguments.
blob_validate_desc = CmdDesc(
    required=[("ligand_id", StringArg)],
    optional=[
        ("rescale_cryoem", BoolArg),
        ("resolution", FloatArg),
    ],
)
