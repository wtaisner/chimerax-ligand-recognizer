import tempfile

import numpy as np
import requests
from chimerax.atomic import AtomicStructure
from chimerax.core.commands import CmdDesc, StringArg, run
from chimerax.geometry import Bounds
from chimerax.map import Volume, VolumeSurface
from chimerax.map_data import ArrayGridData
from chimerax.map_data.ccp4.ccp4_grid import CCP4Grid
from chimerax.mask.maskcommand import mask, ones_mask
from chimerax.pick_blobs.pickblobs import BlobOutlineBox

from .cryoem_utils import cut_ligand, read_map, cut_ligand_v2
from .utils import pretty_print_predictions


def validate_class(session, ligand_id) -> None:
    """ Prepare the ligand for the API and send it to the for validation.

    :param session: ChimeraX session
    :param ligand_id: id of the ligand to be validated
    """
    session.logger.info(msg="Attempting to cut ligand (this may take a while)...")
    models = session.models.list()  # get all models in the session
    cif_model, map_model, residue = None, None, None
    # Loop through models to find the PDB structure and extract the relevant information
    for i, model in enumerate(models):
        if (
                model.opened_data_format and model.opened_data_format.name == "mmCIF"
        ):  # Check if the model is a PDB structure, hopefully

            # Check if there is more than one PDB structure in the session
            if cif_model is not None:
                session.logger.error(
                    msg="Multiple PDB structures found in the session. Make sure only one PDB structure is open.")
            else:
                cif_model = model

            try:
                residue_command = f"select {ligand_id}"
                residue = run(session, residue_command, log=False)
            except Exception:
                session.logger.error(f"Residue {ligand_id} not found in the structure. Please provide a valid id.")
                residue = None
        elif (
                model.opened_data_format
                and model.opened_data_format.name == "CCP4 density map"
        ):  # Check if the model is a density map, hopefully.
            map_model = model

    if map_model is None:
        session.logger.error("Could not find density map. Please open a density map.")
    elif cif_model is None:
        session.logger.error("Could not find PDB structure. Please open a PDB structure.")
    elif residue is None:
        session.logger.error("Could not find ligand. Please provide a valid ligand id.")
    else:
        blob = cut_ligand(map_model, cif_model, residue)
        if blob is None:
            session.logger.error("...could not cut ligand.")
        else:
            url = "https://ligands.cs.put.poznan.pl/api/predict"
            params = {"rescale_cryoem": False}
            with tempfile.NamedTemporaryFile(suffix=".npz", delete=True) as tmp_file:
                np.savez_compressed(tmp_file, blob)
                tmp_file.seek(0)
                r = requests.post(url, files={"file": tmp_file}, data=params)
            results = r.json()
            session.logger.info(msg=pretty_print_predictions(results['predictions']), is_html=True)


# CmdDesc contains the command description. It is used to register the command with ChimeraX.
blob_validate_desc = CmdDesc(
    required=[("ligand_id", StringArg)]
)
blobus_validatus_desc = CmdDesc(
    required=[("ligand_id", StringArg)]
)

def recognize_class(session, ligand_id) -> None:
    """

    :param session:
    :param ligand_id:
    :return:
    """
    models = session.models.list()  # get all models in the session
    print(models)
    cif_model, map_model, residue, blob = None, None, None, None
    # Loop through models to find the PDB structure and extract the relevant information
    for i, model in enumerate(models):
        print(i, type(model), model.opened_data_format)
        if model.opened_data_format and model.opened_data_format.name == "mmCIF":  # Check if the model is a PDB structure, hopefully
            # Check if there is more than one PDB structure in the session
            if cif_model is not None:
                session.logger.error(
                    msg="Multiple PDB structures found in the session. Make sure only one PDB structure is open.")
            else:
                cif_model: AtomicStructure = model

            try:
                residue_command = f"select {ligand_id}"
                residue = run(session, residue_command, log=False)
            except Exception:
                session.logger.error(f"Residue {ligand_id} not found in the structure. Please provide a valid id.")
                residue = None
        elif model.opened_data_format and model.opened_data_format.name == "CCP4 density map":  # Check if the model is a density map, hopefully.
            map_model: Volume = model
        elif isinstance(model, BlobOutlineBox):
            blob: BlobOutlineBox = model

    if map_model is None:
        session.logger.error("Could not find density map. Please open a density map.")
    elif cif_model is None:
        session.logger.error("Could not find PDB structure. Please open a PDB structure.")
    elif residue is None:
        session.logger.error("Could not find ligand. Please provide a valid ligand id.")
    else:
        print("======================================")
        # new_volume: Volume = mask(session, volumes=[map_model], surfaces=map_model.surfaces)[0]
        # new_volume_blob: Volume = mask(session, volumes=[map_model], surfaces=[blob])[0]
        # print(new_volume.bounds(), new_volume_blob.bounds())
        new_volume_ones = ones_mask(session, surfaces=map_model.surfaces)[0]
        print("======================================")
        print(map_model, new_volume_ones)
        grid_data: ArrayGridData = new_volume_ones.data
        print(grid_data.full_matrix().shape)
        # CCP4_data: CCP4Grid = map_model.data
        setattr(new_volume_ones.data, "file_header", map_model.data.file_header)
        print(type(map_model.data), type(new_volume_ones.data))

        blob_v2 = cut_ligand_v2(map_model, new_volume_ones)
        if blob_v2 is None:
            session.logger.error("...could not cut ligand.")
        else:
            url = "https://ligands.cs.put.poznan.pl/api/predict"
            params = {"rescale_cryoem": False}
            with tempfile.NamedTemporaryFile(suffix=".npz", delete=True) as tmp_file:
                np.savez_compressed(tmp_file, blob)
                tmp_file.seek(0)
                r = requests.post(url, files={"file": tmp_file}, data=params)
            results = r.json()
            session.logger.info(msg=pretty_print_predictions(results['predictions']), is_html=True)


blob_recognize_desc = CmdDesc(
    required=[("ligand_id", StringArg)]
)
blobus_recognitus_desc = CmdDesc(
    required=[("ligand_id", StringArg)]
)
