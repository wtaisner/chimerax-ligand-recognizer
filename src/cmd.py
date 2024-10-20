import tempfile

import numpy as np
import requests
from chimerax.core.commands import CmdDesc, StringArg, run
from chimerax.core.session import Session
from chimerax.map import Volume, VolumeSurface
from chimerax.mask.maskcommand import mask

from .cryoem_utils import cut_ligand, cut_and_extract_ligand
from .utils import pretty_print_predictions


def validate_blob(session: Session, blob: np.ndarray) -> None:
    """ Perform the API call and print HTML output.

    :param session: ChimeraX session
    :param blob: Blob data as np.ndarray
    """
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


def validate_class(session: Session, ligand_id: str) -> None:
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

        validate_blob(session, blob)


blob_validate_desc = CmdDesc(
    required=[("ligand_id", StringArg)]
)
blobus_validatus_desc = CmdDesc(
    required=[("ligand_id", StringArg)]
)

def recognize_class(session: Session, map_id: str, surface_id: str) -> None:
    """ Recognize extracted part of a density map
    :param session: ChimeraX Session object
    :param map_id: id of entire density map in ChimeraX
    :param surface_id: id of surface object that one wishes to recognize
    """
    map_tuple = tuple([int(x) for x in map_id[1:].split('.')])
    surface_tuple = tuple([int(x) for x in surface_id[1:].split('.')])

    if not map_tuple in session.models._models:
        session.logger.error(f"Could not find queried map model")
    map_model = session.models._models[map_tuple]
    if not isinstance(map_model, Volume) or not (map_model.opened_data_format and map_model.opened_data_format.name == 'CCP4 density map'):
        session.logger.error(f"Expected the id to refer to CCP4 density map")
    if not surface_tuple in session.models._models:
        session.logger.error(f"Could not find queried surface")
        map_model = None
    blob_model = session.models._models[surface_tuple]
    if not isinstance(blob_model, VolumeSurface):
        session.logger.error(f"Expected the id to refer to surface data")
        blob_model = None

    if map_model is not None and blob_model is not None:
        map_model_ones = map_model.writable_copy(require_copy=True, subregion='all', open_model=False)
        map_model_ones.data.full_matrix()[:] = 1

        blob_mask = mask(session, volumes=[map_model_ones], surfaces=[blob_model], full_map=True)[0]
        setattr(blob_mask.data, "file_header", map_model.data.file_header)
        blob = cut_and_extract_ligand(map_model, blob_mask)

        validate_blob(session, blob)


blob_recognize_desc = CmdDesc(
    required=[("map_id", StringArg), ("surface_id", StringArg)]
)
blobus_recognitus_desc = CmdDesc(
    required=[("map_id", StringArg), ("surface_id", StringArg)]
)
