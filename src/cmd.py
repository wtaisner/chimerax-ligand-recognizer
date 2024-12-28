import tempfile

import numpy as np
import requests
from chimerax.atomic import AtomicStructure
from chimerax.core.commands import CmdDesc, StringArg, run, BoolArg, FloatArg, IntArg
from chimerax.core.errors import UserError
from chimerax.core.models import Surface
from chimerax.core.objects import Objects
from chimerax.core.session import Session
from chimerax.map import Volume, VolumeSurface
from chimerax.map.volumecommand import volume
from chimerax.mask.maskcommand import mask
from scipy.stats import norm

from .cryoem_utils import cut_ligand_from_coords, cut_ligands_by_hand, em_stats, read_map
from .utils import pretty_print_predictions


def get_model(session: Session, model_id: str):
    """Selects a model from ChimeraX session given the model id

    :param session: ChimeraX session
    :param model_id: id of the model (e.g., PDB file, density map)
    :return: ChimeraX model or None
    """
    model_tuple = tuple([int(x) for x in model_id[1:].split('.')])
    if not model_tuple in session.models._models:
        raise UserError(f"Could not find queried id: {model_id}")
    model = session.models._models[model_tuple]
    return model


def validate_blob(session: Session, blob: np.ndarray) -> None:
    """ Perform the API call and print HTML output.

    :param session: ChimeraX session
    :param blob: Blob data as np.ndarray
    """
    if blob is None:
        session.logger.warning("...could not cut ligand.")
    else:
        url = "https://ligands.cs.put.poznan.pl/api/predict"
        params = {"rescale_cryoem": False}
        with tempfile.NamedTemporaryFile(suffix=".npz", delete=True) as tmp_file:
            np.savez_compressed(tmp_file, blob)
            tmp_file.seek(0)
            r = requests.post(url, files={"file": tmp_file}, data=params)
        results = r.json()
        session.logger.info(msg=pretty_print_predictions(results['predictions']), is_html=True)


def validate_class(
        session: Session,
        res_id: str,
        map_id: str | None = None,
        pdb_id: str | None = None,
        flg_xray: bool = False,
        resolution: float | None = None,
        density_threshold: float | None = None,
) -> None:
    """ Prepare the ligand for the API and send it to the for validation.
    :param session: ChimeraX session
    :param res_id: id of the ligand to be validated
    :param map_id: id of the density map
    :param pdb_id: id of the PDB structure
    :param flg_xray: whether the density map comes from xray crystallography or not (if not it is assumed to be cryoem density map)
    :param resolution: resolution of the density map (not recommended - a preferred option is to use PDB file), if pdb_id is given resolution parameter will be taken from PDB file (if found)
    :param density_threshold: threshold value for the density map
    """
    if resolution is not None:
        session.logger.warning(
            "Resolution provided by hand. A recommended option is to use a PDB file with resolution information. See tutorial for details.")

    if pdb_id is None and res_id[0] == '#':
        pdb_id = res_id.split('/')[0]

    cif_model, map_model, residue = None, None, None
    if pdb_id is not None:
        cif_model = get_model(session, pdb_id)
        if not isinstance(cif_model, AtomicStructure) or not (
                cif_model.opened_data_format and (cif_model.opened_data_format.name == "mmCIF" or cif_model.opened_data_format.name == "PDB")):
            raise UserError(f"Expected the id {pdb_id} to refer to PDB structure")

    if map_id is not None:
        map_model = get_model(session, map_id)
        if not isinstance(map_model, Volume) or not (
                map_model.opened_data_format and map_model.opened_data_format.name == 'CCP4 density map'):
            raise UserError(f"Expected the id {map_id} to refer to CCP4 density map")

    residue: Objects = run(session, f"select {res_id}",
                           log=False)  # it will always return an Object, even if it is empty

    if pdb_id is None or map_id is None:
        models = session.models.list()  # get all models in the session
        # Loop through models to find the PDB structure and extract the relevant information
        for model in models:
            if model is None: # skip None objects
                continue
            if model.opened_data_format and (model.opened_data_format.name == "mmCIF" or model.opened_data_format.name == "PDB"):  # Check if the model is a PDB structure, hopefully

                # Check if there is more than one PDB structure in the session
                if cif_model is not None and pdb_id is None:
                    raise UserError(
                        "Multiple PDB structures found in the session. Please provide id of the PBD structure.")
                if cif_model is None and pdb_id is None:
                    cif_model = model

            if model.opened_data_format and model.opened_data_format.name == "CCP4 density map":  # Check if the model is a density map, hopefully.
                if map_model is not None and map_id is None:
                    raise UserError("Multiple density maps found in the session. Please provide id of the density map.")
                if map_model is None and map_id is None:
                    map_model = model

    if map_model is None:
        raise UserError("Could not find density map. Please open a density map or provide a valid density map id.")
    if cif_model is None:
        raise UserError(
            "Could not find PDB structure. Please open a PDB structure or provide a valid PDB structure id.")
    if residue.num_atoms == 0:
        raise UserError(
            f"Residue {res_id} not found in the structure. Please provide a valid ligand id.")

    blob = cut_ligand_from_coords(map_model, cif_model, residue, flg_xray, density_threshold=density_threshold, resolution=resolution)
    validate_blob(session, blob)


blob_validate_desc = CmdDesc(
    required=[("res_id", StringArg)],
    optional=[("map_id", StringArg), ("pdb_id", StringArg), ("flg_xray", BoolArg), ("density_threshold", FloatArg), ("resolution", FloatArg)]
)
blobus_validatus_desc = CmdDesc(
    required=[("res_id", StringArg)],
    optional=[("map_id", StringArg), ("pdb_id", StringArg), ("flg_xray", BoolArg), ("density_threshold", FloatArg), ("resolution", FloatArg)]
)


def recognize_class(
        session: Session,
        map_id: str | None = None,
        surface_id: str | None = None,
        pdb_id: str | None = None,
        flg_xray: bool = False,
        resolution: float | None = None,
        density_threshold: float | None = None
) -> None:
    """ Recognize extracted part of a density map.

     This method creates intermediate mask object,
        characterized by the "copy masked" suffix, which are removed at the end of the method.
    :param session: ChimeraX Session object
    :param map_id: id of entire density map in ChimeraX
    :param pdb_id: id of the PDB structure in ChimeraX
    :param surface_id: id of surface object that one wishes to recognize, if not given defaults to the surface of density map object
    :param flg_xray: whether the density map comes from xray crystallography or not (if not it is assumed to be cryoem density map)
    :param resolution: resolution of the density map (not recommended - a preferred option is to use PDB file), if pdb_id is given resolution parameter will be taken from PDB file (if found)
    :param density_threshold: threshold value for the density map
    """
    if resolution is not None:
        session.logger.warning(
            "Resolution provided by hand. A recommended option is to use a PDB file with resolution information. See tutorial for details.")

    map_model, blob_model, cif_model = None, None, None
    if map_id is not None:
        map_model = get_model(session, map_id)
        if not isinstance(map_model, Volume) or not (
                map_model.opened_data_format and map_model.opened_data_format.name == 'CCP4 density map'):
            raise UserError(f"Expected the id {map_id} to refer to CCP4 density map")

    if surface_id is not None:
        blob_model = get_model(session, surface_id)
        if not isinstance(blob_model, (Surface, VolumeSurface)):
            raise UserError(f"Expected the id {surface_id} to refer to surface data")

    if pdb_id is not None:
        cif_model = get_model(session, pdb_id)
        if not isinstance(cif_model, AtomicStructure) or not (
                cif_model.opened_data_format and (cif_model.opened_data_format.name == "mmCIF" or cif_model.opened_data_format.name == "PDB")):
            raise UserError(f"Expected the id {pdb_id} to refer to PDB structure")

    if map_model is None or (cif_model is None and resolution is None):
        models = session.models.list()  # get all models in the session
        # Loop through models to find the PDB structure and extract the relevant information
        for model in models:
            if model.opened_data_format and model.opened_data_format.name == "CCP4 density map":  # Check if the model is a density map, hopefully.
                if map_model is not None and map_id is None:
                    raise UserError("Multiple density maps found in the session. Please provide id of the density map.")
                if map_model is None and map_id is None:
                    map_model = model

            if resolution is None:
                if model.opened_data_format and (model.opened_data_format.name == "mmCIF" or model.opened_data_format.name == "PDB"):  # Check if the model is a PDB structure, hopefully
                    # Check if there is more than one PDB structure in the session
                    if cif_model is not None and pdb_id:
                        raise UserError(
                            "Multiple PDB structures found in the session. Please provide id of the PBD structure.")
                    if cif_model is None and pdb_id is None:
                        cif_model = model

    if surface_id is None and map_model is not None:
        blob_model = map_model.surfaces[0]

    if cif_model is not None:
        resolution_cif, _ = em_stats(cif_model)
        resolution = resolution_cif if resolution_cif is not None else resolution

    if resolution is None:
        raise UserError(
            "Could not find resolution. Please make sure that it is either defined in the PDB file, or has been passed manually.")
    if map_model is None:
        raise UserError("Could not find density map. Please open density map or provide a valid density map id.")
    if blob_model is None:
        raise UserError("Could not find surface. Please provide a valid surface id.")

    if map_model is not None and blob_model is not None and resolution is not None:
        map_model_ones = map_model.writable_copy(require_copy=True, subregion='all', open_model=False,
                                                 unshow_original=False)
        map_model_ones.data.full_matrix()[:] = 1

        blob_mask = mask(session, volumes=[map_model_ones], surfaces=[blob_model], full_map=True)[0]
        setattr(blob_mask.data, "file_header", map_model.data.file_header)

        blob = cut_ligands_by_hand(
            map_model=map_model,
            mask_model=blob_mask,
            resolution=resolution,
            xray=flg_xray,
            density_threshold=density_threshold
        )

        validate_blob(session, blob)

        # cleanup - remove the intermediate mask object
        models = session.models.list()
        for model in models:
            if model.name.endswith("copy masked"):
                run(session, f"close #{'.'.join(map(str, model.id))}", log=False)


blob_recognize_desc = CmdDesc(
    optional=[("map_id", StringArg), ("surface_id", StringArg), ("pdb_id", StringArg), ("flg_xray", BoolArg),
              ("resolution", FloatArg), ("density_threshold", FloatArg)],
)
blobus_recognitus_desc = CmdDesc(
    optional=[("map_id", StringArg), ("surface_id", StringArg), ("pdb_id", StringArg), ("flg_xray", BoolArg),
              ("resolution", FloatArg), ("density_threshold", FloatArg)],
)


def blob_autothreshold(session: Session, map_id: str | None = None, style: int = 0,
                       density_std_threshold: float = 2.8) -> None:
    """
    Automatically set the level of the density map to the value corresponding to the computed density threshold.
    :param session: ChimeraX session
    :param map_id: id of the density map
    :param style: 0->unchanged, 1->stick, 2->ribbon
    :param density_std_threshold:

    :return: None
    """
    map_model = None
    if map_id is not None:
        map_model = get_model(session, map_id)
        if not isinstance(map_model, Volume) or not (
                map_model.opened_data_format and map_model.opened_data_format.name == 'CCP4 density map'):
            raise UserError(f"Expected the id {map_id} to refer to CCP4 density map")
    if map_model is None:
        models = session.models.list(type=Volume)  # get all Volumes in the session
        for model in models:
            if model.opened_data_format and model.opened_data_format.name == "CCP4 density map":  # Check if the model is a density map, hopefully.
                if map_model is not None and map_id is None:
                    raise UserError("Multiple density maps found in the session. Please provide id of the density map.")
                if map_model is None and map_id is None:
                    map_model = model

    _, map_array, _, _ = read_map(map_model)

    map_median = np.median(map_array)
    map_std = np.std(map_array)
    value_mask = (map_array < map_median - 0.5 * map_std) | (
            map_array > map_median + 0.5 * map_std
    )

    quantile_threshold = norm.cdf(density_std_threshold)
    density_threshold = np.quantile(map_array[value_mask], quantile_threshold)
    if style < 0 or style > 2:
        raise UserError("Style must be 0, 1, or 2!")

    if style == 0:
        volume(session=session, volumes=[map_model], level=[[density_threshold]])
    else:
        volume(session=session, volumes=[map_model], level=[[density_threshold]], style="surface", transparency=0.5,
               step=1)
        if style == 1:
            run(session, "hide ribbons; show atoms; style stick")
        else:
            run(session, "hide atoms; show ribbons")

    session.logger.info(f"Level set to {density_threshold:.4f} for density map #{'.'.join(map(str, map_model.id))}")


blob_autothreshold_desc = CmdDesc(
    optional=[("map_id", StringArg), ("style", IntArg), ("density_std_threshold", FloatArg)]
)
blobus_autothreshold_desc = CmdDesc(
    optional=[("map_id", StringArg), ("style", IntArg), ("density_std_threshold", FloatArg)]
)
