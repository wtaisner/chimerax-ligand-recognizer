import tempfile

import numpy as np
import requests
from chimerax.core.commands import CmdDesc, StringArg, BoolArg, FloatArg, run

from .cryoem_utils import cut_ligand


def pretty_print_predictions(predictions: list[dict[str, float]]) -> None:
    """
    Pretty print the ligand class predictions in a formatted way.

    :param predictions: A list of dictionaries containing 'Class' as the ligand class name
                    and 'Probability' as the associated probability.
    """
    print("Ligand Class Predictions:")
    print("=" * 30)
    for prediction in predictions:
        ligand_class = prediction['Class']
        probability = prediction['Probability']
        print(f"{ligand_class:30} | Probability: {probability:.3f}")
    print("=" * 30)

def validate_class(session, ligand_id) -> None:
    """ Prepare the ligand for the API and send it to the for validation.

    :param session: ChimeraX session
    :param ligand_id: id of the ligand to be validated
    """
    models = session.models.list()  # get all models in the session
    cif_model, map_model, residue = None, None, None
    # Loop through models to find the PDB structure and extract the relevant information
    for i, model in enumerate(models):
        if (
            model.opened_data_format and model.opened_data_format.name == "mmCIF"
        ):  # Check if the model is a PDB structure, hopefully

            # Check if there is more than one PDB structure in the session
            if cif_model is not None:
                print("Multiple PDB structures found in the session. Make sure only one PDB structure is open.")
            else:
                cif_model = model

            try:
                residue_command = f"select {ligand_id}"
                residue = run(session, residue_command)
            except Exception:
                print(f"Residue {ligand_id} not found in the structure. Please provide a valid id.")
                residue = None
        elif (
            model.opened_data_format
            and model.opened_data_format.name == "CCP4 density map"
        ):  # Check if the model is a density map, hopefully.
            map_model = model

    if map_model is None:
        print("Could not find density map. Please open a density map.")
    elif cif_model is None:
        print("Could not find PDB structure. Please open a PDB structure.")
    elif residue is None:
        print("Could not find ligand. Please provide a valid ligand id.")
    else:
        print("Attempting to cut ligand...")
        blob = cut_ligand(map_model, cif_model, residue)
        if blob is None:
            print("...could not cut ligand.")
        else:
            url = "https://ligands.cs.put.poznan.pl/api/predict"
            params = {"rescale_cryoem": False}
            with tempfile.NamedTemporaryFile(suffix=".npz", delete=True) as tmp_file:
                np.savez_compressed(tmp_file, blob)
                tmp_file.seek(0)
                r = requests.post(url, files={"file": tmp_file}, data=params)
            results = r.json()
            pretty_print_predictions(results['predictions'])

# CmdDesc contains the command description. It is used to register the command with ChimeraX.
blob_validate_desc = CmdDesc(
    required=[("ligand_id", StringArg)],
)
