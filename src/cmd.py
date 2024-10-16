import requests
from chimerax.core.commands import CmdDesc, StringArg, BoolArg, FloatArg


def predict_class(session, file_path, rescale_cryoem=False, resolution=1.0):
    # TODO: change file_path to chain or whatever

    models = session.models.list() # get all models in the session
    # Loop through models to find the PDB structure and extract the relevant information
    for i, model in enumerate(models):
        # print(i, model, dir(model)) # throw output of dir into ChatGPT to see if there's anything useful
        print("=====================================")
        if model.opened_data_format and model.opened_data_format.name == "mmCIF": # Check if the model is a PDB structure, hopefully
            pdb_id = model.name  # PDB ID or name of the model
            chains = model.chains  # List of chains in the PDB
            print(f"PDB ID: {pdb_id}")
            for chain in chains:
                print(f"Chain: id: {chain.chain_id}, structure: {chain.structure}, full name: {chain.full_name}")
        elif model.opened_data_format and model.opened_data_format.name == "CCP4 density map": # Check if the model is a density map, hopefully.
            print(dir(model.data))
            print(model.data.mrc_data)
        else: # opened_data_format is None
            print("??")

    # TODO: handle the case where there are multiple PDB structures in the session??


    # API call
    # url = 'https://ligands.cs.put.poznan.pl/api/predict'
    # params = {
    #     'rescale_cryoem': rescale_cryoem,
    #     'resolution': resolution,
    # }
    # print(f'loading files {file_path}')
    # with open(file_path, 'rb') as ligand:
    #     session.logger.info('yep')
    #     r = requests.post(url, files={'file': ligand}, data=params)
    # results = r.json()
    # print(results)


# CmdDesc contains the command description.  For the
# "hello" command, we expect no arguments.
blob_desc = CmdDesc(required=[("file_path", StringArg)],
                    optional=[("rescale_cryoem", BoolArg),
                              ('resolution', FloatArg),
                              ], )
