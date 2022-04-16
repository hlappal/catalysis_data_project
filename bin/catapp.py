# Imports
import matplotlib.pyplot as plt
import numpy as np
import ase.db
import requests
import json
import os

from ase.visualize import view
from sklearn.kernel_ridge import KernelRidge as KRR
from sklearn.neighbors import KNeighborsRegressor as KNN

# Define the project root directory
# ROOT_DIR = os.path.join(os.getcwd(), os.pardir)
ROOT_DIR = os.getcwd()

def main():
    # Download the database, if it does not yeat exist in the root directory
    if not os.path.isfile(f"{ROOT_DIR}/data/catapp.db"):
        url = "https://cmr.fysik.dtu.dk/_downloads/716b1e0826acbb3d80675c116a2cb8a6/catapp.db"
        db_file = requests.get(url)
        with open(f"{ROOT_DIR}/data/catapp.db", "wb") as file:
            file.write(db_file.content)
        print(f"The database downloaded as {ROOT_DIR}/data/catapp.db")

    # Connect to database
    con = ase.db.connect(f"{ROOT_DIR}/data/catapp.db")

    # Save the reactions into a dictionary
    reactions = {}

    # Iterate over the database and save each feature into
    # a row, and then append that row into the dictionary.
    i = 0
    for row in con.select():
        try:
            reaction = {}
            reaction['reactant_a'] = row.a
            reaction['reactant_b'] = row.b
            reaction['product_ab'] = row.ab
            reaction['reactant_surface'] = row.surface
            reaction['reactant_facet'] = row.facet
            try:
                # The adsorption site does not always exist
                reaction['site'] = row.site
            except AttributeError:
                # Assign a string 'None' if value not found
                reaction['site'] = 'None'
            reaction['reaction_energy'] = row.er
            reaction['activation_energy'] = row.ea
            reaction['dft_functional'] = row.xc
            reactions[i] = reaction
            i += 1
        except AttributeError:
            # A crude way to handle unexpected errors:
            # The reaction is simply discarded
            pass

    # Save the reaction dictionary into a Json file
    with open(f"{ROOT_DIR}/data/reactions_catapp.json", "w") as file:
        json.dump(reactions, file)

    print('All reactions read into file')


if __name__ == '__main__':
    main()