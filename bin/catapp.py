#!/usr/bin/env python3

"""
Query and save reaction data from the CatApp database. The CatApp database is
downloaded in full and used locally.
"""

import json
import os
import requests
import ase.db


def catapp(root_dir):
    # If the database file is not found, download it
    db_path = f"{root_dir}/data/catapp.db"
    if not os.path.isfile(db_path):
        print("Database file not found. Downloading... ")
        url = "https://cmr.fysik.dtu.dk/_downloads/716b1e0826acbb3d80675c116a2cb8a6/catapp.db"
        db_file = requests.get(url)
        with open(db_path, "wb") as f:
            f.write(db_file.content)
        print("Done")

    # Connect to the database
    con = ase.db.connect(db_path)

    # Save the reactions into a dictionary
    reactions = {}

    # Iterate over the database and save each feature into
    # a row, and then append that row into the dictionary.
    count = 0
    for row in con.select():
        try:
            reaction = {}
            reaction["reactant_a"] = row.a
            reaction["reactant_b"] = row.b
            reaction["product_ab"] = row.ab
            reaction["reactant_surface"] = row.surface
            reaction["reactant_facet"] = row.facet
            try:
                # The adsorption site does not always exist
                reaction["site"] = row.site
            except AttributeError:
                # Assign a string "None" if value not found
                reaction["site"] = "None"
            reaction["reaction_energy"] = row.er
            reaction["activation_energy"] = row.ea
            reactions[count] = reaction
            count += 1
        except AttributeError:
            # Crude way to handle unexpected errors:
            # simply discard the reaction
            pass

    # Save the reaction dictionary into Json file
    write_file = f"{root_dir}/data/reactions_catapp.json"
    with open(write_file, "w") as f:
        json.dump(reactions, f)

    print(f"Saved {count+1} reactions into file data/reactions_catapp.json")
