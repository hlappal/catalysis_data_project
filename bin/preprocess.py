#!/usr/bin/env python3

"""
Functions for the catalysis reaction data preprocessing. The functions are
quite ad hoc, and suited only for the raw data acquired from the
Catalysis-hub and the CatApp catalysis reaction databases.

Author:  Heikki Lappalainen
E-mail:  heikki.lappalainen@protonmail.com
"""

import json
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from dscribe.descriptors import LMBTR
from ase.io import read
from ase import Atoms


def preprocess(root_dir):
    """
    Process the data into suitable form for the ML models.
    """

    df_cathub_raw = load_json(f"{root_dir}/data/reactions_cathub.json")
    df_catapp_raw = load_json(f"{root_dir}/data/reactions_catapp.json")

    # Get the max number of reactants and products
    max_reactants, max_products = max_reactants_products(
        df_cathub_raw, 0, 0)

    # Parse the reaction data
    df_cathub = parse_reactants_products(df_cathub_raw, "Catalysis-hub",
                                         max_reactants, max_products)
    df_catapp = parse_reactants_products(df_catapp_raw, "CatApp",
                                         max_reactants, max_products)

    # Combine the datasets
    df_raw = df_cathub.append(df_catapp, ignore_index=True, sort=False)

    # Clean the data
    df = clean_data(df_raw)

    # Save the data
    df.to_pickle(f"{root_dir}/data/data.csv")
    print(f"Saved the processed data into file {root_dir}/data/data.csv")

    return df


def load_json(filename):
    """
    Loads data from a Json file into pandas DataFrame.
    """

    with open(filename, "r") as f:
        data = json.load(f)
        df = pd.DataFrame()
        if "catapp" in filename:
            for key in data:
                df = df.append(data[key], ignore_index=True)
        else:
            for key in data:
                df = df.append(data[key]["key_value_pairs"], ignore_index=True)
    print(f"Loaded {len(data)} reactions from file {filename}.")

    return df


def max_reactants_products(df, max_reactants, max_products):
    """
    Finds the maximum number of reactants and products in the reactions.

    Arguments:
      df (DataFrame):       The reaction data.
      max_reactants (int):  The current maximum number of reactants.
      max_products (int):   The current maximum number of products.
    Returns:
      (tuple of ints):      New maximum number of reactants and products,
                             respectively.
    """

    for r in df.reactants:
        rn = len(eval(r))
        if rn > max_reactants:
            max_reactants = rn
    for p in df.products:
        pn = len(eval(p))
        if pn > max_products:
            max_products = pn

    return max_reactants, max_products


def parse_reactants_products(df, db_name, max_reactants, max_products):
    """
    """

    # Drop redundant columns
    drop_names = ["dftCode", "id", "pubId", "username"]
    for name in drop_names:
        try:
            df = df.drop(name, axis=1)
        except KeyError:
            pass

    # Create a list of column names
    col_names = []
    for i in range(max_reactants):
        col_names.append(f"Reactant {i+1}")
    for i in range(max_products):
        col_names.append(f"Product {i+1}")
    col_names.append("Chemical composition")
    col_names.append("Surface composition")
    col_names.append("Facet")
    col_names.append("Adsorption site")
    col_names.append("Coverage")
    col_names.append("Reaction equation")
    col_names.append("DFT functional")
    col_names.append("Reaction energy")
    col_names.append("Activation energy")

    # Fill in the column values
    for i, row in df.iterrows():
        if db_name == "CatApp":
            # Handle CatApp data
            df.at[i, "Reactant 1"] = row.reactant_a
            df.at[i, "Reactant 2"] = row.reactant_b
            df.at[i, "Product 1"] = row.product_ab
            facet = row.reactant_facet
            surface = row.reactant_surface.replace(facet, "")
            facet = facet.split("(")[1].split(")")[0]
            df.at[i, "Facet"] = facet
            df.at[i, "Chemical composition"] = surface
            df.at[i, "Adsorption site"] = row.site
            df.at[i, "Reaction energy"] = row.reaction_energy
            df.at[i, "Activation energy"] = row.activation_energy
        elif db_name == "Catalysis-hub":
            # Handle Catalysis-hub data
            reactants = list(eval(row.reactants))
            for j in range(len(reactants)):
                label = f"Reactant {j+1}"
                df.at[i, label] = reactants[j]
            products = list(eval(row.products))
            for j in range(len(products)):
                label = f"Product {j+1}"
                df.at[i, label] = products[j]
            if row.sites != "None":
                sites = list(eval(row.sites).items())
            if len(sites) > 1:
                site = sites[0]
                try:
                    site = "-".join(site)
                except TypeError:
                    prefix = "-".join(site[1])
                    site = "-".join([site[0], prefix])
            else:
                site = "None"
            df.at[i, "sites"] = site
            if row.coverages != "None":
                coverages = list(eval(row.coverages).items())
                df.at[i, "Coverage"] = row.coverages
        else:
            print(f"Error: No database {db_name} found")

    # Rename columns
    df.rename(columns={
        "chemicalComposition": "Chemical composition",
        "surfaceComposition": "Surface composition",
        "facet": "Facet",
        "sites": "Adsorption site",
        "coverages": "Coverage",
        "Equation": "Reaction equation",
        "dftFunctional": "DFT functional",
        "reactionEnergy": "Reaction energy",
        "activationEnergy": "Activation energy"
    }, inplace=True)

    for name in col_names:
        if name not in df:
            df[name] = "None"
    if db_name == "Catalysis-hub":
        df = df.drop(["reactants", "products"], axis=1)
    df = df[col_names]

    return df


def clean_data(df, drop_duplicates=True, remove_gas=True,
               remove_negatives=True, remove_linears=True, remove_zeros=True):
    """
    """

    # Drop duplicates
    if drop_duplicates:
        df = df.drop_duplicates(ignore_index=True)

    # Unify notational differences
    df = df.replace({"star": "*"}, regex=True)
    df = df.replace(np.nan, "None")
    if remove_gas:
        # Note: In the Catalysis-hub data, gas phase molecules are named with a
        # "gas" prefix, but in the CatApp database they are not. If the flag
        # "remove_gas" is set to True, the Catalysis-hub data will be stripped
        # of this prefix.
        df = df.replace({"gas": ""}, regex=True)

    # Find and drop empty columns
    for col in df.keys():
        if (df[col] == "None").all():
            df = df.drop(col, axis=1)

    # Remove data points with negative activation energy
    if remove_negatives:
        neg_ea = df.iloc[np.where(df["Activation energy"] < 0)]
        for i in neg_ea.index:
            df = df.drop(i)

    # Remove data points where the activation is the same as or the negative of
    # the reaction energy
    if remove_linears:
        linear_ea = df.iloc[np.where(
            df["Reaction energy"] == df["Activation energy"])]
        neg_linear_ea = df.iloc[np.where(
            df["Reaction energy"] == df["Activation energy"])]
        for i in linear_ea.index:
            df = df.drop(i)
        for i in neg_linear_ea.index:
            df = df.drop(i)

    # Remove data points where the activation energy is zero
    if remove_zeros:
        zero_ea = df.iloc[np.where(df["Activation energy"] == 0)]
        for i in zero_ea.index:
            df = df.drop(i)

    return df


def draw_ea_plot(df, root_dir):
    """
    Draws a plot of activation energy as a function of reaction energy,
    and saves it to file.

    Arguments:
      df (pandas DataFrame):  The data from which the plot is drawn.
    """

    # Find the latest image number
    os.chdir(f"{root_dir}/data/images/")
    files = os.listdir()
    # Images are named as "plot_XX.png"
    num = 0
    for file in files:
        inum = int(file[5:7])
        if inum > num:
            num = inum
    num += 1

    # Make plot and save file
    filename = f"{root_dir}/data/images/plot_{num:02}.png"
    plt.plot(df["reactionEnergy"], df["activationEnergy"], "b.")
    plt.xlabel("Reaction energy [eV]")
    plt.ylabel("Activation energy [eV]")
    plt.savefig(filename)


def construct_struct_descriptors(root_dir):
    """
    Construct a list of structural data for each reaction.

    Arguments:
      root_dir (path):    Project root directory
    Returns:
      structures (list):  List of reaction structures, where each item is a
                           list of dictionaries. The dictionaries contain the
                           structures as Ase Atoms objects, and the
                           corresponding energies.
    """

    with open(f"{root_dir}/data/reactions_cathub.json", "r") as f:
        struct_data = json.load(f)

    keys = list(struct_data.keys())
    structures = []
    # Iterate over the reactions
    for key in keys:
        structs = struct_data[key]['structures']
        reaction_structures = []
        # Iterate over the structures available for the current reaction
        for struct in structs:
            struct_dict = {}
            # Save the provided input data into temporary XYZ file
            with open("struct_tmp.xyz", "w") as f:
                f.write(struct["InputFile"])
            # Create an Atoms object from the XYZ file
            atoms = read("struct_tmp.xyz")
            # Save the Atoms object and the corresponding energy to a dictionary
            struct_dict["atoms"] = atoms
            struct_dict["energy"] = struct["energy"]
            # Save the dictionary to a list for the current reaction
            reaction_structures.append(struct_dict)
        # Save the list of structure dictionaries into a global list for all reactions
        structures.append(reaction_structures)

    print(f"Loaded the structures for {len(structures)} reactions.")

    # Build a list of initial-final structure pairs with the corresponding
    # reaction indices
    initial_final_structures, indices = initial_final_structures(structures)

    print(f"Found {len(indices)} initial-final structure pairs.")

    descriptors = []
    targets = []

    for i in range(len(initial_final_structures)):
        descriptors.append(
            build_structural_descriptor(initial_final_structures[i]))
        targets.append(struct_data.iloc[i]["activationEnergy"])

    # Pad the structural descriptors with zeros
    max_len = 0
    for row in descriptors:
        if len(row) > max_len:
            max_len = len(row)
    for i,row in enumerate(descriptors):
        descriptors[i] = list(row) + [0 for i in range(max_len - len(row))]
    descriptors = np.array(descriptors)

    print(f"Built {len(targets)} structural descriptor target pairs.")

    # Save to file
    np.savetxt(f"{root_dir}/data/structural_descriptors.txt", descriptors)
    np.savetxt(f"{root_dir}/data/structural_targets.txt", targets)


def get_initial_final_structures(reaction_structures):
    """
    Get the initial and final structures for the given reaction.
    """

    structure_list = []
    for i,struct in enumerate(reaction_structures):
        # Select only structures that have at least 10 atoms
        # This avoids selecting isolated molecules
        if len(struct["atoms"]) >= 10:
            structure_list.append(struct)
        # Return the first and the last structure only if more than two
        # structures are found
        if len(structure_list) > 1:
            return [structure_list[0], structure_list[-1]]
        else:
            return None


def initial_final_structures(structures):
    """
    Get the initial and final structures for the reactions.

    Parameters:
      structures (list):               List of structural data.
    Returns:
      init_fin_structures (np.array):  The initial and final structures for the
                                        selected reactions.
      indices (list):                  Indices for the selected reactions.
    """

    init_fin_structures = []
    indices = []
    # Get the initial and final structures for each reaction
    for i,struct in enumerate(structures):
        res = get_initial_final_structures(struct)
        if res:
            init_fin_structures.append(res)
            indices.append(i)

    return np.array(init_fin_structures), indices


def build_structural_descriptor(structure):
    """
    """

    # Get the initial and final structures
    initial_structure = Atoms(structure[0]["atoms"])
    final_structure = Atoms(structure[1]["atoms"])

    # Get a list of all the chemical species in the system
    species = []
    species.extend(initial_structure.get_chemical_symbols())
    species.extend(final_structure.get_chemical_symbols())
    species = list(dict.fromkeys(species))

    # Get the system periodicity
    periodic = initial_structure.get_pbc().any()

    # Initialize LMBTR descriptor
    lmbtr = LMBTR(
        species=species,
        k2={
            "geometry": {
                "function": "distance"
            },
            "grid": {
                "min": 0, "max": 5, "n": 50, "sigma": 0.005
            },
            "weighting": {
                "function": "exponential", "scale": 0.5, "cutoff": 1e-3
            },
        },
        k3={
            "geometry": {
                "function": "angle"
            },
            "grid": {
                "min": 0, "max": 180, "n": 50, "sigma": 0.005
            },
            "weighting": {
                "function": "exponential", "scale": 0.5, "cutoff": 1e-3
            }
        },
        periodic=periodic,
        normalization="l2_each",
    )

    # Build LMBTR descriptors from the initial and final structures
    lmbtr_initial = lmbtr.create(initial_structure).flatten()
    lmbtr_final = lmbtr.create(final_structure).flatten()

    # Combine the two LMBTR descriptors into one array
    lmbtr_total = np.concatenate((lmbtr_initial, lmbtr_final), axis=None)

    return lmbtr_total
