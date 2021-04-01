#!/usr/bin/env python3

"""
"""

import json
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt


def preprocess(root_dir):
    """
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

    return df_raw


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
    """

    with open(f"{root_dir}/data/reactions_cathub.json", "r") as f:
        struct_data = json.load(f)
