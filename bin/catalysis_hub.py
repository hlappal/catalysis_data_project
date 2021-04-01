#!/usr/bin/env python3

"""
A script for querying the Catalysis-hub (https://catalysis-hub.org) database.
The data comprises heterogeneous catalytic reactions with atomic structures.
"""

import json
import requests


def catalysis_hub(root_dir):
    """
    Reaction dictionaries are saved into a larger dictionary using
    the unique ids provided in the reaction data.
    """

    reaction_list = {}
    N_fetched = 0
    endcursor = ""
    N = 0
    totalcount = 100000

    # Define selected key values for the reactions
    keyvalues = [
        "chemicalComposition",
        "surfaceComposition",
        "facet",
        "sites",
        "coverages",
        "reactants",
        "products",
        "Equation",
        "reactionEnergy",
        "activationEnergy",
        "dftCode",
        "dftFunctional",
        "username",
        "pubId",
        "id",
    ]

    # The database is queried in batches.
    # The variable totalcount is updated from the arbitrary
    # initial value to the actual number of found reactions.
    while N * 50 + 1 < totalcount:
        # Fetch a reaction data batch
        data = query_reactions(endcursor, keyvalues)

        # Record the reaction data into a dictionary
        for reaction in data["reactions"]["edges"]:
            reaction = reaction["node"]
            reaction_list[reaction["id"]] = parse_reaction(reaction, keyvalues)
        endcursor = data["reactions"]["pageInfo"]["endCursor"]
        totalcount = data["reactions"]["totalCount"]
        N_fetched += totalcount
        count = 50 * (N + 1)
        if count >= totalcount:
            count = totalcount
        print(f"Fetched reactions {50*N+1}-{count}/{totalcount}")
        N += 1

    with open(f"{root_dir}/data/reactions_cathub.json", "w") as outfile:
        json.dump(reaction_list, outfile)

    print(f"Saved {totalcount} reactions into file data/reactions_cathub.json")


def query_reactions(endcursor, keyvalues):
    """
    Run a batch query of 50 reactions. Searching the database in batches
    prevents the backend from overloading. The endcursor defines the starting
    position for each batch query. The query string is a structured string
    that contains all the key values that are included in the search.

    Parameters:
      endcursor (string):  Cursor position for the batch
    Returns:
      data (dict):         The reaction data for the current batch
    """

    # Define the query string
    # (for more info, see: http://docs.catalysis-hub.org/en/latest/)
    query_string = "{"
    query_string += 'reactions('
    query_string += 'first: 50,'
    query_string += f'after: "{endcursor}"'
    query_string += ', activationEnergy: 10, op: "<"'
    query_string += ") {"
    query_string += "  totalCount"
    query_string += "  pageInfo {"
    query_string += "    endCursor"
    query_string += "  }"
    query_string += "  edges {"
    query_string += "    node {"

    # Append the key values with correct structure
    for kv in keyvalues:
        query_string += str("\n" + " "*6 + kv)

    query_string += "\n"
    query_string += "      systems {"
    query_string += "        id"
    query_string += "        Trajdata"
    query_string += "        energy"
    query_string += '        InputFile(format: "xyz")'
    query_string += "        keyValuePairs"
    query_string += "      }"
    query_string += "    }"
    query_string += "  }"
    query_string += "}}"

    # Perform the batch query
    root = 'http://api.catalysis-hub.org/graphql'
    data = requests.post(root, {"query": query_string})
    try:
        data = data.json()["data"]
    except Exception as e:
        print(e)
        print("Error: Something went wrong.")

    return data


def parse_reaction(reaction, keyvalues):
    """
    Parse the reaction data so that all missing values are
    assigned a string value "None". The data for each reaction
    is saved in a dictionary.

    Parameters:
      reaction (dict):       Raw data describing a single reaction
    Returns:
      reaction_dict (dict):  Parsed reaction data
    """

    # "reaction_dict" has two keys: "kv_pairs" and "structures"
    reaction_dict = {}

    # Record the reaction key values into a dictionary
    kv_pairs = {}
    for kv in keyvalues:
        try:
            kv_pairs[kv] = reaction[kv]
        except ValueError:
            kv_pairs[kv] = "None"
    # Some null values need to be parsed explicitly
    for kv in keyvalues:
        if not kv_pairs[kv]:
            kv_pairs[kv] = "None"
    reaction_dict["key_value_pairs"] = kv_pairs

    # Record the structural data into a list of dictionaries
    structures = []
    for structure in reaction["systems"]:
        struct = {}
        struct["energy"] = structure["energy"]
        struct["InputFile"] = structure["InputFile"]
        struct["keyValuePairs"] = structure["keyValuePairs"]
        structures.append(struct)
    reaction_dict["structures"] = structures

    return reaction_dict

