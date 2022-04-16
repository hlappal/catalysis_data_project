# Imports
import numpy as np
import requests
import json
import io
import os
import ase.io
from ase.io import read

# Define the Catalysis-hub API path and the project root directory
GRAPHQL = "http://api.catalysis-hub.org/graphql"
# ROOT_DIR = os.path.join(os.getcwd(), os.pardir)
ROOT_DIR = os.getcwd()

# Define the keyvalues used in the query
KEY_VALUES = [
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

def query_reactions(endcursor):
    """
    The function performs a batch query on the database. The database is
    queried in batches of 50 reactions, and all reactions with an
    activation enregy under 100 eV are selected.
    
    Parameters:
      endcursor:  Cursor to indicate the batch on the query.
    Returns:
      data:       The acquired data on each batch.
    """
    
    # Define the query string
    query_string = "{"
    query_string += f'reactions(first: 50, after: "{endcursor}"'
    query_string += ', activationEnergy: 100, op: "<"'
    query_string += """) {
  totalCount
  pageInfo {
    endCursor
  }
  edges {
    node {"""
  
    # Add the keywords into the query string
    for key_value in KEY_VALUES:
        query_string += str("\n" + " "*6 + key_value)
        
    query_string += """
      systems {
        id
        Trajdata
        energy
        InputFile(format: "xyz")
        keyValuePairs
      }
    }
  }
}}"""
    
    data = requests.post(GRAPHQL, {"query": query_string})
    try:
        # Read the acquired data into a dictionary
        data = data.json()["data"]
    except Exception as e:
        # Handle exceptions in a general manner
        print(e)
        print("Error: Something went wrong. Please check your query string.")
    
    return data

def parse_reaction(reaction):
    """
    The function parses a single reaction. All the missing keyvalues are
    labeled as 'None', and the structural data is saved separately.
    
    Parameters:
      reaction:       A dictionary containing the data for a single reaction.
    Returns:
      reaction_dict:  A parsed data dictionary.
    """
    
    reaction_dict = {}
    key_value_pairs = {}
    
    # Go through the keyvalues
    for key_value in KEY_VALUES:
        try:
            key_value_pairs[key_value] = reaction[key_value]
        except ValueError:
            key_value_pairs[key_value] = "None"
    if key_value_pairs["coverages"] is None:
        key_value_pairs["coverages"] = "None"
    if key_value_pairs["sites"] is None:
        key_value_pairs["sites"] = "None"
    reaction_dict["key_value_pairs"] = key_value_pairs
    
    # Go through the structural data
    structures = []
    for structure in reaction["systems"]:
        struct = {}
        struct["energy"] = structure["energy"]
        struct["InputFile"] = structure["InputFile"]
        struct["keyValuePairs"] = structure["keyValuePairs"]
        structures.append(struct)
    reaction_dict["structures"] = structures
    
    return reaction_dict

def main():
    # Run queries and save results to file
    reaction_list = {}
    N_fetched = 0
    endcursor = ""
    n = 0
    totalcount = 100000

    while n * 50 + 1 < totalcount:
        data = query_reactions(endcursor)
        for reaction in data["reactions"]["edges"]:
            reaction = reaction["node"]
            reaction_list[reaction["id"]] = parse_reaction(reaction)
        endcursor = data["reactions"]["pageInfo"]["endCursor"]
        totalcount = data["reactions"]["totalCount"]
        N_fetched += totalcount
        count = 50 * (n + 1)
        if count >= totalcount:
            count = totalcount
        print(f"Fetched reactions {50*n+1}-{count}/{totalcount}")
        n += 1
    print("Done!")

    with open (f"{ROOT_DIR}/data/reactions_cathub.json", "w") as outfile:
        json.dump(reaction_list, outfile)

if __name__ == "__main__":
    main()