import sys
import os

from charge.servers.server_utils import update_mcp_network, get_hostname
from charge.servers.molecular_property_utils import chemprop_preds_server
from charge.servers import SMILES_utils
from mcp.server.fastmcp import FastMCP
from loguru import logger
from rdkit import Chem
import json

JSON_FILE_PATH = "known_molecules.json"

mcp = FastMCP(
    "SMILES Diagnosis and retrieval MCP Server",
)


@mcp.tool()
def is_already_known(smiles: str) -> bool:
    """
    Check if a SMILES string provided is already known. Only provide
    valid SMILES strings. Returns True if the SMILES string is valid, and
    already in the database, False otherwise.
    Args:
        smiles (str): The input SMILES string.
    Returns:
        bool: True if the SMILES string is valid and known, False otherwise.

    Raises:
        ValueError: If the SMILES string is invalid.
    """
    if not Chem.MolFromSmiles(smiles):
        raise ValueError("Invalid SMILES string.")

    try:
        canonical_smiles = SMILES_utils.canonicalize_smiles(smiles)

        try:
            with open(JSON_FILE_PATH) as f:
                known_mols = json.load(f)
                known_smiles = [mol["smiles"] for mol in known_mols]

        except FileNotFoundError:
            logger.warning(f"{JSON_FILE_PATH} not found. Creating a new one.")
            known_mols = []

    except Exception as e:
        raise ValueError("Error in canonicalizing SMILES string.") from e

    # Check if the SMILES string is already known (in the database)
    # This is a placeholder for the actual database check
    return canonical_smiles in known_smiles


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run a ChARGe MCP Server")
    parser.add_argument(
        "--port", type=int, default=8124, help="Port to run the server on"
    )
    parser.add_argument(
        "--host", type=str, default=None, help="Host to run the server on"
    )
    args = parser.parse_args()

    port = args.port
    host = args.host

    mcp.tool()(SMILES_utils.canonicalize_smiles)
    mcp.tool()(SMILES_utils.verify_smiles)
    mcp.tool()(chemprop_preds_server)

    if host is None:
        _, host = get_hostname()

    update_mcp_network(mcp, host, port)

    print(chemprop_preds_server("CCO", "density"))
    mcp.run(transport="sse")
