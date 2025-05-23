"""
RCSB PDB Connector for Moremi Biokit.

This module provides functions to interact with the RCSB Protein Data Bank (PDB)
APIs (both REST and GraphQL where appropriate). It allows for searching,
downloading PDB files, and fetching metadata associated with PDB entries.

Key functionalities include:
- Downloading PDB files in various formats (PDB, mmCIF, etc).
- Fetching curated summary metadata for PDB entries (for verification or display).
- Performing text-based searches of the PDB.
- Performing sequence-based searches (e.g., using RCSB's sequence search capabilities).
- (Potentially) Performing chemical/ligand-based searches.

This module aims to encapsulate the details of RCSB API endpoints, query construction,
and response parsing, providing a simpler interface for other parts of `moremi_biokit`.
"""

import os
from typing import Optional, Dict, List, Any
from ._utils import make_api_request, APIRequestError
from . import local_pdb_fasta_parser # Import the new local parser module
import requests # For requests.exceptions.JSONDecodeError

# RCSB PDB API Endpoints
RCSB_DATA_API_ENDPOINT = "https://data.rcsb.org/graphql"
RCSB_SEARCH_API_ENDPOINT = "https://search.rcsb.org/rcsbsearch/v2/query"
RSCB_REST_API_ENDPOINT = "https://data.rcsb.org/rest/v1/core/entry"
LOCAL_MOREMI_MICROSERVICE_API_ENDPOINT = "http://localhost:8000/api/v1/sequence"
MOREMI_MICROSERVICE_API_ENDPOINT = "http://167.172.158.230:8000/api/v1/sequence"


def download_pdb_from_rcsb(pdb_id: str, output_dir: str, file_format: str = 'pdb') -> Optional[str]:
    """Download a structure file from RCSB PDB.

    Args:
        pdb_id (str): The 4-character PDB ID (case-insensitive).
        output_dir (str): Directory to save the downloaded file.
        file_format (str, optional): Desired file format ('pdb', 'cif', 'xml', 'fasta', 'pdb.gz', 'cif.gz'). 
                                    Defaults to 'pdb'.

    Returns:
        Optional[str]: The absolute path to the downloaded file on success, None otherwise.
    """
    if not isinstance(pdb_id, str) or len(pdb_id) != 4:
        print(f"Error: Invalid PDB ID '{pdb_id}'. Must be a 4-character string.")
        return None
        
    pdb_id_lower = pdb_id.lower()
    # Construct the download URL
    # Example: https://files.rcsb.org/download/1xyz.pdb
    download_url = f"https://files.rcsb.org/download/{pdb_id_lower}.{file_format}"
    
    # Ensure output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        print(f"Error creating output directory '{output_dir}': {e}")
        return None
        
    output_filename = f"{pdb_id_lower}.{file_format}"
    output_path = os.path.abspath(os.path.join(output_dir, output_filename))

    print(f"Attempting to download {pdb_id} ({file_format}) from RCSB to {output_path}...")
    
    try:
        # Use the utility function for the request
        response = make_api_request(
            url=download_url,
            method="GET",
            stream=True,
            timeout=120, # Can adjust timeout here or rely on default from _utils
            retries=2    # Can adjust retries here or rely on default from _utils
        )

        # Write the content to the file
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                
        print(f"Successfully downloaded {output_filename}.")
        return output_path

    except APIRequestError as api_err:
        # The make_api_request function will print warnings during retries.
        # Here, we print a final error message.
        if api_err.status_code == 404:
            print(f"Error: PDB ID '{pdb_id}' in format '{file_format}' not found on RCSB PDB (URL: {api_err.url}). Final error: {api_err}")
        else:
            print(f"Error downloading from RCSB after retries: {api_err}")
        return None
    except IOError as io_err:
        print(f"File error saving downloaded PDB to {output_path}: {io_err}")
        # Clean up partially downloaded file if it exists
        if os.path.exists(output_path):
            try:
                os.remove(output_path)
            except OSError:
                pass # Ignore error during cleanup
        return None
    except Exception as e:
        print(f"An unexpected error occurred during RCSB download: {e}")
        return None

# (GraphQL query constants and endpoint would be defined here)
#NOTE: An agent would call these tools and use compare these fields (obtained for a PDB ID it selected) against the requirements of its current task.

# TODO: LOOK THROUGH GRAPH QUERY AND MAKE SOME MODIFICATIONS, RAISES ERROR CATCHING AND SO ON
def fetch_pdb_verification_summary(pdb_id: str) -> Optional[Dict[str, Any]]:
    """
    Fetches a curated set of metadata for a given PDB ID from RCSB PDB
    using GraphQL, intended for agent self-verification. If GraphQL fails
    or returns incomplete data, it attempts a fallback to the RCSB REST API
    which might provide a less detailed summary.

    This function aims to retrieve a consistent, structured summary containing
    key information that an autonomous agent can use to confirm if a PDB entry
    matches its intended target criteria before proceeding with downloading
    and using the full PDB file.

    The summary includes:
    - PDB ID (for confirmation)
    - Entry title
    - Keywords
    - Experimental method
    - Resolution (if applicable, e.g., for X-ray or EM)
    - Deposition and release dates
    - Primary citation title, DOI, PubMed ID, and Journal Abbreviation
    - For each polymer entity:
        - Entity ID
        - Description (e.g., protein name)
        - Source organism(s) scientific name
        - Molecule type (e.g., Protein, DNA, RNA)
        - Associated UniProt ID(s)
        - Associated Gene Name(s)
    - Information about non-polymer entities (ligands), if readily available
      in a summary form (e.g., list of ligand IDs).

    Args:
        pdb_id (str): The 4-character PDB identifier (e.g., "4R19", "4HHB").

    Returns:
        Optional[Dict[str, Any]]: A dictionary containing the curated verification
        metadata if the PDB ID is found and data is successfully retrieved.
        The structure of the dictionary aims to be consistent.
        Returns None if the PDB ID is not found, an API error occurs,
        or the response cannot be parsed as expected even after fallback.

    Raises:
        - Potentially requests.exceptions.RequestException for network issues,
          though this function aims to catch them and return None.
          Consider logging errors internally.

    Example of expected output structure (conceptual):
    {
        "pdb_id": "4R19",
        "title": "Crystal Structure of 3D7 strain Plasmodium falciparum AMA1",
        "keywords": ["CELL INVASION"],
        "experimental_method": "X-RAY DIFFRACTION",
        "resolution": 1.8,
        "deposition_date": "2014-08-04T00:00:00Z",
        "initial_release_date": "2015-01-14T00:00:00Z",
        "primary_citation_title": "Structure and Dynamics of Apical Membrane Antigen 1...",
        "primary_citation_doi": "10.1021/bi5012089",
        "primary_citation_pubmed_id": "25584679",
        "primary_citation_journal_abbrev": "Biochemistry",
        "polymer_entities": [
            {
                "entity_id": "1",
                "description": "Apical membrane antigen 1, AMA1",
                "source_organisms": ["Plasmodium falciparum 3D7"],
                "molecule_type": "Protein",
                "uniprot_ids": ["Q7KQK5"],
                "gene_names": ["AMA1", "PF11_0344"]
            }
            // ... more polymer entities if present ...
        ],
        "ligands": ["NAG", "SO4"] // Example: list of ligand IDs
    }


    """
    pdb_id_upper = pdb_id.upper()
    graphql_summary = None # Will store the summary if GraphQL is successful

    # Attempt 1: GraphQL
    graphql_query = f"""
    {{
      entry(entry_id: "{pdb_id_upper}") {{
        rcsb_id # PDB ID
        struct {{
          title # Entry title
        }}
        struct_keywords {{
          pdbx_keywords # Keywords
        }}
        exptl {{
          method # Experimental method
        }}
        rcsb_entry_info {{
          resolution_combined # Resolution
          deposition_date # Deposition date
          initial_release_date # Release date
        }}
        citation {{
          title # Primary citation title
          rcsb_journal_abbrev
          pdbx_database_id_DOI # Primary citation DOI
          pdbx_database_id_PubMed
        }}
        polymer_entities {{
          rcsb_polymer_entity_container_identifiers {{
            entity_id # Entity ID
          }}
          rcsb_entity_source_organism {{
            ncbi_scientific_name # Source organism scientific name
          }}
          rcsb_polymer_entity {{
            pdbx_description # Description (e.g., protein name)
          }}
          entity_poly {{
            rcsb_entity_polymer_type # Molecule type (Protein, DNA, RNA)
          }}
          rcsb_polymer_entity_align {{
            aligned_regions {{
              ref_db_accession 
              ref_db_name
            }}
          }}
          rcsb_polymer_entity_annotation {{
            name 
            type 
            provenance_source
            annotation_lineage {{
                id
                name
            }}
          }}
        }}
        nonpolymer_entities {{
          rcsb_nonpolymer_entity_container_identifiers {{
            nonpolymer_comp_id # Ligand ID
          }}
        }}
      }}
    }}
    """
    try:
        print(f"Attempting GraphQL query for {pdb_id_upper}...")
        response_json = make_api_request(
            url=RCSB_DATA_API_ENDPOINT,
            method="POST",
            json_payload={"query": graphql_query},
            timeout=30
        )

        if response_json:
            # Explicit check for "data: null" with "errors"
            if response_json.get("data") is None and response_json.get("errors"):
                print(f"GraphQL for {pdb_id_upper} returned 'data: null' with errors: {response_json.get('errors')}. Will attempt REST fallback.")
                # graphql_summary remains None, leading to fallback
            elif isinstance(response_json.get("data"), dict) and response_json["data"].get("entry"):
                entry_data_graphql = response_json["data"]["entry"]
                print(f"Successfully fetched data for {pdb_id_upper} using GraphQL.")
                
                # Start parsing, this block can raise parsing errors
                resolution_list = entry_data_graphql.get("rcsb_entry_info", {}).get("resolution_combined")
                resolution_value = resolution_list[0] if isinstance(resolution_list, list) and len(resolution_list) > 0 else None
                citation_data = entry_data_graphql.get("citation", [{}])[0] if entry_data_graphql.get("citation") else {}

                parsed_summary = {
                    "pdb_id": entry_data_graphql.get("rcsb_id"),
                    "title": entry_data_graphql.get("struct", {}).get("title"),
                    "keywords": entry_data_graphql.get("struct_keywords", {}).get("pdbx_keywords", "").split(',') if entry_data_graphql.get("struct_keywords", {}).get("pdbx_keywords") else [],
                    "experimental_method": entry_data_graphql.get("exptl", [{}])[0].get("method") if entry_data_graphql.get("exptl") else None,
                    "resolution": resolution_value,
                    "deposition_date": entry_data_graphql.get("rcsb_entry_info", {}).get("deposition_date"),
                    "initial_release_date": entry_data_graphql.get("rcsb_entry_info", {}).get("initial_release_date"),
                    "primary_citation_title": citation_data.get("title"),
                    "primary_citation_doi": citation_data.get("pdbx_database_id_DOI"),
                    "primary_citation_pubmed_id": citation_data.get("pdbx_database_id_PubMed"),
                    "primary_citation_journal_abbrev": citation_data.get("rcsb_journal_abbrev"),
                    "polymer_entities": [],
                    "ligands": []
                }

                if entry_data_graphql.get("polymer_entities"):
                    for entity in entry_data_graphql["polymer_entities"]:
                        uniprot_ids_set = set()
                        gene_names_set = set()
                        if entity.get("rcsb_polymer_entity_align"):
                            for align_group in entity["rcsb_polymer_entity_align"]:
                                if align_group and align_group.get("aligned_regions"):
                                    for region in align_group["aligned_regions"]:
                                        if region and region.get("ref_db_name") == "UniProt" and region.get("ref_db_accession"):
                                            uniprot_ids_set.add(region["ref_db_accession"])
                        if entity.get("rcsb_polymer_entity_annotation"):
                            for annotation in entity["rcsb_polymer_entity_annotation"]:
                                if annotation:
                                    if annotation.get("type") == "Gene Name" and annotation.get("name"):
                                        gene_names_set.add(annotation["name"].strip())
                                    elif annotation.get("annotation_lineage"):
                                        for lineage_item in annotation["annotation_lineage"]:
                                            if lineage_item and lineage_item.get("name") and annotation.get("type") == "Gene Name":
                                                gene_names_set.add(lineage_item["name"].strip())
                                            elif lineage_item and lineage_item.get("name") and \
                                                 lineage_item["name"].isupper() and 2 <= len(lineage_item["name"]) <= 10 and \
                                                 'gene' in annotation.get("type", "").lower():
                                                 gene_names_set.add(lineage_item["name"].strip())
                        parsed_summary["polymer_entities"].append({
                            "entity_id": entity.get("rcsb_polymer_entity_container_identifiers", {}).get("entity_id"),
                            "description": entity.get("rcsb_polymer_entity", {}).get("pdbx_description"),
                            "source_organisms": [org.get("ncbi_scientific_name") for org in entity.get("rcsb_entity_source_organism", []) if org and org.get("ncbi_scientific_name")],
                            "molecule_type": entity.get("entity_poly", {}).get("rcsb_entity_polymer_type"),
                            "uniprot_ids": sorted(list(uniprot_ids_set)),
                            "gene_names": sorted(list(gene_names_set))
                        })
                if entry_data_graphql.get("nonpolymer_entities"):
                    for np_entity in entry_data_graphql["nonpolymer_entities"]:
                        if np_entity:
                            ligand_id = np_entity.get("rcsb_nonpolymer_entity_container_identifiers", {}).get("nonpolymer_comp_id")
                            if ligand_id:
                                parsed_summary["ligands"].append(ligand_id)
                    parsed_summary["ligands"] = sorted(list(set(parsed_summary["ligands"])))
                graphql_summary = parsed_summary # Assign only after successful parsing

            elif response_json.get("errors"): # Case: No "data:null", but no "entry", yet "errors" exist
                print(f"GraphQL for {pdb_id_upper} returned no entry data but has errors: {response_json.get('errors')}. Will attempt REST fallback.")
                # graphql_summary remains None
            else: # Case: No "entry", no "data:null with errors", no other "errors" field. Could be empty data or unexpected.
                print(f"GraphQL for {pdb_id_upper} returned no entry data or an unexpected structure. Will attempt REST fallback. Response: {response_json}")
                # graphql_summary remains None
        else: # No response_json at all (make_api_request returned None)
            print(f"GraphQL for {pdb_id_upper} returned no response object. Will attempt REST fallback.")
            # graphql_summary remains None

    except APIRequestError as e_graphql:
        print(f"GraphQL API Error for {pdb_id_upper}: {e_graphql}. Will attempt REST fallback.")
        # graphql_summary remains None
    except (KeyError, IndexError, TypeError) as e_parse_graphql:
        import traceback
        print(f"Error parsing GraphQL response for {pdb_id_upper} (data was received but malformed for parser): {e_parse_graphql}\\n{traceback.format_exc()}")
        # Critical decision: if GraphQL data is received but unparseable, do not fallback.
        # This might indicate an issue with our query or an unexpected schema change.
        return None 
    except Exception as e_general_graphql:
        import traceback
        print(f"An unexpected error occurred during GraphQL fetch for {pdb_id_upper}: {e_general_graphql}\\n{traceback.format_exc()}. Will attempt REST fallback.")
        # graphql_summary remains None
    
    if graphql_summary: # Check if GraphQL attempt (and parsing) was successful
        return graphql_summary

    # Attempt 2: REST API Fallback (if graphql_summary is still None)
    print(f"Attempting REST API fallback for {pdb_id_upper}...")
    # rest_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id_upper}" # User updated this
    rest_url = f"{RSCB_REST_API_ENDPOINT}/{pdb_id_upper}"
    try:
        rest_response_json = make_api_request(
            url=rest_url,
            method="GET",
            timeout=30
        )

        if not rest_response_json:
            print(f"REST API fallback for {pdb_id_upper} also failed to return data.")
            return None

        print(f"Successfully fetched data for {pdb_id_upper} using REST API fallback.")
        
        rest_resolution_list = rest_response_json.get("rcsb_entry_info", {}).get("resolution_combined")
        rest_resolution_value = rest_resolution_list[0] if isinstance(rest_resolution_list, list) and len(rest_resolution_list) > 0 else None
        rest_citation_data = rest_response_json.get("citation", [{}])[0] if rest_response_json.get("citation") else {}

        summary_rest = {
            "pdb_id": rest_response_json.get("rcsb_id"),
            "title": rest_response_json.get("struct", {}).get("title"),
            "keywords": rest_response_json.get("struct_keywords", {}).get("pdbx_keywords", "").split(',') if rest_response_json.get("struct_keywords", {}).get("pdbx_keywords") else [],
            "experimental_method": rest_response_json.get("exptl", [{}])[0].get("method") if rest_response_json.get("exptl") else None,
            "resolution": rest_resolution_value,
            "deposition_date": rest_response_json.get("rcsb_entry_info", {}).get("deposition_date"),
            "initial_release_date": rest_response_json.get("rcsb_entry_info", {}).get("initial_release_date"),
            "primary_citation_title": rest_citation_data.get("title"),
            "primary_citation_doi": rest_citation_data.get("pdbx_database_id_DOI"),
            "primary_citation_pubmed_id": rest_citation_data.get("pdbx_database_id_PubMed"),
            "primary_citation_journal_abbrev": rest_citation_data.get("rcsb_journal_abbrev"),
            "polymer_entities": [],
            "ligands": [],
            "entry_source_organisms": [] # Specific to REST fallback for now
        }

        polymer_entity_ids = rest_response_json.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        for entity_id in polymer_entity_ids:
            summary_rest["polymer_entities"].append({
                "entity_id": entity_id,
                "description": "Details require specific polymer_entity REST call",
                "source_organisms": [], 
                "molecule_type": "N/A via entry REST call",
                "uniprot_ids": [],
                "gene_names": []
            })
        
        source_orgs_rest = []
        if rest_response_json.get("rcsb_entity_source_organism"):
            for org in rest_response_json["rcsb_entity_source_organism"]:
                if org and org.get("ncbi_scientific_name"):
                    source_orgs_rest.append(org["ncbi_scientific_name"])
        if source_orgs_rest:
             summary_rest["entry_source_organisms"] = sorted(list(set(source_orgs_rest)))

        nonpolymer_entity_ids = rest_response_json.get("rcsb_entry_container_identifiers", {}).get("nonpolymer_entity_ids", [])
        if nonpolymer_entity_ids:
             summary_rest["ligands"] = sorted(list(set(nonpolymer_entity_ids)))
        
        print(f"Warning: Used REST API fallback for {pdb_id_upper}. Summary may be less detailed.")
        return summary_rest

    except APIRequestError as rest_e:
        print(f"REST API fallback for {pdb_id_upper} also failed: {rest_e}")
        return None
    except (KeyError, IndexError, TypeError) as parse_e_rest:
        import traceback
        print(f"Error parsing REST API fallback response for {pdb_id_upper}: {parse_e_rest}\n{traceback.format_exc()}")
        return None
    except Exception as final_e:
        import traceback
        print(f"An unexpected error occurred during REST API fallback for {pdb_id_upper}: {final_e}\n{traceback.format_exc()}")
        return None


def search_rcsb_by_text(
    query_string: str,
    max_results: int = 25,
    result_type: str = "entry"
) -> Optional[List[Dict[str, Any]]]:
    """
    Performs a text-based search of the RCSB PDB using their search API.

    This function allows querying the PDB with keywords, author names,
    PDB IDs, or more complex query strings supported by the RCSB search syntax.
    It returns a list of matching PDB entries, typically including their IDs
    and a brief summary or relevance score.

    Args:
        query_string (str): The search query string.
            (e.g., "hemoglobin AND human AND resolution:<2.0")
        max_results (int, optional): The maximum number of search results
            to return. Defaults to 25.
        result_type (str, optional): The type of result to search for,
            such as "entry", "polymer_entity", "non_polymer_entity", "assembly". 
            Affects the `results_content_type` and `return_type` in the query.
            Defaults to "entry".

    Returns:
        Optional[List[Dict[str, Any]]]: A list of dictionaries, where each
        dictionary represents a search hit (e.g., {"pdb_id": "XXXX", "score": Y.YY}),
        or None if the search fails or returns no results. The exact structure
        of the hit dictionaries depends on the RCSB Search API response.
    """
    # The RCSB Search API uses a specific JSON structure for queries.
    # We define a terminal query for text search.
    # `results_content_type` can be ["experimental"], ["computational"] or combined.
    # For simplicity, we'll let the API decide based on the query if not specified,
    # or use a broad default. Let's stick to what `result_type` implies for `return_type`.

    # Mapping result_type to API's expected return_type and content_type if needed.
    # Common return types: "entry", "polymer_entity", "non_polymer_entity", "assembly", etc.
    # For simplicity, this implementation will primarily use "entry" for broad searches returning PDB IDs.
    # If a different `result_type` is passed that matches a valid Search API return type, it will be used.
    
    # Validate result_type for common cases, can be expanded
    valid_return_types = ["entry", "polymer_entity", "non_polymer_entity", "assembly", "moldefinition"]
    if result_type not in valid_return_types:
        print(f"Warning: result_type '{result_type}' may not be directly supported by Search API as a return_type. Defaulting to 'entry'.")
        api_return_type = "entry"
    else:
        api_return_type = result_type

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "value": query_string
            }
        },
        "request_options": {
            "pager": {
                "start": 0,
                "rows": max_results
            },
            "scoring_strategy": "combined", # Options: "combined", "text", "sequence", etc.
            # "results_content_type": ["experimental", "computational"] # Broadest search by default
        },
        "return_type": api_return_type
    }

    try:
        response_json = make_api_request(
            url=RCSB_SEARCH_API_ENDPOINT,
            method="POST",
            json_payload=query,
            timeout=60  # Searches can sometimes take longer
        )

        if not response_json:
            print(f"Error: No response from server for text query '{query_string}'.")
            return None

        if response_json.get("errors"):
            print(f"Search API errors for query '{query_string}': {response_json.get('errors')}")
            return None
        
        if "result_set" not in response_json:
            # It's possible to get a 0 results count, which is not an error but an empty set
            if response_json.get("total_count", 0) == 0:
                return [] # No results found, return empty list
            print(f"Error: Malformed response for text query '{query_string}'. 'result_set' not found.")
            print(f"Full response: {response_json}")
            return None

        results = []
        for hit in response_json.get("result_set", []):
            # The identifier is the primary ID of the returned result type (e.g., PDB ID for "entry")
            result_item = {
                "identifier": hit.get("identifier"), 
                "score": hit.get("score")
            }
            # Add more details from the hit if they are consistently available and useful
            # For example, for "entry" type, one might try to extract basic info if present
            if api_return_type == "entry" and hit.get("services"):
                for service in hit.get("services", []):
                    if service.get("service_name") == "text" and service.get("nodes"):
                        entry_node = service.get("nodes", [{}])[0]
                        if entry_node.get("match_context"): # Match context might have highlights
                             pass # Further processing could be added here
            results.append(result_item)
        
        return results

    except APIRequestError as e:
        print(f"API Error during text search for '{query_string}': {e}")
        return None
    except (KeyError, TypeError) as e: # Removed IndexError as get() handles it
        import traceback
        print(f"Error parsing text search response for '{query_string}': {e}\n{traceback.format_exc()}")
        return None
    except Exception as e:
        import traceback
        print(f"An unexpected error occurred during text search for '{query_string}': {e}\n{traceback.format_exc()}")
        return None


def search_rcsb_by_sequence(
    sequence: str,
    identity_cutoff: float = 0.9,
    evalue_cutoff: Optional[float] = 10.0, # Explicitly setting a default as API might require it
    search_target: str = "pdb_protein_sequence",
    max_results: int = 10
) -> Optional[List[Dict[str, Any]]]:
    """
    Performs a sequence similarity search against sequences in the RCSB PDB.

    This function typically uses the RCSB PDB's sequence search tool (which
    might wrap MMseqs2 or similar) to find PDB entries containing sequences
    similar to the input query sequence.

    Args:
        sequence (str): The amino acid or nucleotide query sequence.
        identity_cutoff (float, optional): Minimum sequence identity for a hit.
            Defaults to 0.9 (90%).
        evalue_cutoff (Optional[float], optional): Maximum E-value for a hit.
            If None, a default of 10.0 is used. To indicate no E-value cutoff (if supported by API),
            a very large number might be needed, but typically an E-value is used.
            Defaults to 10.0.
        search_target (str, optional): The database target for the sequence search
            (e.g., "pdb_protein_sequence", "pdb_rna_sequence", "pdb_dna_sequence").
            Defaults to "pdb_protein_sequence".
        max_results (int, optional): The maximum number of hits to return.
            Defaults to 10.

    Returns:
        Optional[List[Dict[str, Any]]]: A list of dictionaries, where each
        dictionary represents a PDB entity hit (e.g., {"pdb_id": "XXXX",
        "entity_id": "Y", "chain_id": "A" (if identifiable from result),
        "identity": 0.95, "evalue": 1e-50}),
        or None if the search fails or no significant hits are found.
    """
    target_map = {
        "pdb_protein_sequence": "protein",
        "pdb_rna_sequence": "rna",
        "pdb_dna_sequence": "dna"
    }
    if search_target not in target_map:
        print(f"Error: Invalid search_target '{search_target}'. Supported targets are: {list(target_map.keys())}")
        return None
    
    rcsb_sequence_type = target_map[search_target]

    # Ensure sequence is uppercase as the API might be case-sensitive or prefer it.
    sequence_upper = sequence.upper()

    # The Search API expects evalue_cutoff. If user passes None, we use the function's default.
    current_evalue_cutoff = evalue_cutoff if evalue_cutoff is not None else 10.0

    query = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "evalue_cutoff": current_evalue_cutoff,
                "identity_cutoff": identity_cutoff,
                "sequence_type": rcsb_sequence_type,
                "value": sequence_upper
            }
        },
        "request_options": {
            "pager": {
                "start": 0,
                "rows": max_results
            },
            "scoring_strategy": "sequence" # Specific for sequence searches
        },
        "return_type": "polymer_entity"
    }

    try:
        response_json = make_api_request(
            url=RCSB_SEARCH_API_ENDPOINT,
            method="POST",
            json_payload=query,
            timeout=180 # Sequence searches can be computationally intensive
        )

        if not response_json:
            print(f"Error: No response from server for sequence query.")
            return None
        
        if response_json.get("errors"):
            print(f"Search API errors for sequence query: {response_json.get('errors')}")
            return None
        
        if "result_set" not in response_json:
            if response_json.get("total_count", 0) == 0:
                return [] # No results found
            print(f"Error: Malformed response for sequence query. 'result_set' not found.")
            print(f"Full response: {response_json}")
            return None

        results = []
        for hit in response_json.get("result_set", []):
            identifier = hit.get("identifier", "_") # e.g., "4HHB_1"
            pdb_id, entity_id_str = identifier.split('_') if '_' in identifier else (identifier, None)
            
            # Extracting alignment details from the service node in the hit.
            # The structure is typically: hit -> services -> nodes -> match_context
            service_data = hit.get("services", [{}])[0] # Assuming one service for sequence search result
            node_data = service_data.get("nodes", [{}])[0] if service_data.get("nodes") else {}
            match_context = node_data.get("match_context", [{}])[0] if node_data.get("match_context") else {}

            result_item = {
                "pdb_id": pdb_id,
                "entity_id": entity_id_str, # This is a string as per PDB conventions
                "identity": match_context.get("alignment_identity"),
                "evalue": match_context.get("evalue"),
                "score": hit.get("score") # Overall relevance score from the search
            }
            results.append(result_item)
        
        return results

    except APIRequestError as e:
        print(f"API Error during sequence search: {e}")
        return None
    except (KeyError, TypeError) as e: # Removed IndexError
        import traceback
        print(f"Error parsing sequence search response: {e}\n{traceback.format_exc()}")
        return None
    except Exception as e:
        import traceback
        print(f"An unexpected error occurred during sequence search: {e}\n{traceback.format_exc()}")
        return None

def fetch_sequence_details_by_pdb_chain_id(pdb_chain_id: str) -> Optional[Dict[str, Any]]:
    """Fetches sequence and details for a specific PDB ID and chain ID from the local microservice.

    Args:
        pdb_chain_id (str): The PDB ID and chain ID (e.g., "101m_A", case-insensitive for input parts).

    Returns:
        Optional[Dict[str, Any]]: A dictionary containing sequence details:
            {
                "id": "pdbid_CHAIN", (e.g., "101m_A", normalized)
                "pdb_id": "101m", (lowercase)
                "chain_id": "A", (uppercase)
                "sequence": "SEQUENCE...",
                "length": sequence_length,
                "description": "Description from API or a default."
            }
        Returns None if the sequence cannot be found or an error occurs.
    """
    if not isinstance(pdb_chain_id, str) or '_' not in pdb_chain_id:
        # Using print for now, consider proper logging for a library
        print(f"Error (fetch_sequence_details_by_pdb_chain_id): Invalid pdb_chain_id format: '{pdb_chain_id}'. Expected 'PDBID_CHAIN'.")
        return None

    try:
        pdb_id_part, chain_id_part = pdb_chain_id.split('_', 1)
        if not pdb_id_part or not chain_id_part: # Ensure both parts are non-empty after split
             raise ValueError("PDB ID or Chain ID part is empty after split.")
    except ValueError as e:
        print(f"Error (fetch_sequence_details_by_pdb_chain_id): Could not parse pdb_chain_id '{pdb_chain_id}': {e}")
        return None

    normalized_pdb_id = pdb_id_part.lower()
    normalized_chain_id = chain_id_part.upper()
    # This is the format the API endpoint expects and also the "id" for the output dict
    api_target_pdb_chain_id = f"{normalized_pdb_id}_{normalized_chain_id}"

    # Microservice endpoint URL
    # For local dev: http://localhost:8000/api/v1/sequence/{pdb_chain_id}
    # api_url = f"{LOCAL_MOREMI_MICROSERVICE_API_ENDPOINT}/{api_target_pdb_chain_id}"
    api_url = f"{MOREMI_MICROSERVICE_API_ENDPOINT}/{api_target_pdb_chain_id}"
    
    
    # Consider using proper logging in a library context
    print(f"Info (fetch_sequence_details_by_pdb_chain_id): Fetching from microservice: {api_url}")

    try:
        response = make_api_request(url=api_url, method="GET")
        # make_api_request handles retries and raises APIRequestError for HTTP errors.
        # If it returns, it means the HTTP request itself was successful (e.g. 200 OK).
        
        api_data = response.json() # Can raise requests.exceptions.JSONDecodeError

        if not isinstance(api_data, dict):
            print(f"Error (fetch_sequence_details_by_pdb_chain_id): API response for '{api_target_pdb_chain_id}' is not a valid JSON object.")
            return None

        sequence = api_data.get("sequence")
        if not sequence or not isinstance(sequence, str):
            print(f"Error (fetch_sequence_details_by_pdb_chain_id): 'sequence' not found or is not a string in API response for '{api_target_pdb_chain_id}'.")
            return None

        # Extract details, ensuring correct casing and types for the output dict
        # Use API provided pdb_id and chain_id if available, otherwise stick to normalized input parts.
        # The API sample returns "pdb_id" and "chain_id" correctly cased.
        out_pdb_id = str(api_data.get("pdb_id", normalized_pdb_id)).lower()
        out_chain_id = str(api_data.get("chain_id", normalized_chain_id)).upper()
        
        out_length = api_data.get("length")
        if not isinstance(out_length, int):
            print(f"Warning (fetch_sequence_details_by_pdb_chain_id): 'length' from API for '{api_target_pdb_chain_id}' is missing or not an integer. Calculating from sequence.")
            out_length = len(sequence)
        
        out_description = api_data.get("description", f"Sequence for {api_target_pdb_chain_id}")

        return {
            "id": api_target_pdb_chain_id, # e.g., "101m_a" -> "101m_A" as per api_target_pdb_chain_id normalization
            "pdb_id": out_pdb_id,
            "chain_id": out_chain_id,
            "sequence": sequence,
            "length": out_length,
            "description": str(out_description) # Ensure description is string
        }

    except APIRequestError as e:
        # Logged by make_api_request during retries, this is the final failure.
        print(f"Error (fetch_sequence_details_by_pdb_chain_id): API request failed for '{api_target_pdb_chain_id}': {e}")
        return None
    except requests.exceptions.JSONDecodeError as e:
        print(f"Error (fetch_sequence_details_by_pdb_chain_id): Failed to decode JSON response from '{api_url}': {e}")
        return None
    except Exception as e: # Catch any other unexpected errors during processing
        print(f"Error (fetch_sequence_details_by_pdb_chain_id): An unexpected error occurred for '{api_target_pdb_chain_id}': {e}")
        return None

# TODO: CHANGE THIS FUNCTIONALITY TO USE NEW DOCKERIZED MOREMI-SERVICE API
def get_all_pdb_sequences_details(pdb_id: str) -> List[Dict[str, Any]]:
    """Fetches all sequence details for a given PDB ID.

    It first attempts to retrieve sequences from the local PDB FASTA database.
    If not found locally (or if the local DB is unavailable), it falls back to
    fetching the FASTA file from the RCSB online service and parsing all sequences.

    Args:
        pdb_id (str): The 4-character PDB ID (case-insensitive).

    Returns:
        List[Dict[str, Any]]: A list of dictionaries, where each dictionary
        contains sequence details for a chain, structured as in
        `get_pdb_chain_sequence_details`. Returns an empty list if no
        sequences can be found or an error occurs.
    """
    pass

