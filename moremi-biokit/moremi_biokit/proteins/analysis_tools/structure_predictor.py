"""
Predicts protein structure using the SWISS-MODEL API.

Submits a sequence to the SWISS-MODEL service, polls for completion,
and retrieves the predicted structure in PDB format along with quality metrics.
"""

import requests
import time
import os
import uuid
import random
from io import BytesIO
import gzip
from datetime import datetime
from typing import Optional, Dict, Any, List
import argparse

# Consider moving to a configuration file or environment variables for production
# For now, keeping as a placeholder if token is not directly provided.
DEFAULT_SWISS_MODEL_API_TOKENS: List[str] = [
    "9d4d090c7ab0283e1f5fb3a92a58fbd6d682dd31", 
    "6d79f42476bf08d1d7777f99446cdc92a755bf8a", 
    "d7a2ddd5a82a33c44e5c4634b48bc77d6d4c35c7",
    "4445722219b06e22c6dd5eee401ba6c02cd14398"
]

# SWISS-MODEL API endpoints
SWISS_MODEL_AUTOMODEL_URL = "https://swissmodel.expasy.org/automodel"
SWISS_MODEL_PROJECT_URL = "https://swissmodel.expasy.org/project/{project_id}/models/summary/"

DEFAULT_POLLING_INTERVAL_SECONDS = 15
DEFAULT_MAX_POLLING_ATTEMPTS = 60 # Max 15 minutes for a typical job

def predict_structure(
    sequence: str,
    api_token: Optional[str] = None,
    project_title_prefix: str = "moremi_biokit_prediction",
    output_directory: str = ".",
    output_pdb_filename_prefix: str = "predicted_structure",
    polling_interval: int = DEFAULT_POLLING_INTERVAL_SECONDS,
    max_polling_attempts: int = DEFAULT_MAX_POLLING_ATTEMPTS
) -> Dict[str, Any]:
    """Predict protein structure using the SWISS-MODEL API.

    Submits a protein sequence to SWISS-MODEL, polls for job completion,
    retrieves model details, and saves the PDB coordinates.

    Args:
        sequence (str): The amino acid sequence for structure prediction.
        api_token (Optional[str], optional): SWISS-MODEL API token. 
            If None, a token from a predefined list will be attempted (not recommended for production).
            It is highly recommended to provide your own token.
            Defaults to None.
        project_title_prefix (str, optional): Prefix for the project title on SWISS-MODEL. 
            A UUID will be appended. Defaults to "moremi_biokit_prediction".
        output_directory (str, optional): Directory to save the output PDB file. 
            Defaults to current directory (".").
        output_pdb_filename_prefix (str, optional): Prefix for the PDB filename. 
            A timestamp and .pdb extension will be added. Defaults to "predicted_structure".
        polling_interval (int, optional): Time in seconds between polling attempts for job status. 
            Defaults to DEFAULT_POLLING_INTERVAL_SECONDS.
        max_polling_attempts (int, optional): Maximum number of polling attempts before timeout. 
            Defaults to DEFAULT_MAX_POLLING_ATTEMPTS.

    Returns:
        Dict[str, Any]: A dictionary containing the prediction results or an error message.
            On success, includes keys like "model_details" (with GMQE, QMEAN etc.), 
            "pdb_file_path", and "pdb_content".
            On failure, includes an "error" key with a descriptive message.
    """
    if not isinstance(sequence, str) or not sequence.strip():
        return {"error": "Input sequence cannot be empty."}
    cleaned_sequence = sequence.strip().upper()

    chosen_api_token = api_token
    if not chosen_api_token:
        if not DEFAULT_SWISS_MODEL_API_TOKENS:
             return {"error": "SWISS-MODEL API token is required and no default tokens are available."}
        chosen_api_token = random.choice(DEFAULT_SWISS_MODEL_API_TOKENS)
        # Consider logging a warning here if using default tokens

    headers = {"Authorization": f"Token {chosen_api_token}"}
    project_uuid = str(uuid.uuid4())
    project_title = f"{project_title_prefix}_{project_uuid}"

    payload = {
        "target_sequences": [cleaned_sequence],
        "project_title": project_title
    }

    try:
        # 1. Submit job
        submit_response = requests.post(SWISS_MODEL_AUTOMODEL_URL, headers=headers, json=payload, timeout=30)
        submit_response.raise_for_status() # Raises HTTPError for bad responses (4XX or 5XX)
        
        project_id = submit_response.json().get("project_id")
        if not project_id:
            return {"error": "Failed to get project_id from SWISS-MODEL submission.", "details": submit_response.text}

        # 2. Poll for job completion
        job_status_url = SWISS_MODEL_PROJECT_URL.format(project_id=project_id)
        for attempt in range(max_polling_attempts):
            time.sleep(polling_interval)
            status_response = requests.get(job_status_url, headers=headers, timeout=30)
            status_response.raise_for_status()
            
            status_data = status_response.json()
            current_job_status = status_data.get("status")
            # print(f"Attempt {attempt + 1}/{max_polling_attempts}: Job {project_id} status: {current_job_status}") # Optional: for verbose logging

            if current_job_status == "COMPLETED":
                break
            elif current_job_status == "FAILED":
                job_errors = status_data.get("errors", [])
                error_messages = ", ".join([e.get("message", "Unknown error") for e in job_errors])
                return {"error": f"SWISS-MODEL job {project_id} failed.", "details": error_messages or status_data}
            elif attempt == max_polling_attempts - 1:
                return {"error": f"SWISS-MODEL job {project_id} timed out after {max_polling_attempts * polling_interval} seconds. Last status: {current_job_status}"}
        else: # Loop completed without break (should be caught by timeout)
            return {"error": f"SWISS-MODEL job {project_id} polling finished unexpectedly. Last status: {current_job_status}"}

        # 3. Retrieve model data if COMPLETED
        # print(f"status_data: {status_data}")
        models = status_data.get("models", [])
        # print(f"models: {models}")
        if not models:
            return {"error": f"SWISS-MODEL job {project_id} completed but no models found.", "details": status_data}
        
        # Assuming the first model is the one of interest
        model_info = models[0]
        model_coordinates_url = model_info.get("coordinates_url")
        if not model_coordinates_url:
            return {"error": f"No coordinates URL found for model in job {project_id}.", "model_info": model_info}

        # 4. Download PDB coordinates
        pdb_response = requests.get(model_coordinates_url, timeout=60)
        pdb_response.raise_for_status()

        pdb_content: str
        if model_coordinates_url.endswith('.gz') or pdb_response.headers.get('Content-Type', '').lower() in ['application/gzip', 'application/x-gzip']:
            try:
                with gzip.GzipFile(fileobj=BytesIO(pdb_response.content)) as gz_file:
                    pdb_content = gz_file.read().decode('utf-8')
            except gzip.BadGzipFile:
                return {"error": "Failed to decompress PDB file: Bad Gzip format."}
            except Exception as e:
                return {"error": f"Error decompressing PDB content: {e}"}
        else:
            pdb_content = pdb_response.text

        if not pdb_content.strip():
             return {"error": "Downloaded PDB content is empty."}

        # 5. Save PDB file
        os.makedirs(output_directory, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        pdb_filename = f"{output_pdb_filename_prefix}_{timestamp}.pdb"
        pdb_file_path = os.path.join(output_directory, pdb_filename)
        
        with open(pdb_file_path, 'w', encoding='utf-8') as f_pdb:
            f_pdb.write(pdb_content)

        return {
            "message": "Structure prediction successful.",
            "project_id": project_id,
            "gmqe": model_info.get("gmqe"),
            "model_details": {
                "model_id": model_info.get("model_id"),
                "status": model_info.get("status"),
                "qmean":model_info.get("qmean_global").get("avg_local_score"),
                "coordinate_url": model_info.get("coordinates_url"),
                "modelcif_url": model_info.get("modelcif_url"),
            },
            "pdb_file_path": pdb_file_path,
            "date_created": status_data.get("date_created"),
            "project_title": status_data.get("project_title"),
            "view_url": status_data.get("view_url"),
            
            # TODO: DELETE FOR NOW BECAUSE OG LARGE PDB
            # "pdb_content": pdb_content # Be cautious with large PDB files in memory
        }

    except requests.exceptions.HTTPError as http_err:
        error_details = http_err.response.text if http_err.response else "No response body"
        return {"error": f"HTTP error occurred: {http_err}. Status: {http_err.response.status_code if http_err.response else 'N/A'}", "details": error_details}
    except requests.exceptions.ConnectionError as conn_err:
        return {"error": f"Connection error: {conn_err}"}
    except requests.exceptions.Timeout as timeout_err:
        return {"error": f"Request timed out: {timeout_err}"}
    except requests.exceptions.RequestException as req_err:
        return {"error": f"An unexpected error occurred with the request: {req_err}"}
    except IOError as io_err:
        return {"error": f"File I/O error: {io_err}"}
    except Exception as e:
        return {"error": f"An unexpected error occurred in predict_structure: {e}"}

def main():
    """Command-line interface for predicting protein structure using SWISS-MODEL."""
    parser = argparse.ArgumentParser(
        description="Predict protein structure using the SWISS-MODEL API.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "sequence", 
        type=str, 
        help="The amino acid sequence for structure prediction."
    )
    parser.add_argument(
        "--api-token", 
        type=str, 
        default=os.environ.get("SWISS_MODEL_API_TOKEN"), 
        help="SWISS-MODEL API token. Reads from SWISS_MODEL_API_TOKEN environment variable if not provided."
    )
    parser.add_argument(
        "--output-dir", 
        type=str, 
        default=".", 
        help="Directory to save the output PDB file."
    )
    parser.add_argument(
        "--filename-prefix", 
        type=str, 
        default="predicted_structure", 
        help="Prefix for the output PDB filename."
    )
    parser.add_argument(
        "--title-prefix",
        type=str,
        default="moremi_biokit_prediction",
        help="Prefix for the project title submitted to SWISS-MODEL."
    )
    parser.add_argument(
        "--polling-interval",
        type=int,
        default=DEFAULT_POLLING_INTERVAL_SECONDS,
        help="Polling interval in seconds to check job status."
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=DEFAULT_MAX_POLLING_ATTEMPTS,
        help="Maximum number of polling attempts before timeout."
    )

    args = parser.parse_args()

    # Basic validation for token presence (either via arg or env var)
    if not args.api_token and not DEFAULT_SWISS_MODEL_API_TOKENS:
        parser.error("SWISS-MODEL API token is required. Provide it via --api-token or set the SWISS_MODEL_API_TOKEN environment variable, or ensure default tokens are available in the script.")
    elif not args.api_token:
        print("Warning: API token not provided via --api-token or environment variable. Attempting to use a default token (not recommended for production).")

    print(f"Starting structure prediction for sequence: {args.sequence[:15]}... Output dir: {args.output_dir}")

    results = predict_structure(
        sequence=args.sequence,
        api_token=args.api_token, # Will be None if not provided and no env var
        project_title_prefix=args.title_prefix,
        output_directory=args.output_dir,
        output_pdb_filename_prefix=args.filename_prefix,
        polling_interval=args.polling_interval,
        max_polling_attempts=args.max_attempts
    )

    print("--- Prediction Results ---")
    if results.get("status") == "error" or "error" in results: # Check both for safety
        print(f"Status: ERROR")
        print(f"Message: {results.get('error', 'Unknown error')}")
        if "details" in results:
            print(f"Details: {results['details']}")
    else:
        print(f"Status: Success")
        print(f"Message: {results.get('message', 'Prediction completed.')}")
        print(f"Project ID: {results.get('project_id')}")
        print(f"PDB File Path: {results.get('pdb_file_path')}")
        print("Model Details:")
        model_details = results.get("model_details", {})
        for key, value in model_details.items():
            if value is not None:
                 print(f"  - {key.replace('_', ' ').capitalize()}: {value}")

if __name__ == '__main__':
    # The example usage block is removed as the main() function now handles CLI execution.
    # If you want to run examples, you would now use the command line.
    main() 