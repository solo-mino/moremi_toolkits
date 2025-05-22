import iedb
from typing import Union
import requests

def predict_immunogenicity(sequence: Union[str, list], method: str = "NetMHCIIpan-4.3" , allele: Union[str, list]="HLA-DRB1*07:01", length: Union[int, list]= None):
    '''
    Query IC50 values for MHC class II peptide binding of either a single amino 
    acid sequence or set of sequences to an MHC class II molecule. Sends POST 
    request to IEDB API. 

    IEDB official web tool can be found here: http://tools.iedb.org/mhcii/

            Parameters:
                    method (str): MHC class II binding prediction method. 
                        Available argument options (format: method-version):
                            - recommended
                            - consensus-2.22
                            - consensus-2.18
                            - netmhciipan-4.0
                            - netmhciipan-3.2
                            - netmhciipan-3.1
                            - smm_align-1.1
                            - nn_align-2.3
                            - nn_align-2.2
                            - comblib-1.0
                            - tepitope-1.0
                        More information on prediction methods can be found 
                        here: http://tools.iedb.org/mhcii/help/#Method

                    sequence ([str, list]): Peptide amino acid sequence or
                        list of sequences.
                    allele ([str, list]): Single or multiple HLA alleles. 
                        Provide single alleles as a string. Multiple alleles
                        should be formatted in a comma separated list. Example
                        below:
                            Ex.) HLA-DRB1*01:01
                        To run alpha and beta chains together format as such:
                            Ex.) DPA1*01/DPB1*04:01"
                    length ([int, list]): Peptide length or list of lengths

            Returns:
                dict[str, Union[str, int, dataframe]]
                    dataframe (pandas.DataFrame): Tabular results formatted
                        as pandas.DataFrame, immunogenicity score, no. of strong binding, no. of weak binding
    '''
    
    valid_methods = ['recommended', 'consensus-2.22', 'consensus-2.18', 'netmhciipan-4.0',
                     'netmhciipan-3.2', 'netmhciipan-3.1', 'NetMHCIIpan-4.3', "smm_align-1.1", "nn_align-2.3", "nn_align-2.2","comblib-1.0", "tepitope-1.0"]
    if method not in valid_methods:
        raise ValueError(f"Invalid prediction method: {method}. Valid options are {valid_methods}")
    try:
        # Send POST request to MHC class II peptide binding prediction tool:
        mhcii_res = iedb.query_mhcii_binding(method=method, sequence=sequence, allele=allele, length=length)
        
    except requests.exceptions.RequestException as e:
        if isinstance(e, ConnectionError):
            error_msg = "Network connection error: Unable to connect to IEDB API."
        elif isinstance(e, TimeoutError):
            error_msg = "Request timed out: API response took too long."
        else:
            error_msg = f"Request failed: {e}"
        return {"error_message": error_msg, "result": None}

    except ValueError as e:
        return {"error_message": f"Error parsing response: {e}", "result": None}
    except Exception as e:
        print(f"An error occurred: {e}")
        return {"error_message": f"Error parsing response: {e}", "result": None}
    
    
    
    # Calculate binding counts based on IC50 values
    if mhcii_res is None:
        return {"error": "No immunogenicity result, API response is empty"}
    
    # Initialize counters
    strong_binding = 0
    moderate_binding = 0
    weak_binding = 0
    
    if not mhcii_res.empty:
        for ic50 in mhcii_res['ic50']:
            try:
                ic50_value = float(ic50)
                if ic50_value < 50:
                    strong_binding += 1
                elif 50 <= ic50_value <= 500:
                    moderate_binding += 1
                else:
                    weak_binding += 1
            except (ValueError, TypeError):
                continue
    
    # Calculate immunogenicity score
    total_predictions = len(mhcii_res)
    if total_predictions > 0:
        # Calculate raw score (proportion of strong binders)
        raw_score = strong_binding / total_predictions
        # Normalize to 0-1 scale
        immunogenic_score = raw_score  # Already normalized since it's a proportion
    else:
        immunogenic_score = 0.0
    
    return {
        "df": mhcii_res,
        "immunogenic_score": immunogenic_score,
        "strong_binding": strong_binding,
        "moderate_binding": moderate_binding,
        "weak_binding": weak_binding,
    }
    
    
    
    