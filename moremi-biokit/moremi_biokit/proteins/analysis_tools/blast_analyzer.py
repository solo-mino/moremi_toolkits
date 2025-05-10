"""
Performs BLAST searches for protein sequences.

This module provides functionality to query NCBI BLAST databases with protein sequences.
"""

from typing import Any, List, Dict, Optional

from Bio.Blast import NCBIWWW, NCBIXML

def perform_blast(sequence: str, database: str ="nr", program: str ="blastp", num_of_alignments: int = 3) -> Optional[Dict[str, Any]]:
    """Perform a BLAST search on a protein sequence.

    Args:
        sequence (str): The protein sequence to search.
        database (str, optional): The BLAST database to use. Defaults to "nr".
        program (str, optional): The BLAST program to use. Defaults to "blastp".
            Other options include blastn, blastx, tblastn, or tblastx.
        num_of_alignments (int, optional): The number of alignments to return. 
            Defaults to 3.

    Returns:
        Optional[Dict[str, Any]]: A dictionary containing the BLAST results 
            (e.g., {"blast_result": [...]}) or None if an error occurs or no 
            significant alignments are found.
            If no alignments are found, it returns a dictionary like:
            {"blast_result": "No significant alignments found.", "alignments": []}

    Raises:
        None: Errors during BLAST execution are caught and printed; None is returned.
    """
    blast_records = None
    try:
        sequence = str(sequence).strip()
        result_handle = NCBIWWW.qblast(program, database, sequence)
        # It's good practice to close the handle, though NCBIXML.parse might do it.
        # Using a try/finally block ensures it's closed.
        try:
            result = NCBIXML.parse(result_handle)
            blast_records = next(result) # Get the first BLAST record
        finally:
            result_handle.close()

    except Exception as e:
        print(f"Error performing BLAST: {e}")
        return None
    
    if not blast_records or not blast_records.alignments:
        return {"blast_result": "No significant alignments found.", "alignments": []}

    blast_result_list: List[Dict[str, Any]] = []
    for alignment in blast_records.alignments[:num_of_alignments]:
        hsp = alignment.hsps[0] # Highest scoring pair
        
        record: Dict[str, Any] = {}
        record['sequence_title'] = alignment.title
        record['length'] = alignment.length
        record['e_value'] = hsp.expect
        # Ensure hsp.align_length is not zero to avoid DivisionByZeroError
        if hsp.align_length > 0:
            record['identity_percentage'] = (hsp.identities / hsp.align_length) * 100
        else:
            record['identity_percentage'] = 0.0
        record['matched_sequence_alignment'] = hsp.sbjct
        
        blast_result_list.append(record)
    
    return {"blast_result": blast_result_list} 