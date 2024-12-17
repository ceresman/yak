"""
BLAST handler for sequence alignment and similarity scoring.
"""
import logging
from pathlib import Path
from typing import Dict, Generator, List, Tuple

import numpy as np
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

logger = logging.getLogger(__name__)

class BlastHandler:
    """Handles BLAST operations for sequence alignment and similarity scoring."""

    def __init__(self, output_dir: str = "data/processed/blast"):
        """Initialize BlastHandler with output directory."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def align_sequences(self, query_file: str, subject_file: str,
                       evalue: float = 1e-10) -> Generator[Dict, None, None]:
        """
        Perform BLAST alignment between query and subject sequences.

        Args:
            query_file: Path to query sequence file (FASTA)
            subject_file: Path to subject sequence file (FASTA)
            evalue: E-value threshold for filtering alignments

        Returns:
            Generator yielding alignment results
        """
        output_xml = self.output_dir / "blast_results.xml"

        try:
            # Configure BLAST command
            blastn_cline = NcbiblastnCommandline(
                query=query_file,
                subject=subject_file,
                outfmt=5,  # XML output format
                evalue=evalue,
                out=str(output_xml)
            )

            logger.info(f"Running BLAST alignment: {str(blastn_cline)}")
            stdout, stderr = blastn_cline()

            if stderr:
                logger.warning(f"BLAST stderr: {stderr}")

            # Parse results
            with open(output_xml) as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                for record in blast_records:
                    for alignment in record.alignments:
                        for hsp in alignment.hsps:
                            yield {
                                'query_id': record.query,
                                'subject_id': alignment.hit_def,
                                'identity': (hsp.identities / hsp.align_length) * 100,
                                'alignment_length': hsp.align_length,
                                'e_value': hsp.expect,
                                'bit_score': hsp.bits,
                                'query_start': hsp.query_start,
                                'query_end': hsp.query_end,
                                'subject_start': hsp.sbjct_start,
                                'subject_end': hsp.sbjct_end
                            }

        except Exception as e:
            logger.error(f"Error during BLAST alignment: {str(e)}")
            raise

    def calculate_similarity_scores(self, alignments: List[Dict]) -> pd.DataFrame:
        """
        Calculate similarity scores from alignment results.

        Args:
            alignments: List of alignment dictionaries

        Returns:
            DataFrame with similarity metrics
        """
        try:
            df = pd.DataFrame(alignments)

            # Calculate additional metrics
            df['similarity_score'] = df.apply(
                lambda row: (row['bit_score'] * row['identity']) / row['alignment_length'],
                axis=1
            )

            # Normalize scores
            df['normalized_score'] = (df['similarity_score'] - df['similarity_score'].min()) / \
                                   (df['similarity_score'].max() - df['similarity_score'].min())

            return df

        except Exception as e:
            logger.error(f"Error calculating similarity scores: {str(e)}")
            raise

    def annotate_genes(self, similarity_df: pd.DataFrame,
                      threshold: float = 0.8) -> Tuple[pd.DataFrame, Dict]:
        """
        Annotate genes based on similarity scores.

        Args:
            similarity_df: DataFrame with similarity scores
            threshold: Similarity score threshold for annotation

        Returns:
            Tuple of (annotated DataFrame, annotation statistics)
        """
        try:
            # Add annotation status
            similarity_df['is_significant'] = similarity_df['normalized_score'] >= threshold

            # Calculate annotation statistics
            stats = {
                'total_alignments': len(similarity_df),
                'significant_alignments': similarity_df['is_significant'].sum(),
                'mean_similarity': similarity_df['normalized_score'].mean(),
                'std_similarity': similarity_df['normalized_score'].std()
            }

            logger.info(f"Gene annotation complete. Statistics: {stats}")
            return similarity_df, stats

        except Exception as e:
            logger.error(f"Error annotating genes: {str(e)}")
            raise
