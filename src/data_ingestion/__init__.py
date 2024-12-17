"""
Data ingestion module for processing genomic data files.
"""
from .fasta_handler import FastaHandler
from .vcf_handler import VcfHandler

__all__ = ['FastaHandler', 'VcfHandler']
