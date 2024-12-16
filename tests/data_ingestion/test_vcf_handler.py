"""
Tests for VCF file handling functionality.
"""
import os
from pathlib import Path
import pytest
import numpy as np
import pandas as pd

from src.data_ingestion.vcf_handler import VCFHandler

@pytest.fixture
def sample_vcf(tmp_path):
    """Create a sample VCF file for testing."""
    vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
chr1\t100\trs1\tA\tT\t100\tPASS\t.\tGT\t0/1\t1/1
chr1\t200\trs2\tG\tC\t100\tPASS\t.\tGT\t0/0\t0/1
"""
    vcf_file = tmp_path / "test.vcf"
    vcf_file.write_text(vcf_content)
    return str(vcf_file)

@pytest.fixture
def vcf_handler():
    """Create VCFHandler instance."""
    return VCFHandler()

def test_vcf_handler_init(vcf_handler):
    """Test VCFHandler initialization."""
    assert isinstance(vcf_handler, VCFHandler)

def test_load_vcf(vcf_handler, sample_vcf):
    """Test loading VCF file."""
    variants = vcf_handler.load_vcf(sample_vcf)
    assert len(variants) == 2
    assert variants[0]['CHROM'] == 'chr1'
    assert variants[0]['POS'] == 100

def test_filter_variants(vcf_handler, sample_vcf):
    """Test variant filtering by quality."""
    variants = vcf_handler.load_vcf(sample_vcf)
    filtered = vcf_handler.filter_variants(variants, min_qual=90)
    assert len(filtered) == 2
    filtered = vcf_handler.filter_variants(variants, min_qual=150)
    assert len(filtered) == 0


def test_extract_genotypes(vcf_handler, sample_vcf):
    """Test genotype extraction."""
    variants = vcf_handler.load_vcf(sample_vcf)
    genotypes = vcf_handler.extract_genotypes(variants)
    assert genotypes.shape == (2, 2)  # 2 variants x 2 samples
    assert genotypes.iloc[0, 1] == 2  # homozygous alt (1/1)
    assert genotypes.iloc[1, 0] == 0  # homozygous ref (0/0)

def test_error_handling(vcf_handler, tmp_path):
    """Test error handling for invalid files."""
    invalid_file = tmp_path / "nonexistent.vcf"
    with pytest.raises(FileNotFoundError):
        vcf_handler.load_vcf(str(invalid_file))
