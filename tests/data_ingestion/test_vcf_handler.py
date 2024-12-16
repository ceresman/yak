"""
Tests for VCF file handler.
"""
import pytest
import pandas as pd
from pathlib import Path
from src.data_ingestion.vcf_handler import VcfHandler

@pytest.fixture
def sample_vcf(tmp_path):
    """Create a sample VCF file for testing."""
    vcf_file = tmp_path / "test.vcf"
    with open(vcf_file, 'w') as f:
        f.write("""##fileformat=VCFv4.2
##INFO=<ID=Gene,Number=1,Type=String,Description="Gene name">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1000\trs1\tA\tG\t100\tPASS\tGene=GENE1
1\t2000\trs2\tC\tT\t20\tPASS\tGene=GENE2
""")
    return vcf_file

@pytest.fixture
def vcf_handler():
    """Create a VcfHandler instance."""
    return VcfHandler(min_quality=30.0)

def test_load_vcf_valid_file(vcf_handler, sample_vcf):
    """Test loading a valid VCF file."""
    vcf_reader = vcf_handler.load_vcf(sample_vcf)
    assert vcf_reader is not None

def test_load_vcf_missing_file(vcf_handler):
    """Test loading a non-existent file."""
    with pytest.raises(FileNotFoundError):
        vcf_handler.load_vcf(Path("nonexistent.vcf"))

def test_filter_variants(vcf_handler, sample_vcf):
    """Test filtering variants based on quality."""
    vcf_reader = vcf_handler.load_vcf(sample_vcf)
    variants = vcf_handler.filter_variants(vcf_reader)
    assert len(variants) == 1  # Only one variant passes quality filter
    assert variants.iloc[0]['ID'] == 'rs1'

def test_organize_by_gene(vcf_handler, sample_vcf):
    """Test organizing variants by gene."""
    vcf_reader = vcf_handler.load_vcf(sample_vcf)
    variants = vcf_handler.filter_variants(vcf_reader)
    gene_variants = vcf_handler.organize_by_gene(variants)
    assert len(gene_variants) == 1  # Only one gene after quality filter
    assert 'GENE1' in gene_variants
