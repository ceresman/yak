"""
Integration tests for the complete yak breeding pipeline.
"""
import os
from pathlib import Path
import pytest
import numpy as np
import pandas as pd

from src.data_ingestion.fasta_handler import FastaHandler
from src.data_ingestion.vcf_handler import VCFHandler
from src.gene_marker.blast_handler import BlastHandler
from src.feature_analysis.linkage_analyzer import LinkageAnalyzer
from src.feature_analysis.selection_analyzer import SelectionAnalyzer
from src.feature_analysis.snp_analyzer import SNPAnalyzer
from src.genetic_diversity.diversity_analyzer import DiversityAnalyzer
from src.visualization.plot_generator import PlotGenerator

@pytest.fixture
def sample_data(tmp_path):
    """Create sample data for pipeline testing."""
    # Create sample FASTA
    fasta_content = """>marker1
ATCGATCGATCG
>marker2
GCTAGCTAGCTA
"""
    fasta_file = tmp_path / "markers.fasta"
    fasta_file.write_text(fasta_content)

    # Create sample VCF
    vcf_content = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\tSAMPLE2
chr1\t100\trs1\tA\tT\t100\tPASS\t.\tGT\t0/1\t1/1
chr1\t200\trs2\tG\tC\t100\tPASS\t.\tGT\t0/0\t0/1
"""
    vcf_file = tmp_path / "variants.vcf"
    vcf_file.write_text(vcf_content)

    return {
        'fasta_file': str(fasta_file),
        'vcf_file': str(vcf_file)
    }

def test_complete_pipeline(sample_data):
    """Test complete pipeline execution."""
    # Initialize components
    fasta_handler = FastaHandler()
    vcf_handler = VCFHandler()
    blast_handler = BlastHandler()
    linkage_analyzer = LinkageAnalyzer()
    selection_analyzer = SelectionAnalyzer()
    snp_analyzer = SNPAnalyzer()
    diversity_analyzer = DiversityAnalyzer()
    plot_generator = PlotGenerator()

    # Load and process marker sequences
    markers = fasta_handler.load_fasta(sample_data['fasta_file'])
    assert len(markers) == 2

    # Load and process variants
    variants = vcf_handler.load_vcf(sample_data['vcf_file'])
    assert len(variants) == 2

    # Extract genotypes
    genotypes = vcf_handler.extract_genotypes(variants)
    assert genotypes.shape == (2, 2)

    # Perform feature analysis
    ld_results = linkage_analyzer.calculate_ld(genotypes)
    assert isinstance(ld_results, pd.DataFrame)

    selection_stats = selection_analyzer.calculate_selection_statistics(genotypes)
    assert isinstance(selection_stats, dict)

    snp_results = snp_analyzer.evaluate_trait_snps(
        genotypes,
        pd.DataFrame({'trait': [1, 0]}, index=['SAMPLE1', 'SAMPLE2'])
    )
    assert isinstance(snp_results, pd.DataFrame)

    # Calculate diversity indices
    diversity_indices = diversity_analyzer.calculate_diversity_indices(genotypes)
    assert isinstance(diversity_indices, dict)

    # Generate visualization
    plot_file = plot_generator.plot_diversity_indices(
        pd.DataFrame({
            'population': ['pop1', 'pop2'],
            'heterozygosity_observed': [0.5, 0.6],
            'heterozygosity_expected': [0.4, 0.5]
        }),
        'test_diversity.png'
    )
    assert os.path.exists(plot_file)

def test_error_handling_pipeline(tmp_path):
    """Test pipeline error handling with invalid inputs."""
    invalid_file = tmp_path / "nonexistent.fasta"

    fasta_handler = FastaHandler()
    with pytest.raises(FileNotFoundError):
        fasta_handler.load_fasta(str(invalid_file))

    vcf_handler = VCFHandler()
    with pytest.raises(FileNotFoundError):
        vcf_handler.load_vcf(str(invalid_file))
