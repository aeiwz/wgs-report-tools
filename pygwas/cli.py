import argparse
import sys
from .pygwas import MapGWASSNPs

def build_parser():
    p = argparse.ArgumentParser(
        prog="mapgwas",
        description="Map VCF variants to GWAS catalog and generate an HTML report"
    )
    p.add_argument("--vcf", required=True, help="VCF path (.vcf or .vcf.gz)")
    p.add_argument("--gwas", required=True, help="GWAS CSV file (CSV or CSV.GZ)")
    p.add_argument("--out", required=True, help="Output root directory")
    p.add_argument("--qual-cutoff", type=float, default=20.0, help="QUAL cutoff (default=20)")
    p.add_argument("--keep-nr", action="store_true", help="Keep rows with DISEASE/TRAIT == NR")
    return p

def main(argv=None):
    args = build_parser().parse_args(argv or sys.argv[1:])
    mapper = MapGWASSNPs(
        vcf_file_path=args.vcf,
        gwas_file_path=args.gwas,
        output_file_path=args.out,
        cut_off_qual=args.qual_cutoff,
        filt_nr_disease=not args.keep_nr
    )
    mapper.map_snps()
    mapper.generate_report()
    return 0
