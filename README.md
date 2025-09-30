
# **mapgwas**

MapGWASSNPs — tool to:
1. Read a VCF  
2. Intersect with a GWAS catalog file  
3. Produce a polished **interactive HTML report**

---

## **Installation**
```bash
pip install .
# or for dev
pip install -e .


# **Usage Guide**

Quick start:
```bash
mapgwas --vcf input.vcf --gwas gwas.csv.gz --out outdir --qual-cutoff 50
