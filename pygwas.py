import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template
import numpy as np
import plotly.graph_objects as go

class MapGWASSNPs:
    def __init__(self, vcf_file_path: str, gwas_file_path: str, output_file_path: str, cut_off_qual: int = 20, filt_nr_disease: bool = True):
        self.vcf_file = vcf_file_path
        self.gwas_file = gwas_file_path
        self.output_file = output_file_path
        self.cut_off_qual = cut_off_qual
        self.filt_nr_disease = filt_nr_disease
        
        # Replace \ with / in the path
        output_file_path = output_file_path.replace('\\', '/')
        
        # If last character is /, remove it
        if output_file_path[-1] == '/':
            output_file_path = output_file_path[:-1]
        
        # Create the output directory if it does not exist
        if not os.path.exists(output_file_path):
            os.makedirs(output_file_path)
            os.makedirs(output_file_path + "/report", exist_ok=True)
            os.makedirs(output_file_path + "/report/fig", exist_ok=True)
            os.makedirs(output_file_path + "/report/data", exist_ok=True)
            
        # Path to the report directory
        self.report_path = output_file_path + "/report"
        self.report_fig_path = output_file_path + "/report/fig"
        self.report_data_path = output_file_path + "/report/data"
        
    def map_snps(self):
        report_data_path = self.report_data_path  # Path to the report data directory
        
        print("Step 1: Reading VCF file...")
        vcf_columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
        vcf_df = pd.read_csv(self.vcf_file, comment="#", sep=r'\s+', names=vcf_columns)

        print("Step 2: Reading GWAS catalog...")
        gwas_df = pd.read_csv(self.gwas_file)

        print("Step 3: Normalizing chromosome identifiers...")
        vcf_df["CHROM"] = vcf_df["CHROM"].str.replace('^chr', '', regex=True)
        gwas_df["CHR_ID"] = gwas_df["CHR_ID"].astype(str)

        print("Step 4: Ensuring position columns are strings...")
        vcf_df["POS"] = vcf_df["POS"].astype(str)
        gwas_df["CHR_POS"] = gwas_df["CHR_POS"].astype(str)

        print("Step 5: Merging VCF and GWAS data...")
        annotated_df = pd.merge(vcf_df, gwas_df, left_on=["CHROM", "POS"], right_on=["CHR_ID", "CHR_POS"], how="left")

        # Remove unmapped data and filter by quality
        annotated_df.dropna(subset=['DISEASE/TRAIT'], inplace=True)
        annotated_df = annotated_df[annotated_df["QUAL"] >= self.cut_off_qual]

        # Optionally filter out NR diseases
        if self.filt_nr_disease:
            annotated_df = annotated_df[annotated_df["DISEASE/TRAIT"] != "NR"]

        self.annotated_df = annotated_df

        # Save the annotated data to the specified directory
        output_path = os.path.join(report_data_path, f'in-house_report.csv')
        print("Saving annotated data to CSV...")
        annotated_df.to_csv(output_path, index=False)
        
        print(f"Annotated data saved to {output_path}")
        return annotated_df

    def prepare_report_data(self):
        df = self.annotated_df.copy()
        report_data = df[['DISEASE/TRAIT', 'RISK ALLELE FREQUENCY', 'P-VALUE', 'REGION', 'SNPS', 'MAPPED_GENE', 'Groups of Disease/Trait', 'Icon']].groupby(
            'DISEASE/TRAIT', as_index=False
        ).agg({
            'RISK ALLELE FREQUENCY': 'max',
            'P-VALUE': 'min',
            'REGION': 'first',
            'SNPS': 'first',
            'MAPPED_GENE': 'first',
            'Groups of Disease/Trait': 'first',
            'Icon': 'first'
        })

        # Convert to numeric and handle non-numeric values
        report_data['RISK ALLELE FREQUENCY'] = pd.to_numeric(report_data['RISK ALLELE FREQUENCY'], errors='coerce')

        # Multiply by 100 to get percentages
        report_data['RISK ALLELE FREQUENCY (%)'] = report_data['RISK ALLELE FREQUENCY'] * 100

        # Drop rows with missing or invalid values
        report_data.dropna(subset=['RISK ALLELE FREQUENCY (%)'], inplace=True)

        # Sort by risk allele frequency in descending order
        report_data.sort_values(by='RISK ALLELE FREQUENCY (%)', ascending=False, inplace=True)
        
        self.report_data = report_data
        
        # Save the report data to the specified directory
        output_path = os.path.join(self.report_data_path, f'report_data.csv')
        print("Saving report data to CSV...")
        report_data.to_csv(output_path, index=False)
        return report_data

    def generate_individual_bar_charts(self):
        data = self.report_data
        sns.set(style="whitegrid")
        for i, row in data.iterrows():
            fig = go.Figure()
            z_ = np.linspace(0, 100, 100)
            fig.add_trace(go.Heatmap(
                z=[z_],
                colorscale=[
                    [0, '#008AA5'],
                    [1, '#F1423E']
                ],
                showscale=False,
                colorbar=dict(
                    tick0=0,
                    dtick=1
                )
            ))

            fig.update_layout(
                xaxis=dict(
                    range=[0, 100]
                )
            )

            fig.update_layout(
                xaxis_title=f"Your Genetic Risk  {np.round(row['RISK ALLELE FREQUENCY (%)'], decimals=2)} (%)",
            )

            fig.add_trace(go.Scatter(
                x=[row['RISK ALLELE FREQUENCY (%)']],
                y=[0.9],
                mode='markers',
                marker=dict(
                    symbol='triangle-down',
                    size=30,
                    color='#434343'
                )
            ))

            fig.update_layout(
                width=1000,
                height=220
            )

            fig.update_layout(
                yaxis=dict(
                    visible=False,
                ),
                showlegend=False,
                xaxis=dict(
                    showgrid=False,
                    zeroline=False
                ),
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)'
            )

            fig.update_layout(
                font=dict(
                    size=16,  
                )
            )

            fig.update_layout(
                margin=dict(t=60, b=100, l=100, r=100, pad=0)
            )

            sanitized_name = row['DISEASE/TRAIT'].replace('/', '_').replace(' ', '_')
            figure_path = os.path.join(self.report_fig_path, f"report_{sanitized_name}.png")
            fig.write_image(figure_path)

    def generate_html_report(self):
        data = self.report_data
        output_path = os.path.join(self.report_path, 'KKUNPhI_report.html')
        chart_paths = []
        title_text = []
        region = []
        snps = []
        mapped_gene = []
        groups_of_disease_trait = []
        icon = []

        for i, row in data.iterrows():
            sanitized_name = row['DISEASE/TRAIT'].replace('/', '_').replace(' ', '_')
            chart_paths.append(f"fig/report_{sanitized_name}.png")
            title_text.append(row['DISEASE/TRAIT'])
            region.append(row['REGION'])
            snps.append(row['SNPS'])
            mapped_gene.append(row['MAPPED_GENE'])
            groups_of_disease_trait.append(row['Groups of Disease/Trait'])
            icon.append(row['Icon'])

        num_diseases = len(data['DISEASE/TRAIT'].unique())
        summary_text = f"This report includes a total of {num_diseases} unique diseases/traits analysed."

        html_template = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>GWAS Report</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 40px; }
                h1 { text-align: center; }
                img { display: block; margin: 20px auto; }
                .summary { margin-top: 50px; font-size: 18px; text-align: center; }
                h2 { text-align: left; margin-top: 50px; font-size: 24px; color: #434343;}
            </style>
        </head>
        <body>
            <h1>KKUNPhl WGS Analysis Report</h1>
            
            {% for chart_, title_, region_, snps_, mapped_gene_ in data_source %}
                <h2>{{ title_ }}</h2>
                <p><b>Gene region:</b> {{ region_ }}</p>
                <p><b>SNPs:</b> {{ snps_ }}</p>
                <p><b>Mapped Gene:</b> {{ mapped_gene_ }}</p>
                <img src="{{ chart_ }}" alt="Risk Allele Frequency Chart">
                <hr>
                <br>
            {% endfor %}
            
            <div class="summary">
                <p>{{ summary }}</p>
            </div>
        </body>
        </html>
        """

        data_source = zip(chart_paths, title_text, region, snps, mapped_gene)

        template = Template(html_template)
        rendered_html = template.render(
            data_source=data_source,
            summary=summary_text
        )

        with open(output_path, 'w') as f:
            f.write(rendered_html)

    def generate_report(self):
        self.prepare_report_data()
        self.generate_individual_bar_charts()
        self.generate_html_report()
        print(f"Report saved to {os.path.join(self.report_path, 'KKUNPhI_report.html')}")

if __name__ == '__main__':
    vcf_file = "Barcode10/medaka-variant-out/medaka.annotated.vcf"
    gwas_file = "gwas_catalog_grouped.csv"
    output_file = "test"

    mapper = MapGWASSNPs(vcf_file, gwas_file, output_file, cut_off_qual=0, filt_nr_disease=True)
    annotated_df = mapper.map_snps()
    mapper.generate_report()