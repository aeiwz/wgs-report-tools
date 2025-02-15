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
        
        # Normalize path format
        output_file_path = output_file_path.replace('\\', '/').rstrip('/')

        # Create necessary directories
        self.report_path = os.path.join(output_file_path, "report")
        self.report_data_path = os.path.join(self.report_path, "data")


        os.makedirs(self.report_data_path, exist_ok=True)

    def map_snps(self):
        print("Step 1: Reading VCF file...")
        vcf_columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
        if vcf_file.endswith('.gz'):
            vcf_df = pd.read_csv(self.vcf_file, compression='gzip', sep='\t', comment="#", names=vcf_columns)
        else:
            vcf_df = pd.read_csv(self.vcf_file, comment="#", sep='\t', names=vcf_columns)
        # Function to classify variant type
        def classify_variant(ref, alt):
            ref_len = len(ref)
            alt_len = len(alt)

            if ref_len == 1 and alt_len == 1:
                return "SNPs"
            elif ref_len < alt_len:
                return "INS"
            elif ref_len > alt_len:
                return "DEL"
            else:
                return "COMPLEX"
            
        
        
        # Apply classification to the DataFrame
        vcf_df["TYPE"] = vcf_df.apply(lambda row: classify_variant(row["REF"], row["ALT"]), axis=1)
        
        self.vcf_report = vcf_df.copy()
        print("VCF file PASS filter count: ", vcf_df[vcf_df["FILTER"] == "PASS"].shape[0])

        print("Step 2: Reading GWAS catalog...")
        #check compress 
        if self.gwas_file.endswith('.gz'):
            gwas_df = pd.read_csv(self.gwas_file, low_memory=False, compression='gzip')  # Prevent dtype warnings
        else:
            gwas_df = pd.read_csv(self.gwas_file, low_memory=False)
        print("GWAS catalog PASS : ", gwas_df.shape)

        print("Step 3: Normalizing chromosome identifiers...")
        #vcf_df["CHROM"] = vcf_df["CHROM"].astype(str).str.replace('^chr', '', regex=True)
        #gwas_df["CHR_ID"] = gwas_df["CHR_ID"].astype(str)
        print('Normalizing chromosome identifiers PASS')

        print("Step 4: Ensuring position columns are strings...")
        vcf_df["POS"] = vcf_df["POS"].astype(str)
        gwas_df["CHR_POS"] = gwas_df["CHR_POS"].astype(str)
        print('Ensuring position columns are strings PASS')

        print("Step 5: Merging VCF and GWAS data...")
        annotated_df = pd.merge(vcf_df, gwas_df, left_on=["CHROM", "POS"], right_on=["CHR_ID", "CHR_POS"], how="inner")
        print('Merging VCF and GWAS data PASS')

        annotated_df.dropna(subset=['DISEASE/TRAIT'], inplace=True)
        annotated_df = annotated_df[annotated_df["QUAL"] >= self.cut_off_qual]

        if self.filt_nr_disease:
            annotated_df = annotated_df[annotated_df["DISEASE/TRAIT"] != "NR"]

        self.annotated_df = annotated_df

        # Save output
        output_path = os.path.join(self.report_data_path, 'in-house_report.csv')
        print("Saving annotated data to CSV...")
        annotated_df.to_csv(output_path, index=False)

        print(f"Annotated data saved to {output_path}")
        return annotated_df

    def prepare_report_data(self):
        df = self.annotated_df.copy()
        report_data = df.groupby('DISEASE/TRAIT', as_index=False).agg({
            'RISK ALLELE FREQUENCY': 'max',
            'P-VALUE': 'min',
            'REGION': 'first',
            'SNPS': 'first',
            'MAPPED_GENE': 'first',
            'Groups of Disease/Trait': 'first',
            'MAPPED_TRAIT_URI': 'first',
            'MAPPED_TRAIT_DESCRIPTION': 'first'
        })

        report_data['RISK ALLELE FREQUENCY'] = pd.to_numeric(report_data['RISK ALLELE FREQUENCY'], errors='coerce')
        report_data['RISK ALLELE FREQUENCY (%)'] = report_data['RISK ALLELE FREQUENCY'] * 100
        report_data.dropna(subset=['RISK ALLELE FREQUENCY (%)'], inplace=True)
        report_data.sort_values(by='RISK ALLELE FREQUENCY (%)', ascending=False, inplace=True)

        self.report_data = report_data

        output_path = os.path.join(self.report_data_path, 'report_data.csv')
        print("Saving report data to CSV...")
        report_data.to_csv(output_path, index=False)
        return report_data



    def generate_html_report(self):
        data = self.report_data
        output_path = os.path.join(self.report_path, 'GWAS_report.html')
        details = []
        embedded_svgs = []
        icon = []
        import plotly.io as pio
        for _, row in data.iterrows():
            # Generate a Plotly figure
            
            import plotly.graph_objects as go
            import numpy as np

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
            #fig = go.Figure(data=[go.Bar(x=[row['DISEASE/TRAIT']], y=[row['RISK ALLELE FREQUENCY (%)']])])
            #fig.update_layout(title=row['DISEASE/TRAIT'], xaxis_title="Trait", yaxis_title="Risk Allele Frequency (%)")

            # Convert figure to SVG
            svg_bytes = pio.to_image(fig, format="svg")
            svg_string = svg_bytes.decode("utf-8")  # Convert bytes to string
            #svg_string = fig.to_html(full_html=False, include_plotlyjs='cdn')
            # Store details and embedded SVG
            with open(f"data/Group of disease traits/{row['Groups of Disease/Trait']}.svg", "r", encoding="utf-8") as file:
                img = file.read()
                
                
            details.append((row['DISEASE/TRAIT'], row['REGION'], row['SNPS'], row['MAPPED_GENE'], row['Groups of Disease/Trait'], row['MAPPED_TRAIT_DESCRIPTION']))
            embedded_svgs.append(svg_string)
            icon.append(str(img))
        
        with open('data/logo/KKUNPhI-01.svg', 'r', encoding='utf-8') as file:
            logo = file.read()

        summary_text = f"This report includes {len(data['DISEASE/TRAIT'].unique())} unique diseases/traits analyzed."
        vcf_report = self.vcf_report
        total_variant = vcf_report.shape[0]
        type_of_varients = vcf_report['TYPE'].value_counts()
        type_of_varients_percentage = np.round(type_of_varients / total_variant * 100, decimals=2)


        import plotly.graph_objects as go

        svg_circle = []

        for i in type_of_varients_percentage.index:    
            # Data for the single progress bar
            value = type_of_varients_percentage[i]  # Percentage of progress
            color = {'COMPLEX':'FF9999', 'DEL':'FF7F3E', 'INS':'3D527D', 'SNPs':'FFB854'}  # Custom color for progress

            #description = "DUIS AUTE IRURE DOL ORIN REPREHENDERIT IN VELIT ESSE CILLUM FUGIAT NULLA PARI"

            # Create the pie chart (donut chart)
            fig2 = go.Figure(go.Pie(
                values=[value, 100 - value],  # Progress and remaining part
                labels=["", ""],  # Hide labels
                hole=0.6,  # Creates a donut shape
                direction="clockwise",
                marker=dict(colors=[color[i], "lightgray"]),  # Custom color and gray for remaining
                textinfo="none",  # Hide Pie chart labels
                showlegend=False
            ))

            # Add percentage text inside the donut chart
            fig2.add_annotation(
                text=f"<b>{value}%</b>", showarrow=False,
                xref="paper", yref="paper",
                x=0.5, y=0.5, font=dict(size=24)
            )

            # Add label below the chart
            fig2.add_annotation(
                text=f'<b>{type_of_varients[i]:,.0f} Positions</b>', showarrow=False,
                xref="paper", yref="paper",
                x=0.5, y=-0.38, font=dict(size=22)
            )

            # Layout adjustments
            fig2.update_layout(
                title=f"<b>{i}</b>",
                height=400, width=400,
                showlegend=False
            )

            #Set title to center and below of the chart
            fig2.update_layout(
                title_x=0.5,
                title_y=0.1,
                font=dict(size=18)
            )

            #Background transparency
            fig2.update_layout(
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)'
            )




            # Display the figure
            svg_bytes = pio.to_image(fig2, format="svg")
            svg_string = svg_bytes.decode("utf-8")  # Convert bytes to string
            svg_circle.append(svg_string)

        

        # HTML Template
        html_template = r"""
            <!DOCTYPE html>
            <html lang="en">
            <head>
                <meta charset="UTF-8">
                <meta name="viewport" content="width=device-width, initial-scale=1.0">
                <title>GWAS Report</title>
                <link href="https://fonts.googleapis.com/css2?family=Poppins:wght@300;400;600&display=swap" rel="stylesheet">
                <style>
                    body {
                        font-family: 'Poppins', sans-serif;
                        background: #f7f9fc;
                        color: #333;
                        margin: 40px;
                        line-height: 1.6;
                    }
                    .header-container {
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        margin-bottom: 20px;
                        flex-wrap: wrap;
                        flex-direction: row;
                    }
                    h1 {
                        text-align: center;
                        flex-grow: 8;
                        background: linear-gradient(to right, #3c79aa, #FF6B6B);
                        color: white;
                        padding: 20px;
                        border-radius: 8px;
                        box-shadow: 0 4px 10px rgba(0, 0, 0, 0.2);
                        font-size: 45px;
                        margin: 0;
                    }
                    h2 {
                        margin-top: 40px;
                        color: #2D3B71;
                        border-bottom: 3px solid #2D3B71;
                        padding-bottom: 5px;
                        display: inline-block;
                        font-size: 32px;
                    }
                    p {
                        margin: 10px 0;
                        font-size: 22px;
                        color: rgb(87, 107, 145);
                    }
                    p2 {
                        margin: 10px 0;
                        font-size: 16px;
                        color: gray;
                        font-style: italic;
                    }
                    b {
                        color: #2D3B71;
                    }
                    .icon-text {
                        display: flex;
                        align-items: center;
                        margin-top: 30px;
                        flex-direction: row;
                    }
                    .icon-text svg {
                        width: 250px;
                        height: auto;
                        border-radius: 50%;
                        margin-right: 20px;
                        display: flex;
                        flex-direction: row;
                        }
                    .chart-container {
                        text-align: left;
                        background: white;
                        padding: 20px;
                        border-radius: 8px;
                        box-shadow: 0 4px 10px rgba(0, 0, 0, 0.1);
                        margin: 20px 0;
                        width: 100%;  /* Adjust based on the window size */
                        max-width: 98%; /* Set a max-width to prevent it from getting too big */
                    }
                    
                    .chart svg {
                        width: 95%;  /* Make SVG scale with container */
                        height: auto; /* Maintain aspect ratio */
                    }
                    .download-btn {
                        display: inline-block;
                        background: rgb(87, 107, 145);
                        color: white;
                        padding: 8px 12px;
                        font-size: 14px;
                        border-radius: 5px;
                        border: none;
                        cursor: pointer;
                        transition: 0.3s;
                        box-shadow: 0 4px 10px rgba(0, 0, 0, 0.1);
                    }
                    .download-btn:hover {
                        background: rgb(87, 107, 145);
                    }
                    .logo-container {
                        text-align: center;
                        margin-top: 0px;
                        margin-bottom: 0px;
                        flex-grow: 2;
                        flex: 250px;
                        padding: 80px;
                    }
                    hr {
                        border: none;
                        height: 2px;
                        background: #ddd;
                        margin: 40px 0;
                    }
                    .title {
                        font-size: 200px;
                        font-weight: bold;
                        color: #2D3B71;
                        text-align: center;
                        font-family: 'Poppins', sans-serif;
                        margin-top: 250px;
                        margin-bottom: 20px;
                    }
                    .subtitle {
                        font-size: 34px;
                        color: black;
                        text-align: center;
                    }
                    .variants {
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        flex-direction: row;
                        flex-flow: space-evenly;
                        margin: 20px;
                        flex-wrap: wrap;
                    }
                    .variant {
                        display: flex;
                        align-items: center;
                        gap: 5px;
                        margin-bottom: 5px;
                        align-self: auto | flex-start | flex-end | center | baseline | stretch;
                        flex: 1 1 auto;
                        flex-wrap: wrap;
                        width: 300px;
                    }
                    .position_infomation {
                        display: flex;
                        align-items: center;
                        justify-content: flex-start;
                        gap: 20px;
                    }

                    .icon-container {
                        flex-shrink: 0; /* Prevents icon from shrinking */
                    }

                    .position_inside {
                        display: flex;
                        flex-direction: column;
                        gap: 5px;
                    }

                    }
                </style>
                                <script>
                    function resizeChart() {
                        var chart = document.getElementById('chart_1');
                        Plotly.relayout(chart, {
                            width: window.innerWidth * 0.9,  // 90% of the window width
                            height: window.innerHeight * 0.3 // Adjust height dynamically
                        });
                    }
                
                    window.addEventListener('resize', resizeChart);
                    window.onload = resizeChart; // Run once when the page loads
                </script>
                
            </head>
            <body>
                <div class="header-container">
                    <div class="logo-container">
                        <div class="logo" id="logo_01">
                            {{ logo_ | safe }}
                        </div>
                    </div>
                        
                        <h1>Whole Genome Analysis Report</h1>
                </div>
                
                <div>
                    <h2>Disclaimer</h2>
                        <p>The service provided by Khon Kaen University National
                            Phenome Institute (KKUNPhI) are currently for research use
                            only. As they have not been submitted for review to any
                            regulatory agency/notified body for clinical diagnostics, caution
                            must be excercisrd when describing the application of KKUNPhI
                            service in the clinical research field.
                        </p>
                    <p2>* This report is for research use only</p2>
                    </div>
                    <div>
                    <h2>Overview</h2>
                        <p>This report presents the findings of whole genome sequencing
                            analysis. It provides detailed associations between genetic
                            variants and various diseases and traits, highlighting the mapped
                            genes, relevant SNPs, and associated gene regions.
                        </p>
                    </div>
                    <div>
                        <h2>Termonology</h2>
                            <p><b>Gene:</b> Gene is a segment of DNA that serves as a blueprint for producing
                                proteins or functional RNA molecules, which carry out vital biological
                                functions in the body. Genes are the basic units of heredity, passed
                                from parents to offspring, and play a crucial role in determining traits
                                and regulating cellular processes.
                            </p>
                            <p><b>Single Nucleotide Polymorphism (SNP):</b> A SNP is a variation in a single nucleotide in the DNA sequence that
                                occurs at a specific position in the genome and is common in a
                                population. SNPs can influence traits, disease susceptibility, and drug
                                response.
                            </p>
                            <p><b>Insertion (INS):</b> An insertion is the addition of one or more nucleotides to a DNA sequence.</p>
                            <p><b>Deletion (DEL):</b> A deletion is the removal of one or more nucleotides from a DNA sequence.</p>
                            <p><b>Complex (COMPLEX):</b> A complex variant is a combination of insertions, deletions, and substitutions</p>
                            <p><b>Alleles:</b> An allele is a variant form of a gene found at a specific position (locus)
                                on a chromosome. Each allele is inherited, one from each parent.
                            </p>
                            <p><b>Mapped Gene:</b> Genes mapped near or overlapping the SNPs.
                            </p>
                            <p><b>Chromosomal region:</b> The genomic region associated with the trait or disease.
                            </p>
                            <p><b>Risk Allele Frequency (%):</b> The frequency of the risk allele in the population.
                            </p>

                            
                            <div class="title">{{count_variant}}</div>
                            <div class="subtitle">Variants have been found in your genome.</div>

                            <div class="variants">
                            <div class="variant">
                                {{variant_1 | safe}}
                            </div>
                            <div class="variant">
                                {{variant_2 | safe}}
                            </div>
                            <div class="variant">
                                {{variant_3 | safe}}
                            </div>
                            <div class="variant">
                                {{variant_4 | safe}}
                            </div>
                        </div>
                    </div>
                {% for (title_, region_, snps_, mapped_gene_, group_trait_, description_trait_), svg_, icon_ in data_source %}
                    <div class="chart-container">
                    <h2>{{ title_ }}</h2>
                    <p>{{ description_trait_ }}</p>

                    <div class="position_infomation">
                        <div class="icon-text">
                            <div class="icon-container">
                                {{icon_ | safe }}
                            </div>
                            <div class="position_inside">
                                <p><b>Region:</b> {{ region_ }}</p>
                                <p><b>SNPs ID:</b> {{ snps_ }}</p>
                                <p><b>Mapped Gene:</b> {{ mapped_gene_ }}</p>
                                <p><b>Group of disease/trait:</b> {{ group_trait_ }}</p>
                            </div>
                        </div>
                    </div>

                    
                        <div class="chart" id="chart_{{ loop.index }}">
                            {{ svg_ | safe }}
                        </div>
                        <button class="download-btn" onclick="downloadChart('{{ loop.index }}', '{{ title_ }}')">Download Chart</button>
                    </div>
                   
                    <hr>
                {% endfor %}

                <script>
                    function downloadChart(chartId, title) {
                        let svgElement = document.querySelector("#chart_" + chartId + " svg");
                        if (!svgElement) {
                            alert("SVG not found!");
                            return;
                        }

                        let serializer = new XMLSerializer();
                        let svgString = serializer.serializeToString(svgElement);
                        
                        let canvas = document.createElement("canvas");
                        let ctx = canvas.getContext("2d");
                        let img = new Image();
                        let svgBlob = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
                        let url = URL.createObjectURL(svgBlob);

                        img.onload = function () {
                            canvas.width = img.width;
                            canvas.height = img.height;
                            ctx.drawImage(img, 0, 0);
                            URL.revokeObjectURL(url);

                            let pngUrl = canvas.toDataURL("image/png");

                            let downloadLink = document.createElement("a");
                            downloadLink.href = pngUrl;
                            downloadLink.download = title.replace(/\s+/g, "_") + ".png";
                            document.body.appendChild(downloadLink);
                            downloadLink.click();
                            document.body.removeChild(downloadLink);
                        };

                        img.onerror = function () {
                            alert("Failed to load the SVG. Please check for unsupported elements.");
                        };

                        img.src = url;
                    }
                </script>

            </body>
            </html>
        """

        template = Template(html_template)
        rendered_html = template.render(data_source=zip(details, embedded_svgs, icon), count_variant=f'{total_variant:,.0f}',
                                        variant_1 = svg_circle[0], variant_2 = svg_circle[1], variant_3 = svg_circle[2], variant_4 = svg_circle[3],
                                        logo_=logo)

        with open(output_path, 'w') as f:
            f.write(rendered_html)

        #print(f"Report saved to {output_path}")

    def generate_report(self):
        self.prepare_report_data()
        self.generate_html_report()
        print(f"Report saved to {os.path.join(self.report_path, 'GWAS_report.html')}")

if __name__ == '__main__':
    vcf_file = "data/medaka.annotated.vcf.gz"
    #gwas_file = "https://raw.githubusercontent.com/aeiwz/wgs-report-tools/main/data/gwas_catalog_grouped.csv"
    gwas_file = "data/gwas_database_with_description_expanded.csv.gz"
    output_file = "test"

    mapper = MapGWASSNPs(vcf_file, gwas_file, output_file, cut_off_qual=30, filt_nr_disease=True)
    mapper.map_snps()
    mapper.generate_report()