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
        self.report_fig_path = os.path.join(self.report_path, "fig")
        self.report_data_path = os.path.join(self.report_path, "data")

        os.makedirs(self.report_fig_path, exist_ok=True)
        os.makedirs(self.report_data_path, exist_ok=True)

    def map_snps(self):
        print("Step 1: Reading VCF file...")
        vcf_columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
        vcf_df = pd.read_csv(self.vcf_file, comment="#", sep=r'\s+', names=vcf_columns)

        print("Step 2: Reading GWAS catalog...")
        gwas_df = pd.read_csv(self.gwas_file, low_memory=False)  # Prevent dtype warnings

        print("Step 3: Normalizing chromosome identifiers...")
        vcf_df["CHROM"] = vcf_df["CHROM"].astype(str).str.replace('^chr', '', regex=True)
        gwas_df["CHR_ID"] = gwas_df["CHR_ID"].astype(str)

        print("Step 4: Ensuring position columns are strings...")
        vcf_df["POS"] = vcf_df["POS"].astype(str)
        gwas_df["CHR_POS"] = gwas_df["CHR_POS"].astype(str)

        print("Step 5: Merging VCF and GWAS data...")
        annotated_df = pd.merge(vcf_df, gwas_df, left_on=["CHROM", "POS"], right_on=["CHR_ID", "CHR_POS"], how="left")

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
            'Icon': 'first'
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
            details.append((row['DISEASE/TRAIT'], row['REGION'], row['SNPS'], row['MAPPED_GENE'], row['Groups of Disease/Trait']))
            embedded_svgs.append(svg_string)

        summary_text = f"This report includes {len(data['DISEASE/TRAIT'].unique())} unique diseases/traits analyzed."

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
                    h1 {
                        text-align: center;
                        background: linear-gradient(to right, #3c79aa, #FF6B6B);
                        color: white;
                        padding: 20px;
                        border-radius: 8px;
                        box-shadow: 0 4px 10px rgba(0, 0, 0, 0.2);
                    }
                    h2 {
                        margin-top: 40px;
                        color: #FF003B;
                        border-bottom: 3px solid #FF003B;
                        padding-bottom: 5px;
                        display: inline-block;
                    }
                    p {
                        margin: 10px 0;
                        font-size: 16px;
                        color: rgb(87, 107, 145);
                    }
                    b {
                        color: #333;
                    }
                    .chart-container {
                        text-align: center;
                        background: white;
                        padding: 20px;
                        border-radius: 8px;
                        box-shadow: 0 4px 10px rgba(0, 0, 0, 0.1);
                        margin: 20px 0;
                        width: 100%;  /* Adjust based on the window size */
                        max-width: 95%; /* Set a max-width to prevent it from getting too big */
                    }
                    
                    .chart svg {
                        width: 95%;  /* Make SVG scale with container */
                        height: auto; /* Maintain aspect ratio */
                    }
                    .download-btn {
                        display: inline-block;
                        background: #FF003B;
                        color: white;
                        padding: 8px 12px;
                        font-size: 14px;
                        border-radius: 5px;
                        border: none;
                        cursor: pointer;
                        transition: 0.3s;
                    }
                    .download-btn:hover {
                        background: #E60033;
                    }
                    .logo-container {
                        text-align: center;
                    }
                    hr {
                        border: none;
                        height: 2px;
                        background: #ddd;
                        margin: 40px 0;

                </style>
                                <script>
                    function resizeChart() {
                        var chart = document.getElementById('chart_7');
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
                <div class="logo-container">
                    <div class="logo" id="logo_01">
                        <?xml version="1.0" encoding="UTF-8" standalone="no"?>
                        <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
                        <svg width="50%" height="50%" viewBox="0 0 1000 1000" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve" xmlns:serif="http://www.serif.com/" style="fill-rule:evenodd;clip-rule:evenodd;stroke-linejoin:round;stroke-miterlimit:2;">
                            <g id="Layer-1" serif:id="Layer 1">
                                <g transform="matrix(1,0,0,1,225.698,724.339)">
                                    <path d="M0,-29.194L3.462,-29.194L3.462,-13.637L17.058,-29.194L21.395,-29.194L10.593,-16.807L22.688,0L18.685,0L8.383,-14.263L3.462,-8.632L3.462,0L0,0L0,-29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,256.102,695.145)">
                                    <path d="M0,29.194L0,0L3.462,0L3.462,12.845L20.853,12.845L20.853,0L24.315,0L24.315,29.194L20.853,29.194L20.853,15.973L3.462,15.973L3.462,29.194L0,29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,295.973,701.443)">
                                    <path d="M0,16.598C2.029,18.795 4.601,19.893 7.716,19.893C10.829,19.893 13.401,18.802 15.432,16.619C17.461,14.437 18.476,11.663 18.476,8.298C18.476,4.935 17.461,2.161 15.432,-0.022C13.401,-2.205 10.829,-3.296 7.716,-3.296C4.601,-3.296 2.029,-2.212 0,-0.043C-2.03,2.154 -3.045,4.935 -3.045,8.298C-3.045,11.663 -2.03,14.429 0,16.598M-2.461,-2.379C0.208,-5.186 3.6,-6.59 7.716,-6.59C11.83,-6.59 15.229,-5.186 17.913,-2.379C20.596,0.43 21.938,3.989 21.938,8.298C21.938,12.609 20.603,16.167 17.934,18.975C15.237,21.785 11.83,23.188 7.716,23.188C3.6,23.188 0.201,21.785 -2.482,18.975C-5.165,16.167 -6.506,12.609 -6.506,8.298C-6.506,3.989 -5.158,0.43 -2.461,-2.379" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,326.961,724.339)">
                                    <path d="M0,-29.194L3.462,-29.194L20.853,-5.547L20.853,-29.194L24.315,-29.194L24.315,0L20.853,0L3.462,-23.647L3.462,0L0,0L0,-29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,375.965,724.339)">
                                    <path d="M0,-29.194L3.462,-29.194L3.462,-13.637L17.058,-29.194L21.396,-29.194L10.594,-16.807L22.688,0L18.685,0L8.383,-14.263L3.462,-8.632L3.462,0L0,0L0,-29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,423.469,705.696)">
                                    <path d="M0,8.091L-5.964,-6.632L-11.928,8.091L0,8.091ZM-7.924,-10.552L-4.004,-10.552L7.799,18.643L4.254,18.643L1.251,11.219L-13.179,11.219L-16.182,18.643L-19.727,18.643L-7.924,-10.552Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,439.026,695.145)">
                                    <path d="M0,29.194L0,0L20.52,0L20.52,3.128L3.462,3.128L3.462,12.845L18.851,12.845L18.851,15.973L3.462,15.973L3.462,26.066L20.52,26.066L20.52,29.194L0,29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,469.054,724.339)">
                                    <path d="M0,-29.194L3.462,-29.194L20.853,-5.547L20.853,-29.194L24.315,-29.194L24.315,0L20.853,0L3.462,-23.647L3.462,0L0,0L0,-29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,521.311,724.631)">
                                    <path d="M0,-29.486L0,-10.384C0,-8.327 0.806,-6.632 2.419,-5.297C4.031,-3.962 6.089,-3.295 8.592,-3.295C11.066,-3.295 13.117,-3.962 14.743,-5.297C16.37,-6.632 17.183,-8.327 17.183,-10.384L17.183,-29.486L20.645,-29.486L20.645,-10.426C20.645,-7.395 19.504,-4.9 17.225,-2.94C14.945,-0.98 12.067,0 8.592,0C5.088,0 2.203,-0.98 -0.062,-2.94C-2.329,-4.9 -3.462,-7.395 -3.462,-10.426L-3.462,-29.486L0,-29.486Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,552.132,724.339)">
                                    <path d="M0,-29.194L3.462,-29.194L20.853,-5.547L20.853,-29.194L24.315,-29.194L24.315,0L20.853,0L3.462,-23.647L3.462,0L0,0L0,-29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,0,419.484)">
                                    <rect x="586.832" y="275.661" width="3.462" height="29.194" style="fill:rgb(168,59,36);"/>
                                </g>
                                <g transform="matrix(1,0,0,1,598.092,724.339)">
                                    <path d="M0,-29.194L3.545,-29.194L13.596,-3.962L23.606,-29.194L27.151,-29.194L15.557,0L11.594,0L0,-29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,633.042,695.145)">
                                    <path d="M0,29.194L0,0L20.52,0L20.52,3.128L3.462,3.128L3.462,12.845L18.851,12.845L18.851,15.973L3.462,15.973L3.462,26.066L20.52,26.066L20.52,29.194L0,29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,675.708,721.211)">
                                    <path d="M0,-22.938L-9.175,-22.938L-9.175,-12.261L0,-12.261C1.501,-12.261 2.738,-12.768 3.712,-13.784C4.685,-14.798 5.172,-16.07 5.172,-17.6C5.172,-19.157 4.685,-20.436 3.712,-21.437C2.738,-22.438 1.501,-22.938 0,-22.938M-12.637,3.128L-12.637,-26.066L0,-26.066C2.502,-26.066 4.56,-25.274 6.173,-23.689C7.785,-22.104 8.592,-20.074 8.592,-17.6C8.592,-15.57 8.021,-13.818 6.882,-12.345C5.741,-10.871 4.226,-9.898 2.336,-9.425L9.384,3.128L5.672,3.128L-1.251,-9.175L-9.175,-9.175L-9.175,3.128L-12.637,3.128Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,692.056,698.523)">
                                    <path d="M0,22.438L2.419,20.019C4.338,21.881 6.854,22.813 9.968,22.813C12.025,22.813 13.673,22.389 14.91,21.541C16.147,20.693 16.766,19.56 16.766,18.142C16.766,17.03 16.363,16.091 15.557,15.327C14.75,14.562 13.749,13.986 12.554,13.596C11.358,13.207 10.058,12.755 8.654,12.241C7.25,11.726 5.95,11.198 4.755,10.656C3.559,10.114 2.558,9.279 1.752,8.153C0.945,7.027 0.542,5.671 0.542,4.086C0.542,1.696 1.39,-0.195 3.086,-1.585C4.782,-2.975 7.062,-3.67 9.926,-3.67C11.872,-3.67 13.652,-3.343 15.265,-2.691C16.877,-2.037 18.128,-1.14 19.018,0L16.641,2.418C16.029,1.557 15.119,0.876 13.909,0.375C12.7,-0.125 11.33,-0.375 9.801,-0.375C7.966,-0.375 6.506,0 5.422,0.751C4.338,1.501 3.795,2.543 3.795,3.878C3.795,5.129 4.351,6.159 5.464,6.964C6.575,7.771 7.924,8.417 9.509,8.904C11.094,9.391 12.685,9.925 14.285,10.509C15.883,11.093 17.239,12.004 18.351,13.242C19.463,14.479 20.019,16.042 20.019,17.933C20.019,20.463 19.102,22.459 17.267,23.918C15.432,25.378 12.943,26.108 9.801,26.108C5.658,26.108 2.391,24.884 0,22.438" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,0,419.484)">
                                    <rect x="720.792" y="275.661" width="3.462" height="29.194" style="fill:rgb(168,59,36);"/>
                                </g>
                                <g transform="matrix(1,0,0,1,732.053,724.339)">
                                    <path d="M0,-29.194L24.523,-29.194L24.523,-26.066L13.972,-26.066L13.972,0L10.51,0L10.51,-26.066L0,-26.066L0,-29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,761.789,724.339)">
                                    <path d="M0,-29.194L3.879,-29.194L12.846,-14.972L21.813,-29.194L25.691,-29.194L14.597,-11.51L14.597,0L11.136,0L11.136,-11.552L0,-29.194Z" style="fill:rgb(168,59,36);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,85.3204,627.655)">
                                    <path d="M0,38.848L-25.564,15.731L-25.564,38.848L-33.004,38.848L-33.004,0L-32.207,0L-6.645,23.49L-6.645,0L0.795,0L0.795,38.848L0,38.848Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,105.619,640.25)">
                                    <path d="M0,13.658L11.107,13.658L5.58,1.913L0,13.658ZM-14.668,26.253L5.208,-12.595L5.952,-12.595L25.828,26.253L17.75,26.253L14.721,20.141L-3.614,20.141L-6.59,26.253L-14.668,26.253Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,164.235,659.913)">
                                    <path d="M0,-25.722L-11.107,-25.722L-11.107,6.59L-18.547,6.59L-18.547,-25.722L-29.707,-25.722L-29.707,-32.312L0,-32.312L0,-25.722Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(-1,0,0,1,348.772,294.158)">
                                    <rect x="170.666" y="333.497" width="7.44" height="38.848" style="fill:rgb(43,41,74);"/>
                                </g>
                                <g transform="matrix(1,0,0,1,206.138,660.287)">
                                    <path d="M0,-26.467C-1.768,-26.467 -3.421,-26.103 -4.959,-25.377C-6.497,-24.651 -7.84,-23.685 -8.99,-22.481C-10.139,-21.276 -11.049,-19.867 -11.721,-18.256C-12.392,-16.644 -12.728,-14.969 -12.728,-13.234C-12.728,-11.462 -12.392,-9.78 -11.721,-8.186C-11.049,-6.591 -10.139,-5.183 -8.99,-3.96C-7.84,-2.738 -6.497,-1.763 -4.959,-1.037C-3.421,-0.311 -1.768,0.052 0,0.052C1.767,0.052 3.42,-0.311 4.958,-1.037C6.496,-1.763 7.839,-2.738 8.989,-3.96C10.138,-5.183 11.048,-6.591 11.72,-8.186C12.392,-9.78 12.728,-11.462 12.728,-13.234C12.728,-14.969 12.392,-16.644 11.72,-18.256C11.048,-19.867 10.138,-21.276 8.989,-22.481C7.839,-23.685 6.496,-24.651 4.958,-25.377C3.42,-26.103 1.767,-26.467 0,-26.467M-0.027,6.961C-2.79,6.961 -5.394,6.429 -7.839,5.366C-10.283,4.304 -12.418,2.86 -14.242,1.035C-16.068,-0.789 -17.511,-2.924 -18.574,-5.368C-19.636,-7.814 -20.168,-10.435 -20.168,-13.234C-20.168,-15.997 -19.636,-18.61 -18.574,-21.072C-17.511,-23.534 -16.068,-25.678 -14.242,-27.503C-12.418,-29.327 -10.283,-30.772 -7.839,-31.834C-5.394,-32.897 -2.79,-33.428 -0.027,-33.428C2.772,-33.428 5.394,-32.897 7.839,-31.834C10.283,-30.772 12.417,-29.327 14.242,-27.503C16.067,-25.678 17.511,-23.534 18.574,-21.072C19.636,-18.61 20.168,-15.997 20.168,-13.234C20.168,-10.435 19.636,-7.814 18.574,-5.368C17.511,-2.924 16.067,-0.789 14.242,1.035C12.417,2.86 10.283,4.304 7.839,5.366C5.394,6.429 2.772,6.961 -0.027,6.961" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,266.165,627.655)">
                                    <path d="M0,38.848L-25.564,15.731L-25.564,38.848L-33.004,38.848L-33.004,0L-32.207,0L-6.645,23.49L-6.645,0L0.795,0L0.795,38.848L0,38.848Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,286.464,640.25)">
                                    <path d="M0,13.658L11.107,13.658L5.58,1.913L0,13.658ZM-14.668,26.253L5.208,-12.595L5.952,-12.595L25.828,26.253L17.75,26.253L14.721,20.141L-3.614,20.141L-6.59,26.253L-14.668,26.253Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,345.133,634.191)">
                                    <path d="M0,25.722L0,32.312L-28.007,32.312L-28.007,-6.59L-20.567,-6.59L-20.567,25.722L0,25.722Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,385.841,645.883)">
                                    <path d="M0,2.339C1.098,2.339 2.054,2.152 2.87,1.781C3.684,1.409 4.349,0.904 4.863,0.267C5.376,-0.371 5.757,-1.098 6.005,-1.913C6.253,-2.727 6.377,-3.577 6.377,-4.464C6.377,-5.348 6.218,-6.225 5.899,-7.094C5.58,-7.962 5.119,-8.732 4.517,-9.406C3.914,-10.078 3.17,-10.628 2.285,-11.054C1.399,-11.478 0.372,-11.691 -0.797,-11.691L-9.885,-11.691L-9.885,2.339L0,2.339ZM0,-18.281C1.913,-18.281 3.711,-17.917 5.394,-17.191C7.076,-16.465 8.538,-15.481 9.778,-14.242C11.018,-13.001 12.001,-11.54 12.728,-9.857C13.454,-8.175 13.817,-6.377 13.817,-4.464C13.817,-2.551 13.498,-0.778 12.861,0.851C12.223,2.481 11.284,3.906 10.044,5.129C8.804,6.352 7.236,7.308 5.341,7.999C3.445,8.69 1.257,9.035 -1.222,9.035L-9.885,9.035L-9.885,20.62L-17.325,20.62L-17.325,-18.281L0,-18.281Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,438.877,666.503)">
                                    <path d="M0,-38.901L0,0L-7.44,0L-7.44,-15.464L-26.625,-15.464L-26.625,0L-34.065,0L-34.065,-38.901L-26.625,-38.901L-26.625,-22.054L-7.44,-22.054L-7.44,-38.901L0,-38.901Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,455.032,659.913)">
                                    <path d="M0,-25.722L0,-15.996L18.813,-15.996L18.813,-9.247L0,-9.247L0,-0.107L21.204,-0.107L21.204,6.59L-7.44,6.59L-7.44,-32.312L21.204,-32.312L21.204,-25.722L0,-25.722Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,517.158,627.655)">
                                    <path d="M0,38.848L-25.564,15.731L-25.564,38.848L-33.004,38.848L-33.004,0L-32.207,0L-6.645,23.49L-6.645,0L0.796,0L0.796,38.848L0,38.848Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,544.924,660.287)">
                                    <path d="M0,-26.467C-1.768,-26.467 -3.421,-26.103 -4.959,-25.377C-6.497,-24.651 -7.84,-23.685 -8.99,-22.481C-10.139,-21.276 -11.049,-19.867 -11.721,-18.256C-12.392,-16.644 -12.728,-14.969 -12.728,-13.234C-12.728,-11.462 -12.392,-9.78 -11.721,-8.186C-11.049,-6.591 -10.139,-5.183 -8.99,-3.96C-7.84,-2.738 -6.497,-1.763 -4.959,-1.037C-3.421,-0.311 -1.768,0.052 0,0.052C1.767,0.052 3.42,-0.311 4.958,-1.037C6.496,-1.763 7.839,-2.738 8.989,-3.96C10.138,-5.183 11.048,-6.591 11.72,-8.186C12.391,-9.78 12.728,-11.462 12.728,-13.234C12.728,-14.969 12.391,-16.644 11.72,-18.256C11.048,-19.867 10.138,-21.276 8.989,-22.481C7.839,-23.685 6.496,-24.651 4.958,-25.377C3.42,-26.103 1.767,-26.467 0,-26.467M-0.027,6.961C-2.79,6.961 -5.394,6.429 -7.839,5.366C-10.283,4.304 -12.418,2.86 -14.242,1.035C-16.068,-0.789 -17.511,-2.924 -18.574,-5.368C-19.636,-7.814 -20.168,-10.435 -20.168,-13.234C-20.168,-15.997 -19.636,-18.61 -18.574,-21.072C-17.511,-23.534 -16.068,-25.678 -14.242,-27.503C-12.418,-29.327 -10.283,-30.772 -7.839,-31.834C-5.394,-32.897 -2.79,-33.428 -0.027,-33.428C2.772,-33.428 5.394,-32.897 7.839,-31.834C10.283,-30.772 12.417,-29.327 14.242,-27.503C16.067,-25.678 17.511,-23.534 18.574,-21.072C19.636,-18.61 20.168,-15.997 20.168,-13.234C20.168,-10.435 19.636,-7.814 18.574,-5.368C17.511,-2.924 16.067,-0.789 14.242,1.035C12.417,2.86 10.283,4.304 7.839,5.366C5.394,6.429 2.772,6.961 -0.027,6.961" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,611.113,666.503)">
                                    <path d="M0,-38.848L0,0L-7.44,0L-7.44,-19.663L-19.61,-5.049L-31.727,-19.663L-31.727,0L-39.167,0L-39.167,-38.848L-38.688,-38.848L-19.61,-16.421L-0.478,-38.848L0,-38.848Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,627.375,659.913)">
                                    <path d="M0,-25.722L0,-15.996L18.813,-15.996L18.813,-9.247L0,-9.247L0,-0.107L21.204,-0.107L21.204,6.59L-7.44,6.59L-7.44,-32.312L21.204,-32.312L21.204,-25.722L0,-25.722Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(-1,0,0,1,1359.44,294.158)">
                                    <rect x="676.001" y="333.497" width="7.44" height="38.848" style="fill:rgb(43,41,74);"/>
                                </g>
                                <g transform="matrix(1,0,0,1,726.276,627.655)">
                                    <path d="M0,38.848L-25.564,15.731L-25.564,38.848L-33.004,38.848L-33.004,0L-32.207,0L-6.645,23.49L-6.645,0L0.796,0L0.796,38.848L0,38.848Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,734.298,632.864)">
                                    <path d="M0,28.538L4.039,22.586C5.987,23.933 7.998,25.093 10.071,26.066C12.143,27.042 14.313,27.528 16.581,27.528C17.785,27.528 18.848,27.378 19.769,27.076C20.69,26.775 21.47,26.386 22.108,25.907C22.745,25.429 23.224,24.88 23.542,24.26C23.861,23.641 24.021,23.029 24.021,22.427C24.021,21.647 23.782,20.983 23.303,20.434C22.825,19.885 22.187,19.433 21.39,19.078C20.593,18.725 19.654,18.423 18.574,18.175C17.492,17.928 16.368,17.696 15.199,17.484C14.419,17.343 13.516,17.175 12.489,16.979C11.461,16.785 10.398,16.519 9.3,16.182C8.202,15.846 7.121,15.403 6.058,14.853C4.995,14.305 4.039,13.613 3.189,12.781C2.338,11.949 1.656,10.957 1.143,9.805C0.629,8.654 0.372,7.299 0.372,5.739C0.372,3.473 0.832,1.586 1.754,0.08C2.675,-1.426 3.853,-2.612 5.288,-3.482C6.723,-4.349 8.317,-4.96 10.071,-5.315C11.824,-5.668 13.516,-5.846 15.146,-5.846C16.775,-5.846 18.228,-5.739 19.504,-5.527C20.779,-5.315 21.966,-5.004 23.064,-4.597C24.162,-4.19 25.234,-3.676 26.279,-3.056C27.324,-2.436 28.449,-1.701 29.654,-0.851L25.509,4.995C23.666,3.72 21.895,2.746 20.194,2.072C18.494,1.399 16.704,1.062 14.827,1.062C14.189,1.062 13.463,1.142 12.648,1.302C11.833,1.461 11.062,1.719 10.336,2.072C9.61,2.427 9.008,2.879 8.53,3.428C8.051,3.977 7.812,4.642 7.812,5.421C7.812,6.2 8.095,6.855 8.662,7.387C9.229,7.918 9.955,8.361 10.841,8.716C11.726,9.07 12.701,9.371 13.764,9.619C14.827,9.867 15.854,10.115 16.846,10.363C17.59,10.541 18.467,10.718 19.477,10.894C20.487,11.072 21.531,11.319 22.613,11.639C23.693,11.957 24.765,12.382 25.828,12.914C26.891,13.445 27.838,14.128 28.671,14.96C29.503,15.793 30.176,16.819 30.69,18.042C31.203,19.265 31.461,20.744 31.461,22.479C31.461,24.499 31.106,26.253 30.398,27.741C29.689,29.228 28.679,30.46 27.369,31.434C26.058,32.409 24.499,33.143 22.692,33.64C20.885,34.136 18.883,34.384 16.687,34.384C14.668,34.384 12.887,34.232 11.346,33.932C9.805,33.631 8.397,33.224 7.121,32.71C5.846,32.196 4.641,31.585 3.507,30.876C2.373,30.168 1.204,29.388 0,28.538" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,800.196,659.913)">
                                    <path d="M0,-25.722L-11.107,-25.722L-11.107,6.59L-18.547,6.59L-18.547,-25.722L-29.707,-25.722L-29.707,-32.312L0,-32.312L0,-25.722Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(-1,0,0,1,1620.69,294.158)">
                                    <rect x="806.626" y="333.497" width="7.44" height="38.848" style="fill:rgb(43,41,74);"/>
                                </g>
                                <g transform="matrix(1,0,0,1,851.85,659.913)">
                                    <path d="M0,-25.722L-11.107,-25.722L-11.107,6.59L-18.547,6.59L-18.547,-25.722L-29.707,-25.722L-29.707,-32.312L0,-32.312L0,-25.722Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,873.905,627.656)">
                                    <path d="M0,39.538C-2.516,39.538 -4.801,39.148 -6.855,38.369C-8.911,37.591 -10.664,36.492 -12.117,35.074C-13.57,33.658 -14.694,31.939 -15.491,29.92C-16.288,27.9 -16.687,25.65 -16.687,23.17L-16.687,0L-9.247,0L-9.247,23.133C-9.247,25.078 -8.928,26.687 -8.29,27.962C-7.653,29.235 -6.865,30.243 -5.925,30.985C-4.987,31.728 -3.977,32.242 -2.896,32.524C-1.816,32.808 -0.85,32.949 0,32.949C0.85,32.949 1.815,32.808 2.896,32.524C3.977,32.242 4.986,31.728 5.926,30.985C6.864,30.243 7.653,29.235 8.29,27.962C8.928,26.687 9.247,25.078 9.247,23.133L9.247,0L16.687,0L16.687,23.17C16.687,25.65 16.297,27.9 15.518,29.92C14.738,31.939 13.622,33.658 12.17,35.074C10.717,36.492 8.963,37.591 6.909,38.369C4.854,39.148 2.551,39.538 0,39.538" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,927.313,659.913)">
                                    <path d="M0,-25.722L-11.107,-25.722L-11.107,6.59L-18.547,6.59L-18.547,-25.722L-29.707,-25.722L-29.707,-32.312L0,-32.312L0,-25.722Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,940.121,659.913)">
                                    <path d="M0,-25.722L0,-15.996L18.813,-15.996L18.813,-9.247L0,-9.247L0,-0.107L21.204,-0.107L21.204,6.59L-7.44,6.59L-7.44,-32.312L21.204,-32.312L21.204,-25.722L0,-25.722Z" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,587.706,322.434)">
                                    <path d="M0,54.721C24.541,76.068 21.572,113.058 10.957,140.625C-4.099,105.428 -41.911,92.471 -75.608,82.15C-123.401,64.508 -176.251,58.687 -219.392,29.903C-238.017,18.213 -252.13,-1.061 -255.921,-22.86C-259.422,-47.349 -247.503,-72.281 -227.001,-85.904C-227.634,-62.576 -221.199,-37.262 -200.744,-23.509C-141.058,19.366 -60.692,13.757 0,54.721" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,404.513,531.736)">
                                    <path d="M0,-145.266C47.757,-102.333 119.749,-105.946 168.518,-64.815C186.123,-47.33 183.041,-19.657 171.311,0C163.333,-9.626 155.218,-19.903 143.415,-25.027C105.831,-43.392 63.286,-49.765 26.769,-70.636C9.969,-80.372 -6.539,-96.001 -7.228,-116.674C-8.303,-126.91 -3.225,-135.99 0,-145.266" style="fill:rgb(244,245,248);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(1,0,0,1,541.583,555.23)">
                                    <path d="M0,-24.468C17.116,-16.874 17.996,13.968 -0.827,20.361C-19.668,27.829 -40.621,22.105 -57.308,11.969C-72.893,3.798 -80.888,-3.949 -88.677,-23.494C-92.765,-30.188 -91.589,-45.199 -85.732,-50.783C-77.819,-52.297 -70.113,-47.966 -62.583,-45.953C-41.969,-38.128 -20.695,-32.135 0,-24.468" style="fill:rgb(165,26,49);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(0.0331851,-0.999449,-0.999449,-0.0331851,521.16,268.315)">
                                    <path d="M-24.757,-23.948C-38.436,-23.947 -49.527,-12.858 -49.527,0.821C-49.527,14.502 -38.437,25.593 -24.757,25.592C-11.076,25.593 0.014,14.503 0.014,0.823C0.014,-12.858 -11.076,-23.949 -24.757,-23.948" style="fill:rgb(165,26,49);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(0.033064,-0.999453,-0.999453,-0.033064,563.382,229.491)">
                                    <path d="M-18.107,-17.516C-28.115,-17.514 -36.224,-9.404 -36.223,0.6C-36.221,10.609 -28.11,18.718 -18.107,18.716C-8.097,18.714 0.012,10.603 0.01,0.6C0.008,-9.409 -8.103,-17.519 -18.107,-17.516" style="fill:rgb(43,41,74);fill-rule:nonzero;"/>
                                </g>
                                <g transform="matrix(0.0331865,-0.999449,-0.999449,-0.0331865,553.828,329.601)">
                                    <path d="M-8.947,-8.655C-13.889,-8.653 -17.896,-4.645 -17.896,0.293C-17.896,5.237 -13.889,9.246 -8.947,9.244C-4.005,9.246 0.002,5.239 0.001,0.295C0.001,-4.648 -4.004,-8.652 -8.947,-8.655" style="fill:rgb(244,245,248);fill-rule:nonzero;"/>
                                </g>
                            </g>
                        </svg>

                    </div>
                </div>
                    <div>
                    <h1>Whole Genome Analysis Report</h1>
                    <p>This report is for research use only</p>
                    <p>{{ summary }}</p>
                    <h2>Disclaimer</h2>
                        <p>The service provided by Khon Kaen University National
                            Phenome Institute (KKUNPhI) are currently for research use
                            only. As they have not been submitted for review to any
                            regulatory agency/notified body for clinical diagnostics, caution
                            must be excercisrd when describing the application of KKUNPhI
                            service in the clinical research field.
                        </p>
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
                            <p><b>Alleles:</b> An allele is a variant form of a gene found at a specific position (locus)
                                on a chromosome. Each allele is inherited, one from each parent.
                            </p>
                            <p><b>Mapped Gene:</b> Genes mapped near or overlapping the SNP.
                            </p>
                            <p><b>Chromosomal region:</b> The genomic region associated with the trait or disease.
                            </p>
                            <p><b>Risk Allele Frequency (%):</b> The frequency of the risk allele in the population.
                            </p>

                    </div>
                {% for (title_, region_, snps_, mapped_gene_, group_trait_), svg_ in data_source %}
                    <h2>{{ title_ }}</h2>
                    <p><b>Region:</b> {{ region_ }}</p>
                    <p><b>SNPs:</b> {{ snps_ }}</p>
                    <p><b>Mapped Gene:</b> {{ mapped_gene_ }}</p>
                    <p><b>Group of disease/trait:</b> {{ group_trait_ }}</p>
                    
                    <div class="chart-container">
                        <div class="chart" id="chart_{{ loop.index }}">
                            {{ svg_ | safe }}
                        </div>
                        <button class="download-btn" onclick="downloadChart('{{ loop.index }}', '{{ title_ }}')">Download PNG</button>
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
        rendered_html = template.render(data_source=zip(details, embedded_svgs), summary=summary_text)

        with open(output_path, 'w') as f:
            f.write(rendered_html)

        print(f"Report saved to {output_path}")

    def generate_report(self):
        self.prepare_report_data()
        self.generate_html_report()
        print(f"Report saved to {os.path.join(self.report_path, 'GWAS_report.html')}")

if __name__ == '__main__':
    vcf_file = "VIP01/Barcode10/medaka-variant-out/medaka.annotated.vcf"
    gwas_file = "VIP01/gwas_catalog_grouped.csv"
    output_file = "VIP01/test"

    mapper = MapGWASSNPs(vcf_file, gwas_file, output_file, cut_off_qual=0, filt_nr_disease=True)
    mapper.map_snps()
    mapper.generate_report()