import os
import re
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from jinja2 import Template

class MapGWASSNPs:
    def __init__(self, vcf_file_path: str, gwas_file_path: str, output_file_path: str,
                 cut_off_qual: int = 20, filt_nr_disease: bool = True):
        self.vcf_file = vcf_file_path
        self.gwas_file = gwas_file_path
        self.output_root = output_file_path.replace('\\', '/').rstrip('/')

        # Output folders
        self.report_path = os.path.join(self.output_root, "report")
        self.report_data_path = os.path.join(self.report_path, "data")
        os.makedirs(self.report_data_path, exist_ok=True)

        self.cut_off_qual = cut_off_qual
        self.filt_nr_disease = filt_nr_disease

        # Will be filled later
        self.vcf_report = None
        self.annotated_df = None
        self.report_data = None

    # ---------- helpers ----------
    @staticmethod
    def _classify_variant(ref: str, alt: str) -> str:
        # ALT can be comma-separated (multi-allelic) – choose first for type check
        alt_first = str(alt).split(',')[0]
        rl, al = len(str(ref)), len(alt_first)
        if rl == 1 and al == 1:
            return "SNPs"
        elif rl < al:
            return "INS"
        elif rl > al:
            return "DEL"
        else:
            return "COMPLEX"

    @staticmethod
    def _to_numeric_safe(s: pd.Series) -> pd.Series:
        """Coerce to numeric; handle weird scientific formats like '1 x 10-4'."""
        # Normalize "a x 10^b" or "a x 10-b" into "aE b"
        cleaned = s.astype(str).str.replace(r'×', 'x', regex=False)
        cleaned = cleaned.str.replace(
            r'^\s*([+-]?\d*\.?\d+)\s*[xX]\s*10\s*[\^]?\s*([+-]?\d+)\s*$',
            lambda m: f"{m.group(1)}e{m.group(2)}",
            regex=True
        )
        # Remove commas and stray spaces, turn NR/NA/– to NaN
        cleaned = cleaned.str.replace(',', '', regex=False).str.strip()
        cleaned = cleaned.replace(
            {r'^(NR|NA|N/?A|None|nan|—|-|\.?)$': np.nan}, regex=True
        )
        return pd.to_numeric(cleaned, errors='coerce')

    # ---------- pipeline ----------
    def map_snps(self):
        print("Step 1: Reading VCF file...")
        vcf_columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
        if self.vcf_file.endswith('.gz'):
            vcf_df = pd.read_csv(self.vcf_file, compression='gzip', sep='\t', comment="#", names=vcf_columns)
        else:
            vcf_df = pd.read_csv(self.vcf_file, comment="#", sep='\t', names=vcf_columns)

        # Normalize basic types
        vcf_df["CHROM"] = vcf_df["CHROM"].astype(str)
        vcf_df["POS"] = vcf_df["POS"].astype(str)
        vcf_df["QUAL"] = self._to_numeric_safe(vcf_df["QUAL"])

        # Filter QUAL >= cutoff ONLY ONCE (affects both merge and stats)
        vcf_df = vcf_df.loc[vcf_df["QUAL"] >= float(self.cut_off_qual)].copy()
        print(f'Number of variants PASS at Quality ≥ {self.cut_off_qual}: {vcf_df.shape[0]:,}')

        # Variant type
        vcf_df["TYPE"] = vcf_df.apply(lambda r: self._classify_variant(r["REF"], r["ALT"]), axis=1)

        self.vcf_report = vcf_df.copy()
        print("VCF file FILTER=='PASS' count (after QUAL filter):", vcf_df[vcf_df["FILTER"] == "PASS"].shape[0])

        print("Step 2: Reading GWAS catalog...")
        if self.gwas_file.endswith('.gz'):
            gwas_df = pd.read_csv(self.gwas_file, low_memory=False, compression='gzip')
        else:
            gwas_df = pd.read_csv(self.gwas_file, low_memory=False)
        print("GWAS catalog shape:", gwas_df.shape)

        print("Step 3: Normalizing identifiers...")
        gwas_df["CHR_ID"] = gwas_df["CHR_ID"].astype(str)
        gwas_df["CHR_POS"] = gwas_df["CHR_POS"].astype(str)
        print("Identifier normalization PASS")

        print("Step 4: Merge on chromosome/position...")
        annotated_df = pd.merge(
            vcf_df, gwas_df,
            left_on=["CHROM", "POS"],
            right_on=["CHR_ID", "CHR_POS"],
            how="inner"
        )
        print("Merge PASS; rows:", annotated_df.shape[0])

        # Keep essential + clean numerics before filters/agg
        if "DISEASE/TRAIT" in annotated_df.columns:
            annotated_df.dropna(subset=['DISEASE/TRAIT'], inplace=True)
        else:
            raise KeyError("Column 'DISEASE/TRAIT' not found in GWAS file.")

        # Optional: filter out "NR" traits
        if self.filt_nr_disease:
            annotated_df = annotated_df[annotated_df["DISEASE/TRAIT"].astype(str) != "NR"]

        # Clean numeric GWAS columns used later
        for col in ["RISK ALLELE FREQUENCY", "P-VALUE"]:
            if col in annotated_df.columns:
                annotated_df[col] = self._to_numeric_safe(annotated_df[col])

        # Persist CSV
        out_csv = os.path.join(self.report_data_path, 'in-house_report.csv')
        print("Saving annotated data to CSV...")
        annotated_df.to_csv(out_csv, index=False)
        print(f"Annotated data saved to {out_csv}")

        self.annotated_df = annotated_df
        return annotated_df

    def prepare_report_data(self):
        if self.annotated_df is None:
            raise RuntimeError("annotated_df is empty. Run map_snps() first.")

        df = self.annotated_df.copy()

        # Build an order key to pick "representative" rows per trait (lowest p-value, then highest RAF)
        df["P_SORT"] = self._to_numeric_safe(df["P-VALUE"]) if "P-VALUE" in df.columns else np.nan
        df["RAF_SORT"] = self._to_numeric_safe(df["RISK ALLELE FREQUENCY"]) if "RISK ALLELE FREQUENCY" in df.columns else np.nan
        df.sort_values(by=["P_SORT", "RAF_SORT"], ascending=[True, False], inplace=True)

        # For each trait, take the first row after sorting
        keep_cols = [
            'DISEASE/TRAIT', 'CHR_ID', 'CHR_POS', 'TYPE',
            'RISK ALLELE FREQUENCY', 'P-VALUE',
            'REGION', 'SNPS', 'MAPPED_GENE',
            'Groups of Disease/Trait', 'MAPPED_TRAIT_URI', 'MAPPED_TRAIT_DESCRIPTION'
        ]
        keep_cols = [c for c in keep_cols if c in df.columns]
        rep = df[keep_cols].drop_duplicates(subset=['DISEASE/TRAIT']).copy()

        # Compute standardized percentage column
        rep['RAF (%)'] = self._to_numeric_safe(rep['RISK ALLELE FREQUENCY']) * 100.0
        rep.dropna(subset=['RAF (%)'], inplace=True)
        rep.sort_values(by='RAF (%)', ascending=False, inplace=True)

        # Save
        out_csv = os.path.join(self.report_data_path, 'report_data.csv')
        print("Saving report data to CSV...")
        rep.to_csv(out_csv, index=False)

        self.report_data = rep
        return rep

    def generate_html_report(self):
        if self.report_data is None:
            raise RuntimeError("report_data is empty. Run prepare_report_data() first.")

        data = self.report_data
        output_path = os.path.join(self.report_path, 'GWAS_report.html')

        # ---------- Sunburst (built once) ----------
        sun_cols_all = ["TYPE", "Groups of Disease/Trait", "CHR_ID", "REGION", "SNPS", "DISEASE/TRAIT"]
        sun_cols = [c for c in sun_cols_all if c in data.columns]
        df_sun = data.copy()

        # Treat empty strings as missing, then DROP rows with missing ancestors (TYPE, group, CHR_ID, REGION)
        df_sun[sun_cols] = df_sun[sun_cols].replace("", np.nan)
        required_ancestors = [c for c in ["TYPE", "Groups of Disease/Trait", "CHR_ID", "REGION"] if c in df_sun.columns]
        if required_ancestors:
            df_sun = df_sun.dropna(subset=required_ancestors)

        # Ensure color column numeric
        color_col = "RAF (%)"
        if color_col in df_sun.columns:
            df_sun[color_col] = pd.to_numeric(df_sun[color_col], errors="coerce")
        else:
            df_sun[color_col] = np.nan

        if df_sun.empty or not sun_cols:
            fig_sun = px.sunburst(pd.DataFrame({c: [] for c in sun_cols_all if c in data.columns}),
                                  path=[c for c in sun_cols_all if c in data.columns])
        else:
            fig_sun = px.sunburst(
                df_sun,
                path=sun_cols,
                color=color_col,
                color_continuous_scale='Ice'
            )
        fig_sun.update_layout(title="", showlegend=True, height=800, width=800,
                              plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
        fig_sun.update_coloraxes(showscale=False)
        sun_plot_json = fig_sun.to_json()

        # ---------- Per-trait mini gauge (SVG) ----------
        try:
            import plotly.io as pio
            HAVE_KALEIDO = True
        except Exception:
            HAVE_KALEIDO = False

        embedded_svgs = []
        details = []
        icons = []

        for _, row in data.iterrows():
            # single horizontal heat "gauge" with pointer at RAF%
            val = float(row['RAF (%)'])
            fig = go.Figure()
            z_ = np.linspace(0, 100, 100)
            fig.add_trace(go.Heatmap(
                z=[z_],
                colorscale=[[0, '#008AA5'], [1, '#F1423E']],
                showscale=False
            ))
            fig.add_trace(go.Scatter(
                x=[val], y=[0.9], mode='markers',
                marker=dict(symbol='triangle-down', size=30, color='#434343')
            ))
            fig.update_layout(
                width=1000, height=220, showlegend=False,
                xaxis=dict(range=[0, 100], showgrid=False, zeroline=False, title=f"Your Genetic Risk  {val:.2f} (%)"),
                yaxis=dict(visible=False),
                plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)',
                font=dict(size=16), margin=dict(t=60, b=100, l=100, r=100, pad=0)
            )

            if HAVE_KALEIDO:
                try:
                    svg_bytes = pio.to_image(fig, format="svg")
                    svg_string = svg_bytes.decode("utf-8")
                except Exception:
                    svg_string = "<div>Chart rendering requires kaleido. Please install: pip install -U kaleido</div>"
            else:
                svg_string = "<div>Chart rendering requires kaleido. Please install: pip install -U kaleido</div>"

            embedded_svgs.append(svg_string)

            # Optional icon per group
            icon_svg = ""
            if "Groups of Disease/Trait" in row and isinstance(row["Groups of Disease/Trait"], str):
                icon_path = os.path.join("data", "Group of disease traits", f"{row['Groups of Disease/Trait']}.svg")
                if os.path.exists(icon_path):
                    try:
                        with open(icon_path, "r", encoding="utf-8") as f:
                            icon_svg = f.read()
                    except Exception:
                        icon_svg = ""
            icons.append(icon_svg)

            # Detail tuple
            details.append((
                row.get('DISEASE/TRAIT', ''),
                row.get('REGION', ''),
                row.get('SNPS', ''),
                row.get('MAPPED_GENE', ''),
                row.get('Groups of Disease/Trait', ''),
                row.get('MAPPED_TRAIT_DESCRIPTION', '')
            ))

        # ---------- Variants donut cards (SNP/INS/DEL/COMPLEX) ----------
        type_counts = self.vcf_report['TYPE'].value_counts()
        total_variant = int(self.vcf_report.shape[0])
        type_pct = (type_counts / max(total_variant, 1) * 100).round(2)

        # Ensure all four keys exist
        for k in ["SNPs", "INS", "DEL", "COMPLEX"]:
            if k not in type_pct.index:
                type_pct.loc[k] = 0.0
                type_counts.loc[k] = 0

        # Order
        type_pct = type_pct[["SNPs", "INS", "DEL", "COMPLEX"]]
        type_counts = type_counts[["SNPs", "INS", "DEL", "COMPLEX"]]

        donut_svgs = []
        try:
            import plotly.io as pio  # already imported above, but keep local for clarity
        except Exception:
            pio = None

        color_map = {'COMPLEX': '#FF9999', 'DEL': '#FF7F3E', 'INS': '#3D527D', 'SNPs': '#FFB854'}
        for typ in type_pct.index:
            value = float(type_pct.loc[typ])
            fig2 = go.Figure()
            fig2.add_trace(go.Pie(values=[value, max(0.0, 100 - value)], hole=0.6,
                                  marker=dict(colors=[color_map[typ], "rgba(0,0,0,0)"]),
                                  direction="clockwise", textinfo="none", showlegend=False))
            fig2.add_trace(go.Pie(values=[max(0.0, 100 - value), value], hole=0.7,
                                  marker=dict(colors=["lightgray", "rgba(0,0,0,0)"]),
                                  textinfo="none", showlegend=False))
            fig2.add_annotation(text=f"<b>{value:.2f}%</b>", showarrow=False, xref="paper", yref="paper", x=0.5, y=0.5, font=dict(size=24))
            fig2.add_annotation(text=f'<b>{int(type_counts.loc[typ]):,} Positions</b>', showarrow=False, xref="paper", yref="paper", x=0.5, y=-0.38, font=dict(size=22))
            fig2.update_layout(title=f"<b>{typ}</b>", height=400, width=400, showlegend=False,
                               title_x=0.5, title_y=0.1, font=dict(size=18),
                               plot_bgcolor='rgba(0,0,0,0)', paper_bgcolor='rgba(0,0,0,0)')
            if pio is not None:
                try:
                    svg_bytes = pio.to_image(fig2, format="svg")
                    donut_svgs.append(svg_bytes.decode("utf-8"))
                except Exception:
                    donut_svgs.append("<div>Install kaleido to render donut charts (pip install -U kaleido)</div>")
            else:
                donut_svgs.append("<div>Install kaleido to render donut charts (pip install -U kaleido)</div>")

        # Summary text
        if 'Groups of Disease/Trait' in df_sun.columns:
            df_sun_summary = (df_sun['Groups of Disease/Trait']
                              .value_counts(normalize=True)
                              .mul(100).round(2)
                              .to_dict())
            total_disease_trait = int(df_sun['Groups of Disease/Trait'].value_counts().sum())
        else:
            df_sun_summary, total_disease_trait = {}, 0

        # Logo (optional)
        logo_svg = ""
        logo_path = os.path.join('data', 'logo', 'KKUNPhI-01.svg')
        if os.path.exists(logo_path):
            try:
                with open(logo_path, 'r', encoding='utf-8') as f:
                    logo_svg = f.read()
            except Exception:
                logo_svg = ""

        # HTML Template (your original, unchanged)
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
                        font-size: clamp(2rem, 3vw, 6rem);
                        margin: 0;
                        }
                    h2 {
                        margin-top: 40px;
                        color: #2D3B71;
                        border-bottom: 3px solid #2D3B71;
                        padding-bottom: 5px;
                        display: inline-block;
                        font-size: clamp(2rem, 2.5vw, 3rem)
                        }
                    h3 {
                        font-size: clamp(1.5rem, 2vw, 2rem)
                        color: #2D3B71;
                        }

                    p {
                        margin: 10px 0;
                        font-size: clamp(0.5rem, 1.75vw, 2rem)
                        color: rgb(87, 107, 145);
                        }
                    p2 {
                        margin: 10px 0;
                        font-size: clamp(0.5rem, 0.75vw, 1.5rem)
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
                    .chart-container h2 {
                        margin-top: 0;
                        text-align: left;
                        flex-grow: 1;
                        font-size: clamp(2rem, 2.5vw, 3rem)
                        front-weight: bold;
                        text-transform: capitalize;
                        flex-align: left;
                        }
                    .chart-container-inside {
                        text-align: center;
                        background: white;
                        padding: 20px;
                        margin: 20px 0;
                        width: 100%;  /* Adjust based on the window size */
                        max-width: 98%; /* Set a max-width to prevent it from getting too big */
                        }

                    .chart-container p {
                        text-align: left;
                        color: rgb(87, 107, 145);
                        }
                    .position_inside {
                        display: flex;
                        flex-direction: column;
                        gap: 5px;
                        text-align: left;
                        }
                    .chart svg {
                        width: clamp(20rem, 80vw, 95rem);  /* Make SVG scale with container */
                        height: auto; /* Maintain aspect ratio */
                        text-align: center;
                        }
                    .download-btn {
                        display: inline-block;
                        background: rgb(87, 107, 145);
                        color: white;
                        padding: 8px 12px;
                        font-size: clamp(0.5rem, 0.75vw, 1.5rem)
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
                        font-size: clamp(8rem, 14vw, 15rem);
                        font-weight: bold;
                        color: #2D3B71;
                        text-align: center;
                        font-family: 'Poppins', sans-serif;
                        margin-top: 50px;
                        margin-bottom: 20px;
                        }
                    .subtitle {
                        font-size: calc(1.5em + 1vw);
                        color: rgb(103, 103, 103);
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
                        width: clamp(20rem, 28vw, 30rem);
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
                    h5 {
                        font-size: clamp(1rem, 2vw, 3rem);
                        font-weight: bold;
                        color: #2D3B71;
                        text-align: center;
                        font-family: 'Poppins', sans-serif;
                        margin-top: 50px;
                        margin-bottom: 20px;
                        }
                    .chart_overview {
                        text-align: center;
                        margin: 50px;
                        display: flex;
                        justify-content: space-between;
                        align-items: center;
                        flex-direction: row;}
                        gap: 5px;
                        margin-bottom: 5px;
                        align-self: center;
                        flex: 1 1 auto;
                        flex-wrap: wrap;
                        width: clamp(20rem, 28vw, 30rem);
                        }
                    #sun_plot {
                        width: clamp(16rem, 19vw, 28rem);
                        height: clamp(16rem, 19vw, 28rem);
                        }

                    .summary_text {
                        display: flex;
                        flex-direction: column;
                        gap: 5px;
                        text-align: center;
                        }
                    .summary_text h3 {
                        font-size: clamp(1.5rem, 2vw, 2rem)
                        color: #2D3B71;
                        }
                    .summary_text p {
                        font-size: clamp(1rem, 1.5vw, 2rem)
                        color: rgb(103, 103, 103);
                        }
                    .plot-container {
                        width: clamp(1rem, 1.5vw, 2rem); /* Make plot container use full page width */
                        max-width: 100%; /* Prevent plot from overflowing */
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        margin: 0 auto;
                        }
                    .plot-container-inside {
                        width: 100%;
                        max-width: 100%;
                        display: flex;
                        justify-content: center;
                        align-items: center;
                        margin: 0 auto;

                        }
                    .plot-container-inside h3 {
                        font-size: clamp(1.5rem, 2vw, 2rem)
                        color: #2D3B71;
                        text-align: center;
                        }
                    .plot-container-inside p {
                        font-size: clamp(1rem, 1.5vw, 2rem)
                        color: rgb(103, 103, 103);
                        text-align: center;
                        }
                    .plot-container-inside .position_inside {
                        display: flex;

                        }
                    .disease-number {
                        font-size: clamp(2rem, 2.5vw, 3rem)
                        font-weight: bold;
                        color: #2D3B71;
                        }
                    .disease-text {
                        font-size: clamp(1rem, 1.5vw, 2rem)
                        color: rgb(103, 103, 103);
                        }
                    .print-btn {
                        display: block;
                        margin: 20px auto;
                        padding: 10px 20px;
                        font-size: 1rem;
                        background-color: #3c79aa;
                        color: white;
                        border: none;
                        cursor: pointer;
                        border-radius: 5px;
                        box-shadow: 0 4px 10px rgba(0, 0, 0, 0.1);
                        transition: 0.3s;
                        }
                    .print-btn:hover {
                        background-color: #2D3B71;
                        }
                    @media print {
                        @page {
                            size: A4 portrait; /* Set to A4 size */
                            margin: 10mm; /* Adjust margin for better fit */
                        }
                        html, body {
                            width: 210mm;
                            height: 297mm;
                            margin: 0;
                            padding: 0;
                        }
                        body {
                            width: 100%;
                        }
                        .header-container,
                        .logo-container,
                        .variants,
                        .plot-container {
                            width: 100%;
                            max-width: 100%;
                            box-sizing: border-box;
                        }
                        h1 {
                            font-size: 48px; /* Title size optimized for A4 */
                            text-align: center;
                        }
                        h5, .print-btn, .download-btn, #sun_plot {
                            display: none; /* Hide unnecessary elements */
                        }
                        section {
                            width: 100%;
                            max-width: 100%;
                            box-sizing: border-box;
                            page-break-before: always;
                            page-break-inside: avoid;
                        }
                        img, table, svg, canvas {
                            max-width: 100%;
                            height: auto;
                        }
                        .variants {
                            display: flex;
                            justify-content: center;
                            align-items: center;
                            flex-direction: row;
                            flex-wrap: wrap;
                            gap: 10px;
                        }
                        .variant {
                            display: flex;
                            align-items: center;
                            gap: 5px;
                            margin-bottom: 5px;
                            align-self: center;
                            width: 40px;
                        }
                        sub-title {
                            font-size: 30px;
                            text-align: center;
                        }
                    }
                        /* Mobile-first responsive design */
                    @media (max-width: 768px) {
                        body {
                            margin: 20px;
                            font-size: 16px;
                        }

                        .header-container {
                            flex-direction: column;
                            text-align: center;
                        }

                        h1 {
                            font-size: clamp(1.5rem, 5vw, 3rem);
                            padding: 15px;
                        }

                        h2, h3, h5 {
                            font-size: clamp(1.2rem, 4vw, 2rem);
                            text-align: center;
                        }
                        .title {
                            font-size: clamp(3rem, 5vw, 8rem);
                        }
                        p {
                            font-size: clamp(0.8rem, 3vw, 1.5rem);
                            padding: 10px;
                        }

                        .chart-container,
                        .chart-container-inside,
                        .plot-container,
                        .plot-container-inside {
                            width: 100%;
                            max-width: 100%;
                            padding: 10px;
                        }

                        .icon-text {
                            flex-direction: column;
                            align-items: center;
                        }

                        .icon-text svg {
                            width: 80%;
                            margin-bottom: 10px;
                        }

                        .variants {
                            flex-direction: column;
                            align-items: center;
                        }

                        .variant {
                            width: 70%;
                        }

                        .download-btn, .print-btn {
                            font-size: 9px;
                            padding: 3px 6px;
                            width: 20%;
                            text-align: center;
                        }

                        .logo-container {
                            padding: 40px;
                        }

                        .summary_text, .position_inside {
                            text-align: center;
                        }

                        .chart svg {
                            width: 100%;
                            height: auto;
                        }

                        #sun_plot {
                            width: 95%;
                            height: auto;
                        }

                        hr {
                            margin: 20px 0;
                        }
                        .position_inside {
                            display: flex;
                            flex-direction: column;
                            gap: 3px;
                            padding: 3px;
                        }
                        .position_infomation {
                            display: flex;
                            flex-direction: column;
                            gap: 3px;
                            padding: 3px;
                            flex-align: left;
                            }
                    }

                </style>
                <script>

                    function resizePlot() {
                        let plotContainer = document.getElementById('plot-container');
                        let plotDiv = document.getElementsByClassName('plotly-graph-div')[0];

                        if (plotDiv) {
                            Plotly.relayout(plotDiv, {
                                width: plotContainer.clientWidth,
                                height: window.innerHeight * 0.8
                            });
                        }
                    }

                window.addEventListener('resize', resizePlot);
                window.onload = resizePlot; // Run once when the page loads
                    function resizeChart() {
                        var chart = document.getElementById('chart_1');
                        Plotly.relayout(chart, {
                            width: window.innerWidth * 0.9,  // 90% of the window width
                            height: window.innerHeight * 0.3 // Adjust height dynamically
                        });
                    }

                    window.addEventListener('resize', resizeChart);
                    window.onload = resizeChart; // Run once when the page loads

                    function printReport() {
                        window.print();
                    }
                </script>

            </head>

            <body>
            <section>
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
            </section>
            <section>
                <div>
                    <h2>Overview</h2>
                        <p>This report presents the findings of whole genome sequencing
                            analysis. It provides detailed associations between genetic
                            variants and various diseases and traits, highlighting the mapped
                            genes, relevant SNPs, and associated gene regions.
                        </p>
                <div class="overview">
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
            </section>
            <section>
                <div>
                <h5>Overview of Disease/Trait in Your Genome</h5>
                </div>
                <div class="chart_overview" id="sun_plot"></div>
            <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
            <script id="sunburst-data" type="application/json">
                {{ sun_plot_ | safe }}
            </script>
            <script>
                // Read JSON data from the script tag
                var jsonData = JSON.parse(document.getElementById("sunburst-data").textContent);

                // Plot the Sunburst chart
                Plotly.newPlot('sun_plot', jsonData.data, jsonData.layout);
                Plotly.relayout('sun_plot', autosize = true);
                //update figure size base on windows size
                window.addEventListener('resize', function() {
                    Plotly.relayout('sun_plot', {
                        //adjust width and height based on grid-over-view-container
                        width: document.querySelector('.chart_overview').offsetWidth,
                        height: document.querySelector('.chart_overview').offsetHeight
                        });
                    });
                //resize colorbar
                window.onload = resizePlot; // Run once when the page loads

            </script>
            </section>
            <section>
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
                            <p><b>Mapped Gene:</b> Genes mapped near or overlapping the SNPs.</p>
                            <p><b>Chromosomal region:</b> The genomic region associated with the trait or disease.</p>
                            <p><b>Risk Allele Frequency (%):</b> The frequency of the risk allele in the population.</p>
                </div>
            </section>
                <hr>
                {% for (title_, region_, snps_, mapped_gene_, group_trait_, description_trait_), svg_, icon_ in data_source %}
                <section>
                    <div class="chart-container">
                        <h2>{{ title_ }}</h2>
                        <p>{{ description_trait_ }}</p>
                        <div class="chart-container-inside">
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

                    </div>
                    <hr>
                </section>
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


                <div>
                    <button class="print-btn" onclick="printReport()">Print Report</button>
                </div>
            </body>
            </html>
        """

        # Render
        template = Template(html_template)
        rendered_html = template.render(
            data_source=zip(details, embedded_svgs, icons),
            count_variant=f'{total_variant:,.0f}',
            variant_1=donut_svgs[0], variant_2=donut_svgs[1],
            variant_3=donut_svgs[2], variant_4=donut_svgs[3],
            logo_=logo_svg,
            sun_plot_=sun_plot_json,
            disease_trait_summary=df_sun_summary,
            total_disease_trait_=total_disease_trait
        )

        os.makedirs(self.report_path, exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(rendered_html)

        print(f"Report saved to {output_path}")

    def generate_report(self):
        self.prepare_report_data()
        self.generate_html_report()

if __name__ == '__main__':
    barcode = 'bc01'
    vcf_file = f"medaka/sort-medaka-{barcode}/medaka.sorted.vcf"
    gwas_file = "data/gwas_database_with_description_expanded.csv.gz"
    output_file = f"medaka/{barcode}"

    mapper = MapGWASSNPs(vcf_file, gwas_file, output_file, cut_off_qual=60, filt_nr_disease=True)
    mapper.map_snps()
    mapper.generate_report()