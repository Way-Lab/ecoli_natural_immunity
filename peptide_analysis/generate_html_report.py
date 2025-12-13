#!/usr/bin/env python3
"""
Generate a comprehensive HTML report with embedded graphics and navigation.
"""

import base64
import os
from pathlib import Path


def encode_image_to_base64(image_path):
    """Encode an image file to base64 for embedding in HTML."""
    if not os.path.exists(image_path):
        return None
    with open(image_path, 'rb') as f:
        return base64.b64encode(f.read()).decode('utf-8')


def read_file(file_path):
    """Read a text file."""
    if not os.path.exists(file_path):
        return None
    with open(file_path, 'r') as f:
        return f.read()


def markdown_to_html_simple(md_text):
    """Simple markdown to HTML converter for basic formatting."""
    if not md_text:
        return ""

    html = md_text

    # Headers
    html = html.replace('# ', '<h1>').replace('\n', '</h1>\n', 1)
    for i in range(6, 0, -1):
        html = html.replace('#' * i + ' ', f'<h{i}>')
        # Close headers at newlines

    # Better header handling
    lines = html.split('\n')
    processed = []
    for line in lines:
        if line.startswith('<h') and not line.endswith('>'):
            # Find where heading content ends
            processed.append(line + line[2] + '>')
        else:
            processed.append(line)
    html = '\n'.join(processed)

    # Code blocks
    html = html.replace('`', '<code>').replace('</code></code>', '</code>')

    # Bold
    import re
    html = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', html)

    # Lists
    html = re.sub(r'^\- (.+)$', r'<li>\1</li>', html, flags=re.MULTILINE)
    html = re.sub(r'(<li>.*</li>\n?)+', r'<ul>\g<0></ul>', html)

    # Paragraphs
    html = re.sub(r'\n\n', r'</p><p>', html)
    html = '<p>' + html + '</p>'

    # Clean up
    html = html.replace('<p><h', '<h').replace('</h1></p>', '</h1>')
    html = html.replace('</h2></p>', '</h2>').replace('</h3></p>', '</h3>')
    html = html.replace('</h4></p>', '</h4>').replace('<p></p>', '')
    html = html.replace('<p><ul>', '<ul>').replace('</ul></p>', '</ul>')
    html = html.replace('<p>---</p>', '<hr>')
    html = html.replace('<p>===', '<hr><p>===')

    return html


def create_html_report():
    """Create the main HTML report."""

    # Load images
    conservation_img = encode_image_to_base64('enterobacteriaceae_conservation.png')
    variability_img = encode_image_to_base64('enterobacteriaceae_variability.png')

    # Load text reports
    summary_md = read_file('enterobacteriaceae_peptide_summary.md')
    species_md = read_file('species_variant_analysis.md')
    readme_md = read_file('PEPTIDE_ANALYSIS_README.md')

    # Load CSV data
    import csv
    summary_csv = []
    if os.path.exists('enterobacteriaceae_summary.csv'):
        with open('enterobacteriaceae_summary.csv', 'r') as f:
            reader = csv.DictReader(f)
            summary_csv = list(reader)

    # Start HTML
    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>OmpA Peptide Variability Analysis Report</title>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background: #f5f5f5;
        }

        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
        }

        header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }

        header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }

        header p {
            font-size: 1.2em;
            opacity: 0.9;
        }

        nav {
            background: #2d3748;
            padding: 15px 40px;
            position: sticky;
            top: 0;
            z-index: 100;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        nav a {
            color: white;
            text-decoration: none;
            margin-right: 30px;
            padding: 8px 16px;
            border-radius: 4px;
            transition: background 0.3s;
            display: inline-block;
        }

        nav a:hover {
            background: #4a5568;
        }

        .content {
            padding: 40px;
        }

        .section {
            margin-bottom: 50px;
        }

        h1 {
            color: #2d3748;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 2em;
        }

        h2 {
            color: #4a5568;
            margin-top: 30px;
            margin-bottom: 15px;
            font-size: 1.5em;
        }

        h3 {
            color: #667eea;
            margin-top: 20px;
            margin-bottom: 10px;
            font-size: 1.2em;
        }

        .summary-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }

        .card {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }

        .card h3 {
            color: white;
            font-size: 1em;
            margin: 0 0 10px 0;
            opacity: 0.9;
        }

        .card .value {
            font-size: 2.5em;
            font-weight: bold;
        }

        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        th, td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e2e8f0;
        }

        th {
            background: #2d3748;
            color: white;
            font-weight: 600;
        }

        tr:hover {
            background: #f7fafc;
        }

        .highlight-box {
            background: #edf2f7;
            border-left: 4px solid #667eea;
            padding: 20px;
            margin: 20px 0;
            border-radius: 4px;
        }

        .finding {
            background: white;
            border: 1px solid #e2e8f0;
            padding: 20px;
            margin: 15px 0;
            border-radius: 8px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }

        .finding h4 {
            color: #667eea;
            margin-bottom: 10px;
        }

        .badge {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 600;
            margin-left: 10px;
        }

        .badge-high {
            background: #fc8181;
            color: white;
        }

        .badge-medium {
            background: #f6ad55;
            color: white;
        }

        .badge-low {
            background: #68d391;
            color: white;
        }

        .badge-not-found {
            background: #cbd5e0;
            color: #2d3748;
        }

        .peptide-sequence {
            font-family: 'Courier New', monospace;
            background: #2d3748;
            color: #68d391;
            padding: 10px 15px;
            border-radius: 4px;
            font-size: 1.1em;
            display: inline-block;
            margin: 10px 0;
        }

        .image-container {
            margin: 30px 0;
            text-align: center;
        }

        .image-container img {
            max-width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }

        .image-caption {
            margin-top: 10px;
            font-style: italic;
            color: #718096;
        }

        code {
            background: #edf2f7;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
        }

        ul {
            margin: 15px 0;
            padding-left: 30px;
        }

        li {
            margin: 8px 0;
        }

        .species-grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
            gap: 10px;
            margin: 20px 0;
        }

        .species-item {
            background: #f7fafc;
            padding: 10px;
            border-radius: 4px;
            border-left: 3px solid #667eea;
        }

        footer {
            background: #2d3748;
            color: white;
            padding: 30px 40px;
            text-align: center;
        }

        .download-section {
            background: #edf2f7;
            padding: 20px;
            margin: 30px 0;
            border-radius: 8px;
        }

        .download-section a {
            color: #667eea;
            text-decoration: none;
            font-weight: 600;
        }

        .download-section a:hover {
            text-decoration: underline;
        }

        hr {
            border: none;
            border-top: 2px solid #e2e8f0;
            margin: 30px 0;
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>OmpA Peptide Variability Analysis</h1>
            <p>Comprehensive analysis across E. coli and Enterobacteriaceae</p>
        </header>

        <nav>
            <a href="#summary">Summary</a>
            <a href="#key-findings">Key Findings</a>
            <a href="#visualizations">Visualizations</a>
            <a href="#detailed-results">Detailed Results</a>
            <a href="#species-analysis">Species Analysis</a>
            <a href="#methods">Methods</a>
        </nav>

        <div class="content">
            <section id="summary" class="section">
                <h1>Executive Summary</h1>

                <div class="summary-cards">
                    <div class="card">
                        <h3>Total Sequences</h3>
                        <div class="value">130</div>
                    </div>
                    <div class="card">
                        <h3>Species Analyzed</h3>
                        <div class="value">12</div>
                    </div>
                    <div class="card">
                        <h3>Peptides Tested</h3>
                        <div class="value">13</div>
                    </div>
                    <div class="card">
                        <h3>Peptides Found</h3>
                        <div class="value">6</div>
                    </div>
                </div>

                <div class="highlight-box">
                    <h3>Analysis Overview</h3>
                    <p>This comprehensive analysis examined 13 specific peptide sequences across 130 OmpA protein sequences from <strong>E. coli</strong> and other <strong>Enterobacteriaceae</strong> species including Klebsiella, Enterobacter, Citrobacter, Salmonella, Shigella, Proteus, and Serratia.</p>
                    <p style="margin-top: 15px;"><strong>Key Question:</strong> How variable are these peptide sequences across different species and strains, and which ones are conserved enough for broad-spectrum applications or specific enough for diagnostics?</p>
                </div>
"""

    # Add summary table
    if summary_csv:
        html += """
                <h2>Peptide Conservation Summary</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Peptide</th>
                            <th>Length (aa)</th>
                            <th>Avg Conservation</th>
                            <th>Min Conservation</th>
                            <th>Variable Positions</th>
                            <th>Assessment</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        for row in summary_csv:
            avg_cons = float(row['Avg_Conservation_%'])
            var_pos = int(row['Variable_Positions'])

            if var_pos == 0:
                badge = '<span class="badge badge-low">Highly Conserved</span>'
            elif var_pos <= 2:
                badge = '<span class="badge badge-medium">Moderately Variable</span>'
            else:
                badge = '<span class="badge badge-high">Highly Variable</span>'

            html += f"""
                        <tr>
                            <td><strong>{row['Peptide']}</strong></td>
                            <td>{row['Length']}</td>
                            <td>{avg_cons:.1f}%</td>
                            <td>{row['Min_Conservation_%']}%</td>
                            <td>{var_pos}</td>
                            <td>{badge}</td>
                        </tr>
"""
        html += """
                    </tbody>
                </table>
"""

    html += """
            </section>

            <section id="key-findings" class="section">
                <h1>Key Findings</h1>

                <div class="finding">
                    <h4>🎯 Finding #1: E. coli-Specific Marker Identified</h4>
                    <p><strong>Peptide_1</strong> is an excellent species-specific marker:</p>
                    <div class="peptide-sequence">WSQYHDTGFINNNGPTHEN</div>
                    <ul>
                        <li><strong>Present in 100% of E. coli sequences</strong></li>
                        <li>Also found in closely related Salmonella (100%) and Shigella (100%)</li>
                        <li><strong>Completely absent</strong> from Klebsiella, Citrobacter, Proteus, and Serratia</li>
                        <li>Only 2 variable positions (moderately conserved)</li>
                        <li><strong>Application:</strong> Ideal for E. coli-specific diagnostics or targeted therapeutics</li>
                    </ul>
                </div>

                <div class="finding">
                    <h4>🔬 Finding #2: Two Major E. coli OmpA Variants</h4>
                    <p>Analysis reveals approximately a 50/50 split in E. coli strains:</p>
                    <ul>
                        <li><strong>Variant A (~50%):</strong> Contains GRMPYKGSVENGAYKA (Peptide_2)</li>
                        <li><strong>Variant B (~50%):</strong> Contains GRMPYKGDNINGAYKA (Peptide_2b)</li>
                        <li>Key differences at positions 8 (S↔D), 9 (V↔N), and 10 (E↔I)</li>
                        <li><strong>Implication:</strong> Suggests balanced polymorphism or recent evolutionary divergence</li>
                    </ul>
                </div>

                <div class="finding">
                    <h4>📊 Finding #3: Most Prevalent and Conserved Peptide</h4>
                    <p><strong>Peptide_2</strong> shows the widest distribution:</p>
                    <div class="peptide-sequence">GRMPYKGSVENGAYKA</div>
                    <ul>
                        <li>Found in <strong>60.8% of all sequences</strong> (79/130)</li>
                        <li>Present in E. coli (100%), Klebsiella (100%), Salmonella (100%), Shigella (100%)</li>
                        <li>5 variable positions - allows for strain/species discrimination</li>
                        <li><strong>Application:</strong> Good candidate for broad-spectrum approaches with variant-specific modifications</li>
                    </ul>
                </div>

                <div class="finding">
                    <h4>🧬 Finding #4: Salmonella-Specific Substitutions</h4>
                    <p>Salmonella shows unique amino acid substitutions in the Peptide_1 region:</p>
                    <ul>
                        <li>60%: WSQYHDTGFI<strong>N</strong>N<strong>D</strong>GPTHEN (N→D at position 13)</li>
                        <li>40%: WSQYHDTGFI<strong>H</strong>N<strong>D</strong>GPTHEN (N→H at position 11, N→D at position 13)</li>
                        <li>These variants are <strong>unique to Salmonella</strong> - not found in any other species</li>
                        <li><strong>Application:</strong> Can be used to distinguish Salmonella from E. coli</li>
                    </ul>
                </div>

                <div class="finding">
                    <h4>❌ Finding #5: Seven Peptides Not Found</h4>
                    <p>The following peptides were <strong>not detected</strong> in any OmpA sequence:</p>
                    <ul>
                        <li>Peptide_1S, 2S, 3aS, 3bS, 4S, 3BS2, 2bS</li>
                        <li>These likely represent synthetic/scrambled variants or sequences from other proteins</li>
                        <li>Could serve as <strong>negative controls</strong> in experimental work</li>
                    </ul>
                </div>
            </section>
"""

    # Add visualizations
    if conservation_img or variability_img:
        html += """
            <section id="visualizations" class="section">
                <h1>Visualizations</h1>
"""

        if conservation_img:
            html += f"""
                <div class="image-container">
                    <h2>Conservation Across Positions</h2>
                    <img src="data:image/png;base64,{conservation_img}" alt="Conservation Plot">
                    <p class="image-caption">
                        Bar plots showing conservation percentage at each position for all peptides.
                        Green bars (≥90%) indicate highly conserved positions, orange (75-90%) moderately conserved,
                        and red (<75%) variable positions.
                    </p>
                </div>
"""

        if variability_img:
            html += f"""
                <div class="image-container">
                    <h2>Amino Acid Variability Heatmaps</h2>
                    <img src="data:image/png;base64,{variability_img}" alt="Variability Heatmap">
                    <p class="image-caption">
                        Heatmaps showing the frequency of each amino acid at each position.
                        Darker colors indicate higher frequency. This visualization reveals which
                        specific amino acid substitutions occur at variable positions.
                    </p>
                </div>
"""

        html += """
            </section>
"""

    # Add detailed results section
    html += """
            <section id="detailed-results" class="section">
                <h1>Detailed Peptide Analysis</h1>
"""

    # Parse and add peptide details
    peptide_details = [
        {
            'name': 'Peptide_1',
            'seq': 'WSQYHDTGFINNNGPTHEN',
            'found': '38.5%',
            'variants': 4,
            'var_pos': 2,
            'status': 'Moderately Conserved',
            'badge': 'medium',
            'notes': 'E. coli-specific marker. Found in 100% of E. coli, Salmonella, and Shigella. Absent from Klebsiella, Enterobacter, Citrobacter, Proteus, and Serratia.'
        },
        {
            'name': 'Peptide_2',
            'seq': 'GRMPYKGSVENGAYKA',
            'found': '60.8%',
            'variants': 6,
            'var_pos': 5,
            'status': 'Highly Variable',
            'badge': 'high',
            'notes': 'Most prevalent peptide. Present across multiple species. 5 variable positions allow for strain discrimination.'
        },
        {
            'name': 'Peptide_3a',
            'seq': 'KSNVYGKN',
            'found': '25.4%',
            'variants': 3,
            'var_pos': 3,
            'status': 'Moderately Conserved',
            'badge': 'medium',
            'notes': 'Mutually exclusive with Peptide_3b. Found in different E. coli strains.'
        },
        {
            'name': 'Peptide_3b',
            'seq': 'KANVPGGASFKD',
            'found': '22.3%',
            'variants': 4,
            'var_pos': 4,
            'status': 'Highly Variable',
            'badge': 'high',
            'notes': 'Alternative variant to Peptide_3a. 4 variable positions.'
        },
        {
            'name': 'Peptide_4',
            'seq': 'TNNIGDAHTIGTRPDNGM',
            'found': '53.8%',
            'variants': 3,
            'var_pos': 4,
            'status': 'Highly Variable',
            'badge': 'high',
            'notes': 'Well distributed across species. Good candidate for broad-spectrum applications despite variability.'
        },
        {
            'name': 'Peptide_2b',
            'seq': 'GRMPYKGDNINGAYKA',
            'found': '46.2%',
            'variants': 5,
            'var_pos': 4,
            'status': 'Highly Variable',
            'badge': 'high',
            'notes': 'Alternative variant of Peptide_2. Found in ~50% of E. coli strains (the other 50% have Peptide_2).'
        },
    ]

    for peptide in peptide_details:
        badge_class = f"badge-{peptide['badge']}"
        html += f"""
                <div class="finding">
                    <h4>{peptide['name']}<span class="badge {badge_class}">{peptide['status']}</span></h4>
                    <div class="peptide-sequence">{peptide['seq']}</div>
                    <table style="margin-top: 15px;">
                        <tr>
                            <td style="border: none;"><strong>Found in:</strong></td>
                            <td style="border: none;">{peptide['found']} of sequences</td>
                        </tr>
                        <tr>
                            <td style="border: none;"><strong>Unique variants:</strong></td>
                            <td style="border: none;">{peptide['variants']}</td>
                        </tr>
                        <tr>
                            <td style="border: none;"><strong>Variable positions:</strong></td>
                            <td style="border: none;">{peptide['var_pos']}</td>
                        </tr>
                    </table>
                    <p style="margin-top: 15px;"><strong>Notes:</strong> {peptide['notes']}</p>
                </div>
"""

    # Add NOT FOUND peptides
    html += """
                <h2 style="margin-top: 40px;">Peptides Not Found in OmpA</h2>
                <div class="highlight-box">
                    <p>The following 7 peptides were <strong>not detected</strong> in any of the 130 OmpA sequences:</p>
                    <ul>
                        <li><code>Peptide_1S</code>: TQPNSHNDINTNHGYGWEF</li>
                        <li><code>Peptide_2S</code>: RSMAGKGPKNAGYVYE</li>
                        <li><code>Peptide_3aS</code>: KYNGVNSK</li>
                        <li><code>Peptide_3bS</code>: FKAKDNGVGSPA</li>
                        <li><code>Peptide_4S</code>: MTNTHNTAIGDRNIGPGD</li>
                        <li><code>Peptide_3BS2</code>: GVSFDKAKNAPG</li>
                        <li><code>Peptide_2bS</code>: ANGGGYIADYNKRKPM</li>
                    </ul>
                    <p style="margin-top: 15px;"><strong>Interpretation:</strong> These sequences likely represent synthetic/scrambled variants designed as negative controls, or they may be from non-OmpA proteins.</p>
                </div>
            </section>
"""

    # Species analysis section
    html += """
            <section id="species-analysis" class="section">
                <h1>Species Distribution Analysis</h1>

                <h2>Species Included</h2>
                <div class="species-grid">
                    <div class="species-item"><strong>Escherichia coli</strong><br>14 sequences</div>
                    <div class="species-item"><strong>Klebsiella pneumoniae</strong><br>10 sequences</div>
                    <div class="species-item"><strong>Klebsiella oxytoca</strong><br>10 sequences</div>
                    <div class="species-item"><strong>Enterobacter cloacae</strong><br>8 sequences</div>
                    <div class="species-item"><strong>Citrobacter freundii</strong><br>9 sequences</div>
                    <div class="species-item"><strong>Citrobacter koseri</strong><br>10 sequences</div>
                    <div class="species-item"><strong>Salmonella enterica</strong><br>10 sequences</div>
                    <div class="species-item"><strong>Shigella flexneri</strong><br>8 sequences</div>
                    <div class="species-item"><strong>Shigella sonnei</strong><br>10 sequences</div>
                    <div class="species-item"><strong>Proteus mirabilis</strong><br>10 sequences</div>
                    <div class="species-item"><strong>Serratia marcescens</strong><br>9 sequences</div>
                </div>

                <h2>Key Species-Specific Patterns</h2>

                <div class="finding">
                    <h4>E. coli / Shigella / Salmonella Cluster</h4>
                    <p>These closely related species share Peptide_1 (WSQYHDTGFINNNGPTHEN):</p>
                    <ul>
                        <li>E. coli: 100% (14/14)</li>
                        <li>Shigella spp: 100% (18/18)</li>
                        <li>Salmonella: 100% (10/10) - but with unique substitutions</li>
                    </ul>
                    <p><strong>Salmonella-specific markers:</strong> N→D substitution at position 13 (unique to Salmonella)</p>
                </div>

                <div class="finding">
                    <h4>Klebsiella Species</h4>
                    <p>Klebsiella shows different conservation patterns:</p>
                    <ul>
                        <li><strong>DO NOT</strong> contain Peptide_1 (0% presence)</li>
                        <li><strong>DO</strong> contain Peptide_2 (100% in K. pneumoniae and K. oxytoca)</li>
                        <li>Can be distinguished from E. coli by absence of Peptide_1</li>
                    </ul>
                </div>

                <div class="finding">
                    <h4>Other Enterobacteriaceae</h4>
                    <ul>
                        <li><strong>Citrobacter:</strong> Low presence of most peptides; distinct from E. coli</li>
                        <li><strong>Proteus & Serratia:</strong> Most distantly related; different OmpA variants</li>
                        <li><strong>Enterobacter:</strong> Intermediate patterns</li>
                    </ul>
                </div>
            </section>

            <section id="methods" class="section">
                <h1>Methods & Data</h1>

                <h2>Computational Approach</h2>
                <ol>
                    <li><strong>Sequence Collection:</strong> Downloaded 120 OmpA protein sequences from NCBI using BioPython Entrez, plus 10 existing sequences = 130 total</li>
                    <li><strong>Peptide Mapping:</strong> Searched each of the 13 peptides against all 130 sequences using fuzzy matching (allowing up to 2-3 mismatches)</li>
                    <li><strong>Conservation Analysis:</strong> Calculated position-specific amino acid frequencies and conservation scores</li>
                    <li><strong>Species Analysis:</strong> Grouped sequences by species/genus and analyzed variant distributions</li>
                    <li><strong>Visualization:</strong> Generated conservation plots and variability heatmaps using matplotlib/seaborn</li>
                </ol>

                <h2>Software & Tools</h2>
                <ul>
                    <li>Python 3 with BioPython for sequence analysis</li>
                    <li>NCBI Entrez API for sequence retrieval</li>
                    <li>Custom scripts for peptide variability analysis</li>
                    <li>matplotlib and seaborn for visualization</li>
                </ul>

                <h2>Analysis Scripts</h2>
                <div class="download-section">
                    <p>All analysis scripts are available in the project directory:</p>
                    <ul>
                        <li><code>analyze_peptide_variability.py</code> - Main analysis script</li>
                        <li><code>visualize_peptide_conservation.py</code> - Generates plots and heatmaps</li>
                        <li><code>analyze_species_variants.py</code> - Species-specific analysis</li>
                        <li><code>download_enterobacteriaceae_ompA.py</code> - NCBI sequence downloader</li>
                        <li><code>generate_final_summary.py</code> - Summary report generator</li>
                    </ul>
                </div>

                <h2>Data Files</h2>
                <div class="download-section">
                    <ul>
                        <li><code>combined_ompA_sequences.fasta</code> - All 130 sequences analyzed</li>
                        <li><code>enterobacteriaceae_peptide_analysis_variability_report.txt</code> - Full detailed report</li>
                        <li><code>species_variant_analysis.md</code> - Species-specific patterns</li>
                        <li><code>enterobacteriaceae_summary.csv</code> - Summary statistics</li>
                    </ul>
                </div>

                <h2>Recommendations for Future Work</h2>
                <ol>
                    <li><strong>Structural Mapping:</strong> Map variable positions onto OmpA 3D structure to assess surface accessibility and functional impact</li>
                    <li><strong>Epitope Prediction:</strong> Run B-cell and T-cell epitope prediction on conserved peptides</li>
                    <li><strong>Clinical Correlation:</strong> Correlate Peptide_2 vs Peptide_2b variants with pathogenicity, antibiotic resistance, or environmental vs clinical origin</li>
                    <li><strong>Extended Sampling:</strong> Include more diverse geographic and temporal isolates</li>
                    <li><strong>Experimental Validation:</strong> Test antibody cross-reactivity using the identified conserved peptides</li>
                </ol>
            </section>
        </div>

        <footer>
            <p><strong>OmpA Peptide Variability Analysis Report</strong></p>
            <p>Analysis Date: October 30, 2025</p>
            <p>Total Sequences: 130 | Species: 12 | Peptides Analyzed: 13</p>
        </footer>
    </div>
</body>
</html>
"""

    # Write HTML file
    with open('OmpA_Analysis_Report.html', 'w') as f:
        f.write(html)

    print("HTML report generated: OmpA_Analysis_Report.html")


if __name__ == "__main__":
    create_html_report()
