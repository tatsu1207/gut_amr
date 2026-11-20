#!/usr/bin/env python3
"""Generate manuscript PDF with embedded figures and tables using reportlab."""

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch, cm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    PageBreak, KeepTogether
)
from reportlab.lib import colors
import os
import re

MANUSCRIPT = os.path.expanduser("~/vault/kjm/miseq/manuscript.txt")
OUTFILE = os.path.expanduser("~/vault/kjm/miseq/manuscript_with_figures.pdf")

# Figure paths
FIGURES = {
    "Fig. 1": os.path.expanduser("~/vault/kjm/miseq/dada2_output/figures/fig1_physiology.png"),
    "Fig. 2": os.path.expanduser("~/vault/kjm/miseq/dada2_output/figures/fig2_diversity.png"),
    # Multi-panel figures are combined by compose_multipanel_figures.py (run it first)
    "Fig. 3": os.path.expanduser("~/vault/kjm/miseq/dada2_output/figures/fig3_combined.png"),
    "Fig. 4": os.path.expanduser("~/vault/kjm/hiseq/figures/fig4_combined.png"),
    "Fig. 5": os.path.expanduser("~/vault/kjm/hiseq/figures/fig5_combined.png"),
    "Fig. 6": os.path.expanduser("~/vault/kjm/hiseq/figures/fig6_butyrate_capacity.png"),
}

def _read_legends(path):
    """Full figure legends, taken from the 'Figure legends' section of manuscript.txt.
    npj requires each legend on the same page as its figure; reading them from the
    manuscript keeps the PDF captions identical to the submitted legends."""
    legends, on = {}, False
    for ln in open(path):
        t = ln.strip()
        if t == 'Figure legends':
            on = True; continue
        if on and t == 'Additional files':
            break
        m = re.match(r'(Fig\. \d+)\.', t) if on else None
        if m:
            legends[m.group(1)] = t
    return legends


FIGURE_CAPTIONS = _read_legends(MANUSCRIPT)
assert len(FIGURE_CAPTIONS) == 6, sorted(FIGURE_CAPTIONS)


def parse_manuscript(filepath):
    """Parse manuscript into sections."""
    with open(filepath) as f:
        text = f.read()
    return text


def build_pdf():
    doc = SimpleDocTemplate(
        OUTFILE,
        pagesize=letter,
        topMargin=1*inch, bottomMargin=1*inch,
        leftMargin=1*inch, rightMargin=1*inch,
    )

    styles = getSampleStyleSheet()

    # Custom styles
    styles.add(ParagraphStyle(
        'ManuscriptTitle', parent=styles['Title'],
        fontSize=16, leading=20, spaceAfter=12, alignment=TA_CENTER,
    ))
    styles.add(ParagraphStyle(
        'Authors', parent=styles['Normal'],
        fontSize=11, leading=14, alignment=TA_CENTER, spaceAfter=6,
    ))
    styles.add(ParagraphStyle(
        'Affiliation', parent=styles['Normal'],
        fontSize=9, leading=12, alignment=TA_CENTER, spaceAfter=4,
        textColor=colors.grey,
    ))
    styles.add(ParagraphStyle(
        'SectionHead', parent=styles['Heading1'],
        fontSize=14, leading=18, spaceBefore=18, spaceAfter=8,
        textColor=colors.HexColor('#2c3e50'),
    ))
    styles.add(ParagraphStyle(
        'SubsectionHead', parent=styles['Heading2'],
        fontSize=12, leading=15, spaceBefore=14, spaceAfter=6,
        textColor=colors.HexColor('#34495e'),
    ))
    styles.add(ParagraphStyle(
        'BodyText2', parent=styles['Normal'],
        fontSize=10, leading=14, alignment=TA_JUSTIFY, spaceAfter=8,
        firstLineIndent=0,
    ))
    styles.add(ParagraphStyle(
        'FigCaption', parent=styles['Normal'],
        fontSize=8.5, leading=11, alignment=TA_JUSTIFY, spaceAfter=12,
        textColor=colors.HexColor('#555555'),
    ))
    styles.add(ParagraphStyle(
        'RefStyle', parent=styles['Normal'],
        fontSize=8.5, leading=11, spaceAfter=2,
    ))

    story = []

    # Read manuscript
    with open(MANUSCRIPT) as f:
        lines = f.readlines()

    # Track sections to insert figures
    section_heads = [
        'Abstract', 'Background', 'Results', 'Discussion', 'Conclusions',
        'Methods', 'List of abbreviations', 'Declarations', 'References',
        'Figure legends', 'Additional files',
    ]
    subsection_heads = [
        'HFD-induced obesity and metabolic changes',
        'Antibiotics disrupt gut microbiome diversity',
        'Diet-associated community differences persist',
        'Antibiotic treatment drives Lactococcus',
        'Akkermansia expansion characterizes',
        'Resistome dynamics reveal',
        'Lactococcus abundance correlates',
        'Genome-resolved metagenomics',
        'Plasmid analysis reveals',
        'Mobile genetic element association',
        'Butyrate-synthesis capacity collapses',
        'Ethics statement', 'Animals and experimental design',
        'Serum cholesterol analysis', 'Quantitative analysis',
        'Fecal DNA extraction', '16S rRNA gene amplicon',
        '16S rRNA gene sequence', 'Shotgun metagenomic sequencing',
        'Metagenomic data processing', 'Metagenomic assembly',
        'Metagenomic binning', 'Taxonomic classification',
        'Plasmid and phage', 'Mobile genetic element detection',
        'Butyrate-synthesis gene analysis', 'Use of large language models',
        'Statistical analysis',
        'Code availability',
        'Ethics approval', 'Consent for publication',
        'Competing interests', 'Funding', 'Authors\' contributions',
        'Acknowledgements', 'Availability of data',
    ]

    # Figure insertion points (after which subsection)
    fig_insert = {
        'HFD-induced obesity and metabolic changes': [("Fig. 1", FIGURE_CAPTIONS["Fig. 1"])],
        'Diet-associated community differences persist': [("Fig. 2", FIGURE_CAPTIONS["Fig. 2"])],
        'Akkermansia expansion characterizes': [("Fig. 3", FIGURE_CAPTIONS["Fig. 3"])],
        'Lactococcus abundance correlates': [("Fig. 4", FIGURE_CAPTIONS["Fig. 4"])],
        'Plasmid analysis reveals': [("Fig. 5", FIGURE_CAPTIONS["Fig. 5"])],
        'Butyrate-synthesis capacity collapses': [("Fig. 6", FIGURE_CAPTIONS["Fig. 6"])],
    }

    current_subsection = None
    skip_comment = False
    i = 0

    while i < len(lines):
        line = lines[i].strip()
        i += 1

        if not line:
            continue

        if line.startswith('#'):
            continue

        # Check if this is a section head
        is_section = any(line.startswith(s) for s in section_heads)
        is_subsection = any(line.startswith(s) for s in subsection_heads)

        if is_section:
            if line == 'Results':
                story.append(PageBreak())
            elif line == 'Discussion':
                story.append(PageBreak())
            elif line == 'Methods':
                story.append(PageBreak())
            elif line == 'References':
                story.append(PageBreak())
            story.append(Paragraph(line, styles['SectionHead']))

        elif is_subsection:
            current_subsection = line
            story.append(Paragraph(line, styles['SubsectionHead']))

        elif line.startswith('Intrinsic, chromosomally encoded resistance dominates'):
            story.append(Paragraph(line, styles['ManuscriptTitle']))

        elif line.startswith('Jungman Kim'):
            story.append(Paragraph(line, styles['Authors']))

        elif line.startswith('1 ') or line.startswith('2 ') or line.startswith('*Corresponding') or line.startswith('Department') or line.startswith('Email') or line.startswith('ORCID'):
            story.append(Paragraph(line, styles['Affiliation']))

        elif line.startswith('Keywords:'):
            story.append(Paragraph(f"<b>{line}</b>", styles['BodyText2']))

        elif line.startswith('(') and line[1].isdigit() and ')' in line:
            # Group list items
            story.append(Paragraph(line, styles['BodyText2']))

        elif any(line.startswith(f"{n}.") for n in range(1, 45)):
            # Reference entries
            story.append(Paragraph(line, styles['RefStyle']))

        elif line.startswith('Additional file'):
            story.append(Paragraph(line, styles['RefStyle']))

        else:
            # Regular paragraph - escape XML chars
            safe_line = line.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
            story.append(Paragraph(safe_line, styles['BodyText2']))

        # Insert figures after relevant subsections end
        # Check if next line is empty (end of subsection content)
        if current_subsection and i < len(lines) and lines[i].strip() == '':
            # Check if the NEXT non-empty line is a new subsection/section
            j = i + 1
            while j < len(lines) and lines[j].strip() == '':
                j += 1
            if j < len(lines):
                next_line = lines[j].strip()
                next_is_head = (any(next_line.startswith(s) for s in section_heads) or
                                any(next_line.startswith(s) for s in subsection_heads))
                if next_is_head:
                    # Insert figures for this subsection
                    for sub_key, fig_list in fig_insert.items():
                        if current_subsection and current_subsection.startswith(sub_key):
                            story.append(Spacer(1, 12))
                            for fig_key, caption in fig_list:
                                block = []
                                if fig_key in FIGURES and os.path.exists(FIGURES[fig_key]):
                                    img = Image(FIGURES[fig_key])
                                    # Fit page width; leave room for the legend on the same page
                                    max_w = 6.5 * inch
                                    max_h = 6.0 * inch
                                    ratio = min(max_w / img.imageWidth, max_h / img.imageHeight)
                                    img._restrictSize(img.imageWidth * ratio, img.imageHeight * ratio)
                                    block += [img, Spacer(1, 4)]
                                else:
                                    raise FileNotFoundError(FIGURES.get(fig_key, fig_key))
                                if caption:
                                    safe_caption = caption.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
                                    m = re.match(r'(Fig\. \d+\.)(.*)', safe_caption)
                                    block.append(Paragraph(f"<b>{m.group(1)}</b>{m.group(2)}", styles['FigCaption']))
                                story.append(KeepTogether(block))
                                story.append(Spacer(1, 8))
                            # Remove from dict to avoid re-insertion
                            del fig_insert[sub_key]
                            break
                    current_subsection = None

    # Build PDF
    doc.build(story)
    print(f"PDF saved to: {OUTFILE}")
    print(f"Size: {os.path.getsize(OUTFILE) / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    build_pdf()
