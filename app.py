import streamlit as st
import matplotlib.pyplot as plt
from collections import Counter

# Codon usage tables (frequency per 1000 codons)
CODON_TABLES = {
    'human': {
        'TTT': 17.6, 'TTC': 20.3, 'TTA': 7.7, 'TTG': 12.9,
        'CTT': 13.2, 'CTC': 19.6, 'CTA': 7.2, 'CTG': 39.6,
        'ATT': 16.0, 'ATC': 20.8, 'ATA': 7.5, 'ATG': 22.0,
        'GTT': 11.0, 'GTC': 14.5, 'GTA': 7.1, 'GTG': 28.1,
        'TCT': 15.2, 'TCC': 17.7, 'TCA': 12.2, 'TCG': 4.4,
        'CCT': 17.5, 'CCC': 19.8, 'CCA': 16.9, 'CCG': 6.9,
        'ACT': 13.1, 'ACC': 18.9, 'ACA': 15.1, 'ACG': 6.1,
        'GCT': 18.4, 'GCC': 27.7, 'GCA': 15.8, 'GCG': 7.4,
        'TAT': 12.2, 'TAC': 15.3, 'TAA': 1.0, 'TAG': 0.8,
        'CAT': 10.9, 'CAC': 15.1, 'CAA': 12.3, 'CAG': 34.2,
        'AAT': 17.0, 'AAC': 19.1, 'AAA': 24.4, 'AAG': 31.9,
        'GAT': 21.8, 'GAC': 25.1, 'GAA': 29.0, 'GAG': 39.6,
        'TGT': 10.6, 'TGC': 12.6, 'TGA': 1.6, 'TGG': 13.2,
        'CGT': 4.5, 'CGC': 10.4, 'CGA': 6.2, 'CGG': 11.4,
        'AGT': 12.1, 'AGC': 19.5, 'AGA': 12.2, 'AGG': 12.0,
        'GGT': 10.8, 'GGC': 22.2, 'GGA': 16.5, 'GGG': 16.5
    },
    'ecoli': {
        'TTT': 22.0, 'TTC': 16.0, 'TTA': 13.0, 'TTG': 13.0,
        'CTT': 11.0, 'CTC': 10.0, 'CTA': 4.0, 'CTG': 50.0,
        'ATT': 30.0, 'ATC': 25.0, 'ATA': 7.0, 'ATG': 27.0,
        'GTT': 26.0, 'GTC': 15.0, 'GTA': 15.0, 'GTG': 26.0,
        'TCT': 9.0, 'TCC': 9.0, 'TCA': 7.0, 'TCG': 9.0,
        'CCT': 7.0, 'CCC': 5.0, 'CCA': 8.0, 'CCG': 23.0,
        'ACT': 13.0, 'ACC': 23.0, 'ACA': 7.0, 'ACG': 14.0,
        'GCT': 15.0, 'GCC': 25.0, 'GCA': 21.0, 'GCG': 33.0,
        'TAT': 16.0, 'TAC': 12.0, 'TAA': 2.0, 'TAG': 0.2,
        'CAT': 13.0, 'CAC': 10.0, 'CAA': 15.0, 'CAG': 29.0,
        'AAT': 19.0, 'AAC': 22.0, 'AAA': 33.0, 'AAG': 10.0,
        'GAT': 32.0, 'GAC': 19.0, 'GAA': 40.0, 'GAG': 18.0,
        'TGT': 5.0, 'TGC': 6.0, 'TGA': 1.0, 'TGG': 15.0,
        'CGT': 21.0, 'CGC': 21.0, 'CGA': 4.0, 'CGG': 5.0,
        'AGT': 9.0, 'AGC': 16.0, 'AGA': 2.0, 'AGG': 1.0,
        'GGT': 24.0, 'GGC': 29.0, 'GGA': 8.0, 'GGG': 11.0
    },
    'cho': {
        'TTT': 18.0, 'TTC': 20.0, 'TTA': 8.0, 'TTG': 13.0,
        'CTT': 13.0, 'CTC': 19.0, 'CTA': 7.0, 'CTG': 40.0,
        'ATT': 17.0, 'ATC': 21.0, 'ATA': 8.0, 'ATG': 22.0,
        'GTT': 12.0, 'GTC': 15.0, 'GTA': 7.0, 'GTG': 27.0,
        'TCT': 15.0, 'TCC': 18.0, 'TCA': 12.0, 'TCG': 5.0,
        'CCT': 17.0, 'CCC': 20.0, 'CCA': 17.0, 'CCG': 7.0,
        'ACT': 13.0, 'ACC': 19.0, 'ACA': 15.0, 'ACG': 6.0,
        'GCT': 19.0, 'GCC': 28.0, 'GCA': 16.0, 'GCG': 7.0,
        'TAT': 12.0, 'TAC': 15.0, 'TAA': 1.0, 'TAG': 0.8,
        'CAT': 11.0, 'CAC': 15.0, 'CAA': 12.0, 'CAG': 34.0,
        'AAT': 17.0, 'AAC': 19.0, 'AAA': 25.0, 'AAG': 32.0,
        'GAT': 22.0, 'GAC': 25.0, 'GAA': 29.0, 'GAG': 40.0,
        'TGT': 11.0, 'TGC': 13.0, 'TGA': 1.5, 'TGG': 13.0,
        'CGT': 5.0, 'CGC': 10.0, 'CGA': 6.0, 'CGG': 11.0,
        'AGT': 12.0, 'AGC': 19.0, 'AGA': 12.0, 'AGG': 12.0,
        'GGT': 11.0, 'GGC': 22.0, 'GGA': 17.0, 'GGG': 17.0
    }
}

GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def translate_dna(dna):
    """Translate DNA sequence to protein."""
    protein = ''
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        protein += GENETIC_CODE.get(codon, 'X')
    return protein


def optimize_sequence(protein_seq, codon_table):
    """Optimize protein sequence for given organism."""
    # Group codons by amino acid
    aa_to_codons = {}
    for codon, aa in GENETIC_CODE.items():
        if aa not in aa_to_codons:
            aa_to_codons[aa] = []
        aa_to_codons[aa].append(codon)

    # Build optimized DNA
    optimized = ''
    for aa in protein_seq:
        codons = aa_to_codons.get(aa, [])
        if not codons:
            continue
        # Select codon with highest frequency
        best_codon = max(codons, key=lambda c: codon_table.get(c, 0))
        optimized += best_codon

    return optimized


def find_optimized_codons(input_dna, optimized_dna):
    """Find which codons were changed during optimization."""
    if not input_dna:
        return []

    changed_positions = []
    for i in range(0, min(len(input_dna), len(optimized_dna)) - 2, 3):
        input_codon = input_dna[i:i+3]
        opt_codon = optimized_dna[i:i+3]
        if input_codon != opt_codon:
            changed_positions.append(i // 3)  # Codon position

    return changed_positions


def highlight_sequence(dna, changed_positions):
    """Create HTML with highlighted codons."""
    html = '<div style="font-family: monospace; line-height: 2; word-wrap: break-word;">'

    for i in range(0, len(dna), 3):
        codon_pos = i // 3
        codon = dna[i:i+3]

        if codon_pos in changed_positions:
            html += f'<span style="background-color: #10b981; color: white; padding: 2px 4px; margin: 1px; border-radius: 3px;">{codon}</span>'
        else:
            html += f'<span style="padding: 2px 4px; margin: 1px;">{codon}</span>'

    html += '</div>'
    return html


def analyze_codon_usage(dna):
    """Count codon usage in DNA sequence."""
    codons = [dna[i:i+3] for i in range(0, len(dna) - 2, 3)]
    return Counter(codons)


def plot_codon_usage(usage, title, color):
    """Plot codon usage as bar chart."""
    # Get top 15 codons
    top_codons = usage.most_common(15)
    if not top_codons:
        return None

    codons, counts = zip(*top_codons)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(codons, counts, color=color)
    ax.set_xlabel('Codon')
    ax.set_ylabel('Count')
    ax.set_title(title)
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()

    return fig


# Streamlit UI
st.set_page_config(page_title="Codon Optimizer", layout="wide")

st.title("🧬 Codon Optimizer")

# Input section
col1, col2 = st.columns(2)

with col1:
    input_type = st.selectbox("Sequence Type", ["DNA", "Protein"])

with col2:
    organism = st.selectbox("Organism", ["human", "ecoli", "cho"],
                            format_func=lambda x: {"human": "Human", "ecoli": "E. coli", "cho": "CHO Cells"}[x])

sequence = st.text_area(
    "Enter Sequence",
    placeholder="ATGGCA..." if input_type == "DNA" else "MAKL...",
    height=150
)

if st.button("Optimize", type="primary", use_container_width=True):
    if not sequence:
        st.error("Please enter a sequence!")
    else:
        try:
            # Clean and validate input
            clean_seq = sequence.upper().replace(" ", "").replace("\n", "")

            if input_type == "Protein":
                if not all(c in "ACDEFGHIKLMNPQRSTVWY*" for c in clean_seq):
                    st.error("Invalid protein sequence!")
                    st.stop()
                protein = clean_seq
                input_dna = None
            else:
                if not all(c in "ATGC" for c in clean_seq):
                    st.error("Invalid DNA sequence! Only A, T, G, C allowed.")
                    st.stop()
                if len(clean_seq) % 3 != 0:
                    st.error(
                        f"DNA length must be divisible by 3! Current length: {len(clean_seq)} (remainder: {len(clean_seq) % 3})")
                    st.stop()
                input_dna = clean_seq
                protein = translate_dna(clean_seq)

            # Optimize
            codon_table = CODON_TABLES[organism]
            optimized_dna = optimize_sequence(protein, codon_table)

            # Find changed codons
            changed_positions = find_optimized_codons(input_dna, optimized_dna)

            # Analyze codon usage
            if input_dna:
                input_usage = analyze_codon_usage(input_dna)
            output_usage = analyze_codon_usage(optimized_dna)

            # Display results
            st.success("Optimization successful!")

            st.subheader("Results")

            if input_dna:
                st.markdown("**Original DNA:**")
                st.code(input_dna, language=None)

            st.markdown("**Protein:**")
            st.code(protein, language=None)

            st.markdown("**Optimized DNA:** (🟢 = changed codons)")
            highlighted_html = highlight_sequence(
                optimized_dna, changed_positions)
            st.markdown(highlighted_html, unsafe_allow_html=True)

            # Show statistics
            if input_dna:
                total_codons = len(input_dna) // 3
                num_changed = len(changed_positions)
                pct_changed = (num_changed / total_codons *
                               100) if total_codons > 0 else 0

                col1, col2, col3 = st.columns(3)
                col1.metric("Total Codons", total_codons)
                col2.metric("Optimized Codons", num_changed)
                col3.metric("Change Rate", f"{pct_changed:.1f}%")

            # Visualization
            st.subheader("Codon Usage Visualization")

            if input_dna:
                col1, col2 = st.columns(2)

                with col1:
                    st.markdown("**Original**")
                    fig1 = plot_codon_usage(
                        input_usage, "Original Codon Usage", "#6366f1")
                    if fig1:
                        st.pyplot(fig1)

                with col2:
                    st.markdown("**Optimized**")
                    fig2 = plot_codon_usage(
                        output_usage, "Optimized Codon Usage", "#10b981")
                    if fig2:
                        st.pyplot(fig2)
            else:
                st.markdown("**Optimized**")
                fig = plot_codon_usage(
                    output_usage, "Optimized Codon Usage", "#10b981")
                if fig:
                    st.pyplot(fig)

        except Exception as e:
            st.error(f"Error: {str(e)}")

# Sidebar info
with st.sidebar:
    st.header("ℹ️ Information")
    st.markdown("""
    **Codon Optimizer** selects the preferred codon for each amino acid in the target organism.
    
    **Organisms:**
    - Human
    - E. coli
    - CHO Cells
    
    **Data Source:**
    Codon frequencies based on Kazusa Database
    
    **Legend:**
    - 🟢 Green = Optimized codon (changed from original)
    - White = Unchanged codon
    """)

    st.markdown("---")
    st.markdown("**Example DNA:**")
    st.code("ATGGCATTTAAACGTGAAGAT", language=None)
    st.markdown("**Example Protein:**")
    st.code("MAKLFG", language=None)
