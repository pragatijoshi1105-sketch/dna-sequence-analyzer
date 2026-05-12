import streamlit as st
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ─── Page Configuration ───────────────────────────────────────────────────────
st.set_page_config(
    page_title="DNA Sequence Analyzer",
    page_icon="🧬",
    layout="wide"
)

# ─── App Title ─────────────────────────────────────────────────────────────────
st.title("🧬 DNA Sequence Analyzer")
st.markdown("### A Bioinformatics Tool for Sequence Analysis")
st.markdown("---")

# ─── Sidebar ───────────────────────────────────────────────────────────────────
st.sidebar.title("⚙️ Settings")
st.sidebar.markdown("Paste your DNA sequence and explore its properties.")

example_sequences = {
    "None": "",
    "Short Example": "ATGCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCTAGC",
    "Human BRCA1 Fragment": "ATGCGATCGATCGAATTCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    "GC Rich Sequence": "GCGCGCGCGCGCGCATGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
}

selected_example = st.sidebar.selectbox("Load an example sequence:", list(example_sequences.keys()))

# ─── Input Area ────────────────────────────────────────────────────────────────
st.subheader("📥 Input Your DNA Sequence")

if selected_example != "None":
    default_seq = example_sequences[selected_example]
else:
    default_seq = ""

raw_input = st.text_area(
    "Paste your DNA sequence below (A, T, G, C only):",
    value=default_seq,
    height=150,
    placeholder="Example: ATGCGATCGATCG..."
)

analyze_button = st.button("🔬 Analyze Sequence", type="primary")

# ─── Analysis Logic ────────────────────────────────────────────────────────────
if analyze_button:

    # Clean the input
    sequence = raw_input.upper().strip().replace(" ", "").replace("\n", "")

    # Validate input
    if len(sequence) == 0:
        st.error("❌ Please enter a DNA sequence first.")

    elif not all(base in "ATGC" for base in sequence):
        invalid = set(sequence) - set("ATGC")
        st.error(f"❌ Invalid characters found: {invalid}. Please use only A, T, G, C.")

    else:
        # Create BioPython Seq object
        dna = Seq(sequence)

        st.success(f"✅ Valid sequence detected! Analyzing {len(sequence)} base pairs...")
        st.markdown("---")

        # ── Section 1: Basic Statistics ──────────────────────────────────────
        st.subheader("📊 Basic Sequence Statistics")

        col1, col2, col3, col4 = st.columns(4)

        count_a = sequence.count("A")
        count_t = sequence.count("T")
        count_g = sequence.count("G")
        count_c = sequence.count("C")
        gc = round(gc_fraction(dna) * 100, 2)
        at = round(100 - gc, 2)

        col1.metric("🔵 Adenine (A)", count_a)
        col2.metric("🔴 Thymine (T)", count_t)
        col3.metric("🟢 Guanine (G)", count_g)
        col4.metric("🟡 Cytosine (C)", count_c)

        col5, col6, col7 = st.columns(3)
        col5.metric("📏 Sequence Length", f"{len(sequence)} bp")
        col6.metric("🧪 GC Content", f"{gc}%")
        col7.metric("🧪 AT Content", f"{at}%")

        st.markdown("---")

        # ── Section 2: Charts ─────────────────────────────────────────────────
        st.subheader("📈 Visual Analysis")

        chart_col1, chart_col2 = st.columns(2)

        with chart_col1:
            st.markdown("*Nucleotide Composition*")
            fig1, ax1 = plt.subplots(figsize=(5, 4))
            bases = ["Adenine (A)", "Thymine (T)", "Guanine (G)", "Cytosine (C)"]
            counts = [count_a, count_t, count_g, count_c]
            colors = ["#4A90D9", "#E74C3C", "#2ECC71", "#F1C40F"]
            bars = ax1.bar(bases, counts, color=colors, edgecolor="white", linewidth=1.5)
            ax1.set_ylabel("Count", fontsize=11)
            ax1.set_title("Nucleotide Frequency", fontsize=13, fontweight="bold")
            ax1.set_facecolor("#0E1117")
            fig1.patch.set_facecolor("#0E1117")
            ax1.tick_params(colors="white")
            ax1.yaxis.label.set_color("white")
            ax1.title.set_color("white")
            for spine in ax1.spines.values():
                spine.set_edgecolor("#444")
            for bar, count in zip(bars, counts):
                ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                        str(count), ha="center", va="bottom", color="white", fontweight="bold")
            st.pyplot(fig1)

        with chart_col2:
            st.markdown("*GC vs AT Content*")
            fig2, ax2 = plt.subplots(figsize=(5, 4))
            labels = ["GC Content", "AT Content"]
            values = [gc, at]
            colors2 = ["#2ECC71", "#4A90D9"]
            wedges, texts, autotexts = ax2.pie(
                values, labels=labels, colors=colors2,
                autopct="%1.1f%%", startangle=90,
                textprops={"color": "white", "fontsize": 11}
            )
            for autotext in autotexts:
                autotext.set_color("white")
                autotext.set_fontweight("bold")
            ax2.set_title("GC vs AT Ratio", fontsize=13, fontweight="bold", color="white")
            fig2.patch.set_facecolor("#0E1117")
            st.pyplot(fig2)

        st.markdown("---")

        # ── Section 3: Molecular Biology Operations ───────────────────────────
        st.subheader("🔬 Molecular Biology Operations")

        mol_col1, mol_col2 = st.columns(2)

        with mol_col1:
            st.markdown("*Complement Strand*")
            st.code(str(dna.complement()), language=None)

            st.markdown("*Reverse Complement*")
            st.code(str(dna.reverse_complement()), language=None)

        with mol_col2:
            st.markdown("*mRNA Transcript (DNA → RNA)*")
            st.code(str(dna.transcribe()), language=None)

            st.markdown("*Protein Translation (first ORF)*")
            try:
                protein = str(dna.translate(to_stop=True))
                if len(protein) == 0:
                    st.info("No complete ORF found in this reading frame.")
                else:
                    st.code(protein, language=None)
            except Exception:
                st.info("Could not translate — sequence may not start with ATG.")

        st.markdown("---")

        # ── Section 4: Codon Analysis ─────────────────────────────────────────
        st.subheader("🧩 Codon Analysis")

        codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
        codon_counts = {}
        for codon in codons:
            if len(codon) == 3:
                codon_counts[codon] = codon_counts.get(codon, 0) + 1

        if codon_counts:
            sorted_codons = sorted(codon_counts.items(), key=lambda x: x[1], reverse=True)
            top_codons = sorted_codons[:10]

            fig3, ax3 = plt.subplots(figsize=(10, 4))
            codon_labels = [c[0] for c in top_codons]
            codon_values = [c[1] for c in top_codons]
            ax3.bar(codon_labels, codon_values, color="#9B59B6", edgecolor="white", linewidth=1.2)
            ax3.set_xlabel("Codon", fontsize=11, color="white")
            ax3.set_ylabel("Frequency", fontsize=11, color="white")
            ax3.set_title("Top 10 Most Frequent Codons", fontsize=13, fontweight="bold", color="white")
            ax3.set_facecolor("#0E1117")
            fig3.patch.set_facecolor("#0E1117")
            ax3.tick_params(colors="white")
            for spine in ax3.spines.values():
                spine.set_edgecolor("#444")
            st.pyplot(fig3)

            st.markdown(f"*Total codons found:* {len(codons)}")

        st.markdown("---")

        # ── Section 5: Melting Temperature ───────────────────────────────────
        st.subheader("🌡️ Estimated Melting Temperature (Tm)")

        if len(sequence) < 14:
            tm = 2 * (count_a + count_t) + 4 * (count_g + count_c)
            formula_used = "Wallace Rule (for short sequences < 14 bp)"
        else:
            tm = 64.9 + 41 * (count_g + count_c - 16.4) / len(sequence)
            formula_used = "Marmur-Chargaff Formula (for longer sequences)"

        st.metric("🌡️ Melting Temperature (Tm)", f"{round(tm, 2)} °C")
        st.caption(f"Formula used: {formula_used}")

        st.markdown("---")

        # ── Footer ────────────────────────────────────────────────────────────
        st.markdown(
            """
            <div style='text-align: center; color: grey; font-size: 13px;'>
            Built with ❤️ using Python, BioPython & Streamlit | DNA Sequence Analyzer v1.0
            </div>
            """,
            unsafe_allow_html=True
        )