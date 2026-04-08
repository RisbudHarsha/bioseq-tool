import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")

# -------------------------------------

def transcribe(dna):
    return dna.replace('T', 'U')

def translate(rna):
    codons = {
        'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
        'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
        'UAU':'Y','UAC':'Y','UAA':'','UAG':'',
        'UGU':'C','UGC':'C','UGA':'*','UGG':'W',

        'CUU':'L','CUC':'L','CUA':'L','CUG':'L',
        'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'CGU':'R','CGC':'R','CGA':'R','CGG':'R',

        'AUU':'I','AUC':'I','AUA':'I','AUG':'M',
        'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAU':'N','AAC':'N','AAA':'K','AAG':'K',
        'AGU':'S','AGC':'S','AGA':'R','AGG':'R',

        'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
        'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAU':'D','GAC':'D','GAA':'E','GAG':'E',
        'GGU':'G','GGC':'G','GGA':'G','GGG':'G'
    }

    protein = ''
    for i in range(0, len(rna)-2, 3):
        protein += codons.get(rna[i:i+3], 'X')
    return protein

def gc_content(seq):
    return (seq.count('G') + seq.count('C')) * 100 / len(seq) if seq else 0

# .......... ORF FINDER ..........

def find_orfs(seq):
    start = "ATG"
    stop = ["TAA", "TAG", "TGA"]
    orfs = []

    for frame in range(3):
        i = frame
        while i < len(seq)-2:
            codon = seq[i:i+3]

            if codon == start:
                start_pos = i

                for j in range(i+3, len(seq), 3):
                    next_codon = seq[j:j+3]

                    if next_codon in stop:
                        orfs.append((start_pos, j+3))
                        break

                i = j
            i += 3

    return orfs

# -------- HIGHLIGHT ORF --------

def highlight_orfs(seq, orfs):
    highlighted = ""

    for i in range(len(seq)):
        in_orf = False

        for start, end in orfs:
            if start <= i < end:
                in_orf = True
                break

        if in_orf:
            highlighted += f"<span style='color:red; font-weight:bold'>{seq[i]}</span>"
        else:
            highlighted += seq[i]

    return highlighted

# -------- UI --------

st.title("🧬 BioSeq Tool")

tab1, tab2, tab3, tab4 = st.tabs(["Home", "About", "Tool", "Team"])

# HOME
with tab1:
    st.write("Welcome to DNA analysis tool")
    st.write("Analyze DNA → RNA → Protein and find ORFs")

# ABOUT
with tab2:
    st.write("This tool performs basic bioinformatics analysis.")
    st.write("- DNA to RNA (Transcription)")
    st.write("- RNA to Protein (Translation)")
    st.write("- GC Content calculation")
    st.write("- Nucleotide composition graph")
    st.write("- ORF Finder with highlighting")

    st.write("Technologies Used:")
    st.write("- Python")
    st.write("- Streamlit")
    st.write("- Matplotlib")

    st.write("Users:")
    st.write("- Students and beginners in bioinformatics")

# TOOL
with tab3:
    seq = st.text_area("Enter DNA Sequence").upper()

    if st.button("Run Analysis"):
        if seq:

            # RNA
            rna = transcribe(seq)
            st.subheader("mRNA")
            st.write(rna)

            # Protein
            protein = translate(rna)
            st.subheader("Protein")
            st.write(protein)

            # ORFs
            st.subheader("ORFs")
            orfs = find_orfs(seq)

            if orfs:
                for i, (start, end) in enumerate(orfs):
                    st.write(f"ORF {i+1}: Start {start}, End {end}")
            else:
                st.write("No ORFs found")

            # Highlighted DNA
            st.subheader("Highlighted ORFs in DNA")
            colored_seq = highlight_orfs(seq, orfs)
            st.markdown(colored_seq, unsafe_allow_html=True)

            # Pie Chart
            st.subheader("Nucleotide Composition")
            counts = [seq.count('A'), seq.count('T'), seq.count('G'), seq.count('C')]
            labels = ['A','T','G','C']
            fig, ax = plt.subplots()
            ax.pie(counts, labels=labels, autopct='%1.1f%%')
            st.pyplot(fig)

           # GC Content 
            gc = gc_content(seq)
            st.subheader("GC Content")
            fig2, ax2 = plt.subplots()
            ax2.bar(['GC%'], [gc])
            st.pyplot(fig2)
            

            # DOWNLOAD
            result = f"DNA: {seq}\nRNA: {rna}\nProtein: {protein}\nGC: {gc:.2f}%\n\nORFs:\n"

            for i, (start, end) in enumerate(orfs):
                result += f"ORF {i+1}: {start}-{end}\n"

            st.download_button("Download Results", result)

        else:
            st.warning("Please enter DNA sequence")

# TEAM
with tab4:
    st.write("Name: Harshada Sudhir Risbud")
    st.write("Email: 3522511016@despu.edu.in")
    st.write("LinkedIn: your link")





