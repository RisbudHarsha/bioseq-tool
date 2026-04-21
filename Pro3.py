import streamlit as st
import matplotlib.pyplot as plt

# -------- FUNCTIONS --------

def transcribe(dna):
    return dna.replace('T', 'U')

def translate(rna):
    codons = {
        'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
        'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
        'UAU':'Y','UAC':'Y','UAA':'*','UAG':'*',
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


# -------- NEW: PROTEIN CHARACTERIZATION --------

def protein_characterization(protein):
    aa_count = {}
    for aa in protein:
        aa_count[aa] = aa_count.get(aa, 0) + 1

    length = len(protein)

    # Approx molecular weight (110 Da per amino acid)
    molecular_weight = length * 110

    # Hydrophobic amino acids
    hydrophobic = ['A','V','I','L','M','F','W','Y']
    hydro_count = sum([protein.count(aa) for aa in hydrophobic])

    return length, molecular_weight, aa_count, hydro_count


# -------- UI --------

st.title("🧬 BioSeq Tool")

menu = st.sidebar.radio("Menu", ["Home","About","Tool","Team"])


# -------- HOME --------
if menu == "Home":
    st.write("Welcome to DNA analysis tool")
    st.write("You can convert DNA to RNA and Protein and analyze sequence.")


# -------- ABOUT --------
elif menu == "About":
    st.title("About Tool")

    st.write("This is a simple bioinformatics tool built using Python and Streamlit.")

    st.write("Purpose:")
    st.write("- To perform basic DNA sequence analysis easily")

    st.write("Features:")
    st.write("- DNA to RNA (Transcription)")
    st.write("- RNA to Protein (Translation)")
    st.write("- GC Content calculation")
    st.write("- Nucleotide composition graph")
    st.write("- Protein characterization")

    st.write("Technologies Used:")
    st.write("- Python")
    st.write("- Streamlit")
    st.write("- Matplotlib")

    st.write("Users:")
    st.write("- Students and beginners in bioinformatics")


# -------- TOOL --------
elif menu == "Tool":

    seq = st.text_area("Enter DNA Sequence").upper()

    if st.button("Run"):
        if seq:

            # Transcription
            rna = transcribe(seq)
            st.subheader("mRNA")
            st.write(rna)

            # Translation
            protein = translate(rna)
            st.subheader("Protein")
            st.write(protein)

            # -------- Protein Characterization --------
            length, mw, aa_count, hydro = protein_characterization(protein)

            st.subheader("Protein Characterization")
            st.write(f"Length: {length} amino acids")
            st.write(f"Approx Molecular Weight: {mw} Da")
            st.write(f"Hydrophobic Residues: {hydro}")

            st.write("Amino Acid Composition:")
            st.write(aa_count)

            # Amino Acid Graph
            st.subheader("Amino Acid Composition Graph")
            fig2, ax2 = plt.subplots()
            ax2.bar(aa_count.keys(), aa_count.values())
            st.pyplot(fig2)

            # -------- Nucleotide Composition --------
            st.subheader("Nucleotide Composition")
            counts = [seq.count('A'), seq.count('T'), seq.count('G'), seq.count('C')]
            labels = ['A','T','G','C']
            fig, ax = plt.subplots()
            ax.pie(counts, labels=labels, autopct='%1.1f%%')
            st.pyplot(fig)

            # -------- GC Content --------
            gc = gc_content(seq)
            st.subheader("GC Content")
            st.write(f"{gc:.2f}%")

            # -------- Download --------
            result = f"DNA: {seq}\n"
            result += f"RNA: {rna}\n"
            result += f"Protein: {protein}\n"
            result += f"Protein Length: {length}\n"
            result += f"Molecular Weight: {mw} Da\n"
            result += f"Hydrophobic Residues: {hydro}\n"
            result += f"GC Content: {gc:.2f}%\n"

            st.download_button("Download Results", result, file_name="result.txt")

        else:
            st.write("Please enter sequence")


# -------- TEAM --------
elif menu == "Team":
    st.title("Team")
    st.write("Name: Harshada Risbud")
    st.write("Email: 3522511016@despu.edu.in")





