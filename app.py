import streamlit as st
from Bio import SeqIO
import io
import pandas as pd
import matplotlib.pyplot as plt

st.title("Antibiotic Resistance Gene Detection")

uploaded_file = st.file_uploader("Upload Genome FASTA", type=["fasta","fa","fna"])

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return (g+c)/len(seq)*100


def load_resistance_genes():
    genes = []
    for record in SeqIO.parse("resistance_genes.fasta","fasta"):
        genes.append((record.id,str(record.seq)))
    return genes


def search_resistance(genome_seq,resistance_db):
    detected = []

    for gene_id,gene_seq in resistance_db:

        if gene_seq[:30] in genome_seq:
            detected.append(gene_id)

    return detected


def extract_gene_name(gene_id):
    if "|" in gene_id:
        return gene_id.split("|")[-1]
    return gene_id


if uploaded_file:

    genome_record = SeqIO.read(io.TextIOWrapper(uploaded_file, encoding="utf-8"), "fasta")
    genome_seq = str(genome_record.seq)

    st.success("Genome uploaded successfully")

    st.write("Genome length:",len(genome_seq))
    st.write("GC Content:",gc_content(genome_seq))

    resistance_db = load_resistance_genes()

    # Searching message placeholder
    search_msg = st.empty()
    search_msg.markdown("**Searching for Antibiotic Resistance Genes...**")

    detected_genes = search_resistance(genome_seq,resistance_db)

    # Remove searching message after results
    search_msg.empty()

    # Prediction logic
    if len(detected_genes)==0:
        prediction = "Sensitive"
    else:
        prediction = "Resistant"

    st.write("Predicted Resistance:", prediction)

    if len(detected_genes)==0:

        st.warning("No resistance genes detected")

    else:

        st.subheader("Detected Resistance Genes")

        df = pd.DataFrame(detected_genes,columns=["Detected Resistance Gene"])
        st.dataframe(df)

        # Load ARO index
        aro = pd.read_csv("aro_index.tsv",sep="\t")

        antibiotics = []
        mechanisms = []
        gene_names = []

        for gene in detected_genes:

            gene_name = extract_gene_name(gene)
            gene_names.append(gene_name)

            aro_id = None

            if "ARO:" in gene:
                aro_id = gene.split("ARO:")[1].split("|")[0]

            if aro_id:

                row = aro[aro["ARO Accession"]=="ARO:"+aro_id]

                if len(row)>0:

                    antibiotics.append(row["Drug Class"].values[0])
                    mechanisms.append(row["Resistance Mechanism"].values[0])

                else:

                    antibiotics.append("Unknown")
                    mechanisms.append("Unknown")

        result_df = pd.DataFrame({

        "Gene":gene_names,
        "Antibiotic":antibiotics,
        "Mechanism":mechanisms

        })

        st.subheader("Associated Antibiotics")

        st.dataframe(result_df)

        gene_counts = pd.Series(gene_names).value_counts()

        fig,ax = plt.subplots()

        gene_counts.plot(kind="bar",ax=ax)

        ax.set_xlabel("Resistance Genes")
        ax.set_ylabel("Count")

        st.pyplot(fig)
