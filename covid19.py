import streamlit as st
import pandas as pd
import altair as alt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

"""
# Covid-19 Genome analysis
This is a simple analysis of SARS-COV-2 complete genome samples from Thailand and China. The data is downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202%20(SARS-CoV-2),%20taxid:2697049).

Tools
- BioPython
- Streamlit
"""

covid19th = SeqIO.read('./data/thailand/covid19-dna-complete.fasta', "fasta")
covid19ch = SeqIO.read('./data/china/covid19-dna-complete.fasta', "fasta")
covid_DNA_th = covid19th.seq
covid_DNA_ch = covid19ch.seq
dna_th_200 = covid_DNA_th[:200]
dna_ch_200 = covid_DNA_ch[:200]

checkbox1 = st.sidebar.checkbox('View full TH DNA sequence')
checkbox2 = st.sidebar.checkbox('View full CH DNA sequence')

if checkbox1:
    dnaLength1 = covid_DNA_th
else:
    dnaLength1 = dna_th_200

if checkbox2:
    dnaLength2 = covid_DNA_ch
else:
    dnaLength2 = dna_ch_200

'''
## DNA sequence
'''

sequence = st.text_area("SARS-CoV-2 Thailand DNA Seqeunce", dnaLength1, height=30)
sequence = st.text_area("SARS-CoV-2 China DNA Seqeunce", dnaLength2, height=30)
st.write(f'Thailand DNA Seqeunce consists of {len(covid19th)} genetic bases or letters.')
st.write(f'China DNA Seqeunce consists of {len(covid19ch)} genetic bases or letters.')

'''
## Nucleotides
### Nucleotides frequency
Count the nucleotides frequency in the DNA
'''
DNA_th = covid_DNA_th
DNA_ch = covid_DNA_ch
nucleotides_th={}
nucleotides_ch={}

# get nucleotides frequency
def getNucleotides(dna, nucleotides):
    for n in dna:
        if n in nucleotides:
            nucleotides[n] += 1
        else:
            nucleotides[n] =  1

getNucleotides(DNA_th, nucleotides_th)
getNucleotides(DNA_ch, nucleotides_ch)

col1, col2 = st.beta_columns(2)

col1.write('Thailand')
df1 = pd.DataFrame.from_dict(nucleotides_th, orient='index')
df1 = df1.rename({0: 'Nucleotide Frequency'}, axis='columns')
df1.reset_index(inplace=True)
df1 = df1.rename(columns = {'index': 'Nucleotide'})
col1.dataframe(df1)

col2.write('China')
df2 = pd.DataFrame.from_dict(nucleotides_ch, orient='index')
df2 = df2.rename({0: 'Nucleotide Frequency'}, axis='columns')
df2.reset_index(inplace=True)
df2 = df2.rename(columns = {'index': 'Nucleotide'})
col2.dataframe(df2)

'''
### Nucleotide chart
'''
st.write('Thailand')
p1 = alt.Chart(df1).mark_bar().encode(
    x=alt.X('Nucleotide', sort='ascending'),
    y='Nucleotide Frequency',
    tooltip=['Nucleotide','Nucleotide Frequency']
)
st.altair_chart(p1, use_container_width=True)

st.write('China')
p2 = alt.Chart(df2).mark_bar().encode(
    x=alt.X('Nucleotide', sort='ascending'),
    y='Nucleotide Frequency',
    tooltip=['Nucleotide','Nucleotide Frequency']
)
st.altair_chart(p2, use_container_width=True)

'''
## mRNA
### Transcription
Basically the mRNA is a copy of our DNA.
However, in RNA, a base called uracil (U) replaces thymine (T) as the complementary
nucleotide to adenine (that's the only difference, T is replaced by U).
'''
covid_mRNA_th = covid_DNA_th.transcribe()
covid_mRNA_ch = covid_DNA_ch.transcribe()
sequence = st.text_area("Thailand mRNA Seqeunce", covid_mRNA_th[:100], height=25)
sequence = st.text_area("China mRNA Seqeunce", covid_mRNA_ch[:100], height=25)

'''
In order to see the difference, we align the covid-19 DNA and mRNA sequences.
We can see that the mRNA is an identical copy with the T base replaced by U.
In this example we show COVID-19 from Thailand.
'''

st.markdown(f'Covid-19 TH DNA: `{covid_DNA_th[:50]}`')
st.markdown(f'Covid-19 TH RNA: `{covid_mRNA_th[:50]}`')

'''
## Amino acids
### Translation
Translation is the process that takes the information passed from DNA as
messenger RNA and turns it into a series of amino acids.
'''

def pad_seq(sequence):
    """ Pad sequence to multiple of 3 with N"""
    remainder = len(sequence) % 3
    return sequence if remainder == 0 else sequence + Seq('N' * (3 - remainder))

aa_th = pad_seq(covid_mRNA_th).translate()
aa_ch= pad_seq(covid_mRNA_ch).translate()
proteins_th = aa_th.split('*')
proteins_ch = aa_ch.split('*')
st.text_area("Thailand mRNA Translation", aa_th[:100], height=25)
st.text_area("China mRNA Translation", aa_ch[:100], height=25)
st.markdown(f"TH Covid-19's genome has {len(aa_th)} amino acids and {len(proteins_th)} amino acid chains")
st.markdown(f"CH Covid-19's genome has {len(aa_ch)} amino acids and {len(proteins_ch)} amino acid chains")

for i in proteins_th[:]:
    if len(i) < 50:
        proteins_th.remove(i)

for i in proteins_ch[:]:
    if len(i) < 50:
        proteins_ch.remove(i)

col1, col2 = st.beta_columns(2)
protein_clean_th = [str(i) for i in proteins_th]
df3 = pd.DataFrame({'amino_acids':protein_clean_th})
df3['count'] = df3['amino_acids'].str.len()
df3.head()
col1.write(df3)
col1.markdown(f'{len(proteins_th)} proteins with more than 50 amino acids from the Thailand covid-19 genome')

protein_clean_ch = [str(i) for i in proteins_ch]
d4 = pd.DataFrame({'amino_acids':protein_clean_ch})
d4['count'] = d4['amino_acids'].str.len()
d4.head()
col2.write(d4)
col2.markdown(f'{len(proteins_ch)} proteins with more than 50 amino acids from the China covid-19 genome')

'''
## Comparing Thailand and China strain
Using pairwise2 to compare sequence alignment in order to identify similarities.
'''

alignment_score = pairwise2.align.globalxx(covid_DNA_th, covid_DNA_ch, one_alignment_only=True, score_only=True)
st.write('COV-TH/COV-CH Similarity (%):', alignment_score / len(covid_DNA_th) * 100)