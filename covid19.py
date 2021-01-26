import streamlit as st
import pandas as pd
import altair as alt
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import pairwise2

"""
# Covid-19 Genome analysis
"""
covid19th = SeqIO.read('./data/thailand/covid19-dna-complete.fasta', "fasta")
covid19ch = SeqIO.read('./data/china/covid19-dna-complete.fasta', "fasta")
st.write(f'The genome of the virus causing Covid-19 TH (known as SARS-CoV-2) consists of {len(covid19th)} genetic bases or letters.')
st.write(f'The genome of the virus causing Covid-19 CH (known as SARS-CoV-2) consists of {len(covid19ch)} genetic bases or letters.')

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

sequence = st.text_area("SARS-CoV-2 Thailand DNA Seqeunce", dnaLength1, height=30)
sequence = st.text_area("SARS-CoV-2 China DNA Seqeunce", dnaLength2, height=30)


'''
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
st.json(nucleotides_th)
st.json(nucleotides_ch)

'''
### Data frame
Thailand
'''
df1 = pd.DataFrame.from_dict(nucleotides_th, orient='index')
df1 = df1.rename({0: 'Nucleotide Frequency'}, axis='columns')
df1.reset_index(inplace=True)
df1 = df1.rename(columns = {'index': 'Nucleotide'})
st.write(df1)
'''
China
'''
df2 = pd.DataFrame.from_dict(nucleotides_ch, orient='index')
df2 = df2.rename({0: 'Nucleotide Frequency'}, axis='columns')
df2.reset_index(inplace=True)
df2 = df2.rename(columns = {'index': 'Nucleotide'})
st.write(df2)

'''
### Nucleotide chart
Thailand
'''
p1 = alt.Chart(df1).mark_bar().encode(
    x=alt.X('Nucleotide', sort='ascending'),
    y='Nucleotide Frequency',
    tooltip=['Nucleotide','Nucleotide Frequency']
)
st.altair_chart(p1, use_container_width=True)

'''
China
'''
p2 = alt.Chart(df2).mark_bar().encode(
    x=alt.X('Nucleotide', sort='ascending'),
    y='Nucleotide Frequency',
    tooltip=['Nucleotide','Nucleotide Frequency']
)
st.altair_chart(p2, use_container_width=True)

'''
### Transcription
Basically the mRNA is a copy of our DNA.
However, in RNA, a base called uracil (U) replaces thymine (T) as the complementary
nucleotide to adenine (that's the only difference, T is replaced by U).
'''
covid_mRNA = covid_DNA_th.transcribe()
sequence = st.text_area("mRNA Seqeunce", covid_mRNA[:100], height=25)

'''
In order to see the difference, we align the covid-19 DNA and mRNA sequences.
We can see that the mRNA is an identical copy with the T base replaced by U
'''

st.markdown(f'Covid-19 DNA: `{covid_DNA_th[:50]}`')
st.markdown(f'Covid-19 RNA: `{covid_mRNA[:50]}`')

'''
### Translation
Translation is the process that takes the information passed from DNA as
messenger RNA and turns it into a series of amino acids.
'''
covid_aa = covid_mRNA.translate()
sequence = st.text_area("mRNA Translation", covid_aa[:100], height=25)
st.markdown(f"Covid-19's genome has {len(covid_aa)} amino acids")

Proteins = covid_aa.split('*')
st.markdown(f'We have {len(Proteins)} amino acid chains in the covid-19 genome')

for i in Proteins[:]:
    if len(i) < 20:
        Proteins.remove(i)

st.markdown(f'We have {len(Proteins)} proteins with  more than 20 amino acids in the covid-19 genome')

for i in Proteins[:]:
    if len(i) < 50:
        Proteins.remove(i)

st.markdown(f'We have {len(Proteins)} proteins with  more than 50 amino acids in the covid-19 genome')


protein_clean = [str(i) for i in Proteins]
df3 = pd.DataFrame({'amino_acids':protein_clean})
df3['count'] = df3['amino_acids'].str.len()
df3.head()
st.write(df3)

'''
### Compare Thailand and China strain
'''
# covid_DNA_th
# covid_DNA_ch

# comparison_score = pairwise2.align.globalxx(covid_DNA_th, covid_DNA_ch, one_alignment_only=True, score_only=True)
# st.write(f'{len(covid_DNA_th)} + {comparison_score} = ', len(covid_DNA_th) - comparison_score)
# st.write('SARS/COV Similarity (%):', comparison_score / len(covid_DNA_th) * 100)