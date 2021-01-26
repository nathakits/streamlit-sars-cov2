# SARS-CoV-2 analysis with BioPython and Streamlit
## Introduction
Analyze SARS-CoV-2 genome from Thailand and China.
- Nucleotide frequency
- Transcription (DNA -> mRNA)
- Translation (mRNA -> Amino acids)
- Compare DNA with pairwise2

## Build locally
Instill Anaconda if you haven't already.

Set up and activate development environment:
```bash
conda env update --file environment.yml --name pyviz
conda activate pyviz
```

Run web app
```bash
streamlit run covid19.py
```

Test Streamlit
```bash
streamlit hello
```