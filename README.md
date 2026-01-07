# codon-optimizer
Codon Optimizer tool which takes DNA/protein sequences and optimizes them with regards to their codon bias of the organism of choice. This improves protein expression and cloning efficiency.

Its an web application implemented in two different ways: as React/JavaScript and as Python /Streamlit project.
Each has its pros and cons.
React/JavaScript can be rendered in a browser and be instantly ready. To run it locally Node.js has to be installed and it is more difficult to embed hardcore bioinformatics pipelines.
Python/Streamlit needs a Python Interpreter and streamlit installed to run. The App has to be called app.py. But this version is easily extentable with BioPython or other bioinformatics packages available for Python.

What the Codon Optimizer does:
- reads in pubicy available condon usages tables for 3 most common species (can be extended).
- user inputs their DNA or Protein Sequence and the organism of choice
- validates the input: DNA contains only A,T,C,G and protein only 20 amino acids, Length of DNA is dividable by 3.
- for each codon it translates in an amino acid and finds the most frequent codon in the organism of choice.
- output is the new sequence with the changed (aka optimized) codons in green, some statistics and the codon usage bar chart.

What possible extensions are:
- calculate and display the Codon Adaption Index (CAI). This index measures how well a codon in the source sequence matches with the codon bias of an organism of choice. Usually the CAI should get higher after the sequence was optimized using the algorithm described above.
- for cloning efficiency, a goog distribution of nucleotides across the sequence is important. That is where the GC content comes in. It should not be too high and not be too low, around 50% is acceptable.
- detecting repeats and reducing them.

These multi-factors (codon bias, GC content and repeats) needs to be all considered when optimizing the codons. The are equally important when optimizing a sequence for better expressability and cloning efficiency.
