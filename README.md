I built a Python tool using the Biopython library that takes a DNA or RNA sequence as input and processes it through three steps.
First, it validates the sequence to make sure it only contains valid nucleotide characters. 
Then it transcribes the DNA into RNA, and finally translates the RNA into a protein sequence starting from the AUG start codon. 
To test whether my tool was actually working correctly, I used the real human insulin gene sequence from NCBI and ran it through my tool. 
The protein sequence my tool produced was then compared to the official human insulin protein entry on UniProt. 
The two sequences matched completely, which confirms that my tool is correctly performing transcription and translation and producing biologically accurate results. 
