# Sequence Application
This package provides a Javascript application that attempts to rebuild a DNA segment from a set of FASTA-compliant records.  The application reads the FASTA-compliant records from a specified data-file and outputs a DNA sequence based on the FASTA records.

Comments and questions contact [Glenn Inn](mailto:glenninn@yahoo.com "Glenn Inn").

## About this application ##
- The program is written in Javascript and run using NodeJS.
- The output is human readable and should be machine readable FASTA

## Getting Started ##

- Ensure that you have installed the NodeJS environment on your computer.  You can [download NodeJS](https://nodejs.org/ "NodeJS web site") from here. (The application was build using version v4.5.0)
- Download this github project.

### To Run the Sequencer ###
------------
To run the sequencer, you invoke the following command-line using Node or its package manager *npm*:
> node sequencer.js *data_file options*
> 
> npm run sequencer **--** *data_file options*

To redirect the application output to a file, use command-line redirection:

> node sequencer.js myfile.txt > results.txt

The program currently recognizes the following option(s):

-  -v:  Output in verbose mode ("node sequencer.js myfile.txt -v")


### Output File Format ###
------------

The text contained in the output should be machine-readable by FASTA-compliant applications.  By nature the output is also human-readable (it just will have a few lines that start with ';' semicolons)

The header of the file shows the following information:

    ;****************************************************
    ;* FASTA Resequencer Application
    ;* 2017, G. Inn, comments to: glenninn@yahoo.com
    ;*
    ;****************************************************
    ;Processing file: coding_challenge_data_set.txt
    ;Read [ 50 ] DNA sequences from: coding_challenge_data_set.txt
    ;There are initially < 48 > seed pairs of DNA sequences

If **verbose output** is enabled the output data includes additional summary sections:

**Articulated list** of starting DNA pairs that were tried during the chaining operation
 
 
    ;Building Total DNA from Seed Pair(0): Rosalind_1836:Rosalind_0030 fail.
    ;Building Total DNA from Seed Pair(1): Rosalind_9816:Rosalind_7890 fail.
    ;Building Total DNA from Seed Pair(2): Rosalind_5517:Rosalind_2472 fail.
    ;Building ...
    ;Building Total DNA from Seed Pair(8): Rosalind_0505:Rosalind_9944 ..success!


**Ordered List** showing the ordering of the DNA segments constituting the total DNA segment:

    ;Ordered listing of FASTA segments
    ;---------------------------------
    ;:Rosalind_0505:Rosalind_9944:Rosalind_2165:Rosalind_0372:Rosalind_0746
    ;:Rosalind_5256:Rosalind_7008:Rosalind_0093:Rosalind_2924:Rosalind_5268
    ;:Rosalind_9391 ...
    ;:Rosalind_3706:Rosalind_1836:Rosalind_0030:Rosalind_3606:Rosalind_9985

Lastly, the FASTA record showing the chained DNA

    ;Rebuilt DNA Sequence from FASTA Segments
    ;vv-------------clip here--------------vv
    >coding_challenge_data_set.txt
    TTCATACCTCGGATCTATCAATCCGAAGTACAGTAGTGCGCGATCGAAGAATCCCGCACC
    CCAGGTGTTCCTCCTACAGCGTTGTGACAACTATGTCTTGAGGTACTCGCCACACGCCCG
    TTCAGGTGGAAGAGCCGCAGGCGAGTCAACTGAGCGTACTGGATCTATGCTTCAATACTT
