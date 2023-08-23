# Consensus alignment

## Break down of logic used in the manual creation of minimal consensus sequences

### Initial Setup

Download latest versions of sequences from NCBI/Ensembl/Genbank etc.
Using MUSCLE, align all sequences into one MSA
    In UGene, just default settings
Create a tree to order sequences by similarity
    PHYLIP Neighbor Joining
    DMM F84
    Transition/Transversion ratio 2.00
Break down into chunks of alignments consisting of the lowest two levels of the tree
    Probably case specific or subjective
    For HMR, mouse and rat are usually nearest neighbors in the tree with human being slightly farther out
    Some outliers need to either be saved separately or a judgement call made if they're at least somewhat close to another sequence
    Probably should solidify an actual threshold here for grouping
    So as a practical illustration of a group:

```txt
           |----------- Human
    ...----|    |------ Mouse
           |----|
                |------ Rat
```

### For each neighbor grouping

Re-align with muscle to get as complete a consensus as possible to start with
Compare base by base until a gap is reached

#### Single, simple base variant

This is the case where there is just a position which doesn't have a majority agreement.
Match the residue from whichever strand has the fewest "wins" in this category
For ties, or if this is the first residue, just pick at random

#### Single base variant, one/minority sequences has a gap

If the remaining bases match, just use that residue for consensus
Otherwise, follow the above logic to determine the victor

#### Single base, majority sequences have a gap

Keep this base empty in consensus

#### Short (<15 bp) high variance region

..TGA---AAC..
..TCGTACATC..
..TTACTTA-C..

Factor in any ongoing bias towards any given sequence, but err on the side of just trimming this out of the consensus

#### Tiny region (2-4 bp), mild variation

AA
CT
TC

Use up and downstream factors to decide what base wins. In this example if we look at the alignment 10bp to either side:

1 AACTAGAAATAACTTTGCAAGGA
2 AACTAGAAAACTTCTAACTAAAA
3 AACTAGAAAATCCTTAACAAAAA
C AACTAGAAAA--CTTAACAAAAA

Sequence 2 and 3 have fewer mismatches to the consensus than sequence 1, so for those 2 bases use sequence 1 as a reference

Ties can be split up using a base or two from each reference sequence to fill in the consensus
