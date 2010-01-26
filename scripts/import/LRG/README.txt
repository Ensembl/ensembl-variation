# LRG.rnc version 1.5 Release notes
###################################

UPDATABLE_ANNOTATION section:

- most_recent attribute, indicating if the assembly is the most recent available,
  added to mapping element

- lrg_gene_name element added to hold the HGNC gene identifier (optional)

- fixed_id attribute added to transcript element to link the updatable annotation
  to a transcript in the fixed annotation section (optional)

- alternative_amino_acid_numbering element is now optional

- other_exon_naming element is now optional

- features element is now optional

# LRG.rnc version 1.4 Release notes
###################################

FIXED_ANNOTATION section:

- exons being followed by an intron now has an element denoting intron phase

# LRG.rnc version 1.3 Release notes
###################################

FIXED_ANNOTATION section:

- minor modification of intron phase representation

# LRG.rnc version 1.2 Release notes
###################################

FIXED_ANNOTATION section:

- multiple contacts now permitted in source section.

- transcript elements now have start and end attributes.

- multiple transcript elements now permitted.

- cdna element containing sequence moved to be first element in transcript.

- coding_region element now has start, end and optional codon_start and
  codon_selenocysteine attributes.

- exon elements now follow coding_region element, rather than being contained
  within them.
  
- exons are no longer explicitly numbered in the fixed layer. Legacy
  naming/numbering can be added in the updatable layer using other_exon_naming.

- exons now have three elements within them; lrg_coords, cdna_coords and
  peptide_coords. These represent exon coordinates in three systems, and each
  have start and end attributes. The peptide_coords elements have start_phase
  and end_phase attributes with values of '0', '1' or '2' to represent the phase
  of the exon/intron boundary.



UPDATABLE_ANNOTATION section:

- now contains one or more annotation_set elements to allow for sets of
  annotation from multiple sources (e.g. NCBI, EBI).

- amino_acid_mapping element renamed alternate_amino_acid_numbering.

- alternate_amino_acid_numbering and other_exon_naming elements now have a
  transcript element with name attribute within them, to allow for
  transcript-specific annotation.

- mapping element has been completely reworked. It has chr_name, chr_start and
  chr_end attributes to represent the total region covered by the mapped LRG
  sequence. It then contains one or more mapping_span elements that represent a
  contiguous block of the LRG sequence aligned to the genome. Each mapping_span
  element may contain one or more diff elements that describe short sequence
  differences between the LRG sequence and the reference genome - these have a
  type attribute that can be 'mismatch', 'lrg_ins' or 'genomic_ins'. This new
  element configuration allows for a complete description of how the LRG
  sequence relates to and differs from the reference seqeunce.

- cds elements in features have been renamed transcript, and now occur within
  their respective gene elements. The transcript elements have start and end
  attributes to bring them in line with their equivalents in the fixed
  annotation layer.

- protein_product elements have cds_start and cds_end attributes to represent
  the coordinates of the coding sequence.

- the aforementioned gene, transcript and protein_product elements may now
  contain a partial element that indicates if the feature in question partially
  overlaps the 5' or 3' end of the LRG sequence (i.e. the feature is not fully
  contained within the LRG).
