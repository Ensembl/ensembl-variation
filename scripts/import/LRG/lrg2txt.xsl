<?xml version="1.0" encoding="utf-8"?>

<xsl:transform version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="text"/>

<xsl:template match="/lrg">
LRG definition: <xsl:value-of select="fixed_annotation/id"/>

FIXED ANNOTATION

Organism: <xsl:value-of select="fixed_annotation/organism"/> Taxonomy ID: <xsl:value-of select="fixed_annotation/organism/@taxon"/>

Original source: <xsl:value-of select="fixed_annotation/source/name"/> URL: <xsl:value-of select="fixed_annotation/source/url"/> Address: <xsl:value-of select="fixed_annotation/source/contact/address"/> Email: <xsl:value-of select="fixed_annotation/source/contact/email"/>

Molecule type: <xsl:value-of select="fixed_annotation/mol_type"/>

Creation date: <xsl:value-of select="fixed_annotation/creation_date"/>

Sequence: <xsl:value-of select="fixed_annotation/sequence"/>

<xsl:for-each select="fixed_annotation/transcript">

Transcript: <xsl:value-of select="@name"/> coding region: <xsl:value-of select="coding_region/cds_start"/>-<xsl:value-of select="coding_region/cds_end"/>
Translated sequence: <xsl:value-of select="coding_region/translation/sequence"/>
Exons: <xsl:for-each select="exon"> <xsl:value-of select="@lrg_number"/> (<xsl:value-of select="start"/>-<xsl:value-of select="end"/>)<xsl:if test="position()!=last()">, </xsl:if></xsl:for-each>
CDNA:  <xsl:value-of select="cdna/sequence"/>
Exons (CDNA coordinates): <xsl:for-each select="cdna/exon"> <xsl:value-of select="@lrg_number"/> (<xsl:value-of select="start"/>-<xsl:value-of select="end"/>)<xsl:if test="position()!=last()">, </xsl:if></xsl:for-each>

</xsl:for-each>

UPDATABLE ANNOTATION

Modification date: <xsl:value-of select="updatable_annotation/modification_date"/>

Modification source: <xsl:value-of select="updatable_annotation/source/name"/> URL: <xsl:value-of select="updatable_annotation/source/url"/> Address: <xsl:value-of select="updatable_annotation/source/contact/address"/> Email: <xsl:value-of select="updatable_annotation/source/contact/email"/>

Alternate exon naming (source <xsl:value-of select="updatable_annotation/alternate_exon_naming/source"/>)
Exons: <xsl:for-each select="updatable_annotation/alternate_exon_naming/exon"> current name: <xsl:value-of select="@lrg_number"/> other name: <xsl:value-of select="@other_name"/><xsl:if test="position()!=last()">, </xsl:if></xsl:for-each>

Mapping (assembly <xsl:value-of select="updatable_annotation/mapping/@assembly"/>): <xsl:for-each select="updatable_annotation/mapping/align">
chromosome <xsl:value-of select="@chromosome"/> strand <xsl:value-of select="@strand"/>  LRG start <xsl:value-of select="@lrg_start"/> LRG end <xsl:value-of select="@lrg_end"/> start <xsl:value-of select="@start"/> end <xsl:value-of select="@end"/> </xsl:for-each>

Features: 
Genes: <xsl:for-each select="updatable_annotation/features/gene"> <xsl:value-of select="@name"/> <xsl:if test="synonym"> synonym(s) <xsl:for-each select="synonym"> <xsl:value-of select="."/> </xsl:for-each> </xsl:if> start <xsl:value-of select="@start"/> end <xsl:value-of select="@end"/> Xrefs <xsl:for-each select="db_xref"><xsl:value-of select="@source"/>:<xsl:value-of select="@accession"/><xsl:if test="position()!=last()">, </xsl:if> </xsl:for-each> <xsl:if test="note"> note <xsl:value-of select="note"/>
</xsl:if> </xsl:for-each>

CDS: <xsl:for-each select="updatable_annotation/features/cds"> source <xsl:value-of select="@source"/> <xsl:value-of select="@transcript_id"/> start <xsl:value-of select="@codon_start"/> Xrefs <xsl:for-each select="db_xref"><xsl:value-of select="@source"/>:<xsl:value-of select="@accession"/><xsl:if test="position()!=last()">, </xsl:if> </xsl:for-each> </xsl:for-each>

<xsl:for-each select="updatable_annotation/features/cds/protein_product"> Protein product source=<xsl:value-of select="protein_id/@source"/> ID=<xsl:value-of select="protein_id/@accession"/> <xsl:if test="note"> note <xsl:value-of select="note"/> </xsl:if> Xrefs <xsl:for-each select="db_xref"><xsl:value-of select="@source"/>:<xsl:value-of select="@accession"/> <xsl:if test="position()!=last()">, </xsl:if> </xsl:for-each> </xsl:for-each> 

</xsl:template>

</xsl:transform>