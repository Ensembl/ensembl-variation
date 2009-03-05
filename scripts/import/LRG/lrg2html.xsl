<?xml version="1.0" encoding="utf-8"?>

<xsl:transform version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

  <xsl:output method="html"/>

  <xsl:template match="/lrg">
    <html>
      <head>
	<title><xsl:value-of select="fixed_annotation/id"/></title>
      </head>
      <body>
	<h1><xsl:value-of select="fixed_annotation/id"/></h1>

	<h2>Fixed Annotation</h2>

	<p><strong>Organism: </strong><xsl:value-of select="fixed_annotation/organism"/><br/><strong>Taxonomy ID: </strong><xsl:value-of select="fixed_annotation/organism/@taxon"/></p>

	<p>
		<strong>Original source: </strong><xsl:value-of select="fixed_annotation/source/name"/><br/>
		<strong>URL: </strong>
		<a>
			<xsl:attribute name="href">
				<xsl:value-of select="fixed_annotation/source/url"/>
			</xsl:attribute>
			<xsl:value-of select="fixed_annotation/source/url"/>
		</a><br/>
		<strong>Name: </strong> <xsl:value-of select="fixed_annotation/source/contact/name"/> <br/>
		<strong>Address: </strong><xsl:value-of select="fixed_annotation/source/contact/address"/><br/>
		<strong>Email: </strong><xsl:value-of select="fixed_annotation/source/contact/email"/>
	</p>

	<p><strong>Molecule type: </strong><xsl:value-of select="fixed_annotation/mol_type"/></p>

	<p><strong>Creation date: </strong><xsl:value-of select="fixed_annotation/creation_date"/></p>

	<p><strong>Sequence:</strong> <code><xsl:value-of select="fixed_annotation/sequence"/></code></p>

	<h3>Transcripts</h3>

	<xsl:for-each select="fixed_annotation/transcript">

	  <p><strong>Transcript: </strong><xsl:value-of select="@name"/><br/><strong>Coding region: </strong><xsl:value-of select="coding_region/cds_start"/>-<xsl:value-of select="coding_region/cds_end"/><br/>

	  <p><strong>Translated sequence:</strong> <code><xsl:value-of select="coding_region/translation/sequence"/></code></p>

	  <strong>Exons: </strong>
	  <table border='1'>
	    <tr><th>Name</th><th>Start</th><th>End</th></tr>
	    <xsl:for-each select="exon">
	      <tr>
		<td><xsl:value-of select="@lrg_number"/></td>
		<td><xsl:value-of select="start"/></td>
		<td><xsl:value-of select="end"/></td>
	    </tr></xsl:for-each>
	  </table>

	  <p><strong>CDNA: </strong><code><xsl:value-of select="cdna/sequence"/></code></p>
	  <strong>Exons (cDNA coordinates): </strong>
	  <table border='1'>
	    <tr><th>Name</th><th>Start</th><th>End</th></tr>
	    <xsl:for-each select="cdna/exon">
	      <tr>
		<td><xsl:value-of select="@lrg_number"/></td>
		<td><xsl:value-of select="start"/></td>
		<td><xsl:value-of select="end"/></td>
	    </tr></xsl:for-each>
	  </table>
	  </p>

	</xsl:for-each>

	<hr/>

	<h2>UPDATABLE ANNOTATION</h2>

	<p><strong>Modification date:</strong> <xsl:value-of select="updatable_annotation/modification_date"/></p>

	<p>
		<strong>Original source: </strong><xsl:value-of select="updatable_annotation/source/name"/><br/>
		<strong>URL: </strong>
		<a>
			<xsl:attribute name="href">
				<xsl:value-of select="updatable_annotation/source/url"/>
			</xsl:attribute>
			<xsl:value-of select="updatable_annotation/source/url"/>
		</a><br/>
		<strong>Name: </strong> <xsl:value-of select="updatable_annotation/source/contact/name"/> <br/>
		<strong>Address: </strong><xsl:value-of select="updatable_annotation/source/contact/address"/><br/>
		<strong>Email: </strong><xsl:value-of select="updatable_annotation/source/contact/email"/>
	</p>

	<xsl:if test="updatable_annotation/other_exon_naming/*">
		<p><strong>Alternate exon naming</strong><br/>
		<xsl:for-each select="updatable_annotation/other_exon_naming/source">
		Source <xsl:value-of select="@description"/><br/>
		
		<table border='1'>
			<tr><th>LRG number</th><th>Other name</th></tr>
			<xsl:for-each select="exon">
			<tr>
			<td><xsl:value-of select="@lrg_number"/></td>
			<td><xsl:value-of select="@other_name"/></td>
			</tr>
			</xsl:for-each>
		</table>
		</xsl:for-each></p>
	</xsl:if>

	<xsl:if test="updatable_annotation/mapping/*">
		<h3>Mapping</h3> 
		(assembly <xsl:value-of select="updatable_annotation/mapping/@assembly"/>)<br/>
		<p>
		<table border='1'>
			<tr><th>Chromosome</th><th>Strand</th><th>LRG start</th><th>LRG end</th><th>Start</th><th>End</th></tr>
	
			<xsl:for-each select="updatable_annotation/mapping/align">
			<tr>
			<td><xsl:value-of select="@chromosome"/></td>
			<td><xsl:value-of select="@strand"/></td>
			<td><xsl:value-of select="@lrg_start"/></td>
			<td><xsl:value-of select="@lrg_end"/></td>
			<td><xsl:value-of select="@start"/></td>
			<td><xsl:value-of select="@end"/></td>
			</tr>
			</xsl:for-each>
		</table></p>
	</xsl:if>

	<xsl:if test="updatable_annotation/amino_acid_mapping/*">
		<h3>Amino acid mapping</h3> 
		<xsl:for-each select="updatable_annotation/amino_acid_mapping/source">
		<p>Source <xsl:value-of select="@description"/><br/>
		<table border='1'>
			<tr><th>LRG start</th><th>LRG end</th><th>Start</th><th>End</th></tr>
	
			<xsl:for-each select="align">
			<tr>
			<td><xsl:value-of select="@lrg_start"/></td>
			<td><xsl:value-of select="@lrg_end"/></td>
			<td><xsl:value-of select="@start"/></td>
			<td><xsl:value-of select="@end"/></td>
			</tr>
			</xsl:for-each>
		</table></p>
		</xsl:for-each>
	</xsl:if>
	
	
	<xsl:if test="updatable_annotation/features/*">
		<h3>Features</h3>

		<xsl:if test="updatable_annotation/features/gene/*">
			<h4>Genes</h4> 
		
			<table border="1">
			<tr><th>Name</th><th>Synonym(s)</th><th>Start</th><th>End</th><th>Xrefs</th><th>Note</th></tr>
		
			<xsl:for-each select="updatable_annotation/features/gene">
				<tr>
				<td><xsl:value-of select="@name"/></td>
				<td><xsl:for-each select="synonym"><xsl:value-of select="."/></xsl:for-each></td>
				<td><xsl:value-of select="@start"/></td>
				<td><xsl:value-of select="@end"/></td>
				<td><xsl:for-each select="db_xref"><xsl:value-of select="@source"/>:<xsl:value-of select="@accession"/><xsl:if test="position()!=last()">, </xsl:if> </xsl:for-each></td>
				<td><xsl:value-of select="note"/></td>
				</tr>
			</xsl:for-each>
			</table>
		</xsl:if>

		<xsl:if test="updatable_annotation/features/cds/*">
			<h4>CDS</h4>
			
			<table border="1">
			<tr><th>Source</th><th>Transcript ID</th><th>Codon start</th><th>Codon end</th><th>Xrefs</th></tr>
			<xsl:for-each select="updatable_annotation/features/cds">
				<tr>
				<td><xsl:value-of select="@source"/></td>
				<td><xsl:value-of select="@transcript_id"/></td>
				<td><xsl:value-of select="@codon_start"/></td>
				<td><xsl:for-each select="db_xref"><xsl:value-of select="@source"/>:<xsl:value-of select="@accession"/><xsl:if test="position()!=last()">, </xsl:if> </xsl:for-each></td>
				</tr>
			</xsl:for-each>
			</table>
		
			<xsl:if test="updatable_annotation/features/cds/protein_product/*">
				<h4>Protein product</h4>
			
				<table border="1">
				<tr><th>Source</th><th>ID</th><th>Xrefs</th><th>Note</th></tr>
				
				<xsl:for-each select="updatable_annotation/features/cds/protein_product"> 
					<tr>
					<td><xsl:value-of select="protein_id/@source"/></td>
					<td><xsl:value-of select="protein_id/@accession"/></td>
					<td><xsl:for-each select="db_xref"><xsl:value-of select="@source"/>:<xsl:value-of select="@accession"/> <xsl:if test="position()!=last()">, </xsl:if> </xsl:for-each></td>
					<td><xsl:value-of select="note"/></td> 
					</tr>
				</xsl:for-each> 
				</table>
			</xsl:if>
		</xsl:if>
	</xsl:if>

      </body>
    </html>
  </xsl:template>

</xsl:transform>