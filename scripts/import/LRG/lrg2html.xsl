<?xml version="1.0" encoding="utf-8"?>

<xsl:transform version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:loop="http://informatik.hu-berlin.de/loop">

  <xsl:output method="html" doctype-public="-//W3C//DTD HTML 4.0 Transitional//EN"/>

  <xsl:template match="/lrg">
    <html lang="en">
      <head>
	  <title>
		<xsl:value-of select="fixed_annotation/id"/> -
		<xsl:for-each select="updatable_annotation/annotation_set/features/gene">
		  <xsl:value-of select="@symbol"/>
		  <xsl:if test="position()!=last()"> / </xsl:if>
		</xsl:for-each>
	  </title>
	  <style type="text/css">
		body {
		  font-family: "Lucida Grande", "Helvetica", "Arial", sans-serif;
	      font-size: 80%;
		}
		
		a {
		  text-decoration: none;
		  color: #48a726;
		}
		
		a:hover {
		  text-decoration: underline;
		}
		
		h1 {
		  color: #FFFFFF;
		  background: #003399;
		  border: 1px dotted black;
		  padding: 10px;
		}
		
		h2 {
		  color: #48a726;
		}
		
		h3 {
		  color: #003399;
		}
		
		h4 {
		  color: #003399;
		  font-style: italic;
		}
		
		th {
		  color: #FFFFFF;
		  background: #003399;
		  border: 1px solid #003399;
		  padding: 5px;
		}
		
		table {
		  font-size: 100%;
		  border-collapse: collapse;
		}
		
		td {
		  border: 1px solid black;
		  padding: 5px;
		}
		
		ul{
		  margin-left:10px;
		  padding-left:5px;
		}
		
		.sequence {
		  font-family: monospace;
		  border: 0px;
		  padding: 0px;
		}
		
		.coord {
		  padding: 0px;
		  padding-right: 12px;
		  text-align: right;
		  color: #666666;
		  border: 0px;
		}
		
		.showhide {
		  border: 0px;
		}
		
		.intron {
		  font-family: monospace;
		  color: #000000;
		}
		
		.exon {
		  font-family: monospace;
		  color: #003399;
		}
		
		.outphase {
		  font-family: monospace;
		  color: white;
		  background: red;
		}
		
		.exonkey {
		  color: #003399;
		}
		
		.outphasekey {
		  color: white;
		  background: red;
		}
		
		.hardbreak {
		  -moz-binding: url('./wordwrap.xml#wordwrap');
		  word-wrap: break-word; /* Internet Explorer 5.5+ */
		  border: 0px;
		  padding: 0px;
		  width: 624px;
		}
	  </style>
	  <script type="text/javascript">
		// function to show/hide layers
		function showhide(lyr) {
		  var lyrobj = document.getElementById(lyr);
		  
		  if(lyrobj.style.height == "0px") {
			lyrobj.style.height = "";
			lyrobj.style.display = "";
		  }
		  
		  else {
			lyrobj.style.height = "0px";
			lyrobj.style.display = "none";
		  }
		}
	  </script>
      </head>
      <body>
	<h1>
	  <xsl:value-of select="fixed_annotation/id"/><br/>
	  <xsl:for-each select="updatable_annotation/annotation_set/features/gene">
		- <xsl:value-of select="@symbol"/>
		<xsl:if test="long_name"> : <xsl:value-of select="long_name"/></xsl:if>
		<xsl:if test="position()!=last()"><br/></xsl:if>
	  </xsl:for-each>
	</h1>

	<h2>FIXED ANNOTATION</h2>

	<p><strong>Organism: </strong><xsl:value-of select="fixed_annotation/organism"/><br/><strong>Taxonomy ID: </strong><xsl:value-of select="fixed_annotation/organism/@taxon"/></p>

	<p>
		<strong>Original source: </strong><xsl:value-of select="fixed_annotation/source/name"/><br/>
		<strong>URL: </strong>
		<a>
			<xsl:attribute name="href">
				<xsl:if test="not(substring(fixed_annotation/source/url, 0, 4)='http')">http://</xsl:if>
				<xsl:value-of select="fixed_annotation/source/url"/>
			</xsl:attribute>
			<xsl:attribute name="target">_blank</xsl:attribute>
			<xsl:if test="not(substring(fixed_annotation/source/url, 0, 4)='http')">http://</xsl:if>
			<xsl:value-of select="fixed_annotation/source/url"/>
		</a><br/>
		<strong>Name: </strong> <xsl:value-of select="fixed_annotation/source/contact/name"/> <br/>
		<strong>Address: </strong><xsl:value-of select="fixed_annotation/source/contact/address"/><br/>
		<strong>Email: </strong><xsl:value-of select="fixed_annotation/source/contact/email"/>
	</p>

	<p><strong>Molecule type: </strong><xsl:value-of select="fixed_annotation/mol_type"/></p>

	<p><strong>Creation date: </strong><xsl:value-of select="fixed_annotation/creation_date"/></p>

	<table>
	  <tr valign="top">
		<td style="border:0px;"><h4>- Genomic sequence</h4></td>
		<td style="border:0px;"><a><xsl:attribute name="href">javascript:showhide('sequence');</xsl:attribute>show/hide</a></td>
	  </tr>
	</table>
	
	<div id="sequence" style="height: 0px; display: none;">
	  <!--<table border="0" cellpadding="0" cellspacing="0" class="sequence">
		<loop:for name="i" from="1" to="string-length(fixed_annotation/sequence)" step="100">
		  <tr>
			<td class="coord" unselectable="on"><xsl:value-of select="$i "/></td>
			<td class="sequence"><xsl:value-of select="substring(fixed_annotation/sequence,$i,100)"/></td>
		  </tr>
		</loop:for>
		<tr>
		  <td colspan="2" class="showhide">
			<a><xsl:attribute name="href">javascript:showhide('sequence');</xsl:attribute>^^ hide ^^</a>
		  </td>
		</tr>
	  </table>-->
	  
	  <table>
		<tr valign="middle">
		  <td style="border:0px;"><strong>Key: </strong></td>
		  <td style="border:0px;"><span class="exon">Exons like this</span></td>
		</tr>
		<tr>
		  <td style="border:0px;"> </td>
		  <td style="border:0px;"><span class="intron">Introns like this</span></td>
		</tr>
		<tr>
		  <td style="border:0px;"> </td>
		  <td style="border:0px;" colspan="2">Exons shown from <a><xsl:attribute name="href">#transcript_<xsl:value-of select="fixed_annotation/transcript[position() = 1]/@name"/></xsl:attribute>transcript <xsl:value-of select="fixed_annotation/transcript[position() = 1]/@name"/></a></td>
		</tr>
	  </table>
	  <table>
		<tr>
		  <td width="624" class="sequence">
			<div class="hardbreak">
			  <xsl:variable name="genseq" select="fixed_annotation/sequence"/>
			  <xsl:for-each select="fixed_annotation/transcript[position() = 1]/exon">
				<xsl:choose>
				  <xsl:when test="position()=1">
					<span class="intron">
					  <xsl:attribute name="title">Intron 1-<xsl:value-of select="lrg_coords/@start - 1"/></xsl:attribute>
					  <xsl:value-of select="substring($genseq,1,lrg_coords/@start)"/>
					</span>
				  </xsl:when>
				  <xsl:otherwise>
					<span class="intron">
					  <xsl:attribute name="title">Intron <xsl:value-of select="lrg_coords/@start"/>-<xsl:value-of select="lrg_coords/@end"/></xsl:attribute>
					  <xsl:variable name="start" select="lrg_coords/@start"/>
					  <xsl:for-each select="preceding-sibling::*/lrg_coords">
						<xsl:if test="position()=last()">
						  <xsl:value-of select="substring($genseq, @end + 1, ($start - @end) - 1)"/>
						</xsl:if>
					  </xsl:for-each>
					</span>
				  </xsl:otherwise>
				</xsl:choose>
				<span class="exon">
				  <xsl:attribute name="title">Exon <xsl:value-of select="lrg_coords/@start"/>-<xsl:value-of select="lrg_coords/@end"/></xsl:attribute>
				  <xsl:value-of select="substring($genseq,lrg_coords/@start,(lrg_coords/@end - lrg_coords/@start) + 1)"/>
				</span>
				<xsl:if test="position()=last()">
				  <xsl:if test="lrg_coords/@end &lt; string-length($genseq)">
					<span class="intron">
					  <xsl:attribute name="title">Intron <xsl:value-of select="lrg_coords/@end + 1"/>-<xsl:value-of select="string-length($genseq)"/></xsl:attribute>
					  <xsl:value-of select="substring($genseq,lrg_coords/@end + 1, string-length($genseq) - lrg_coords/@end + 1)"/>
					</span>
				  </xsl:if>
				</xsl:if>
			  </xsl:for-each>
			</div>
		  </td>
		</tr>
		<tr>
		  <td class="showhide">
			<a><xsl:attribute name="href">javascript:showhide('sequence');</xsl:attribute>^^ hide ^^</a>
		  </td>
		</tr>
	  </table>
	</div>

	<h3>Transcripts</h3>

	<xsl:for-each select="fixed_annotation/transcript">

	  <p>
		<a><xsl:attribute name="name">transcript_<xsl:value-of select="@name"/></xsl:attribute></a>
		<strong>Transcript: </strong>
		<xsl:value-of select="@name"/><br/>
		<strong>Coding region: </strong>
		<xsl:value-of select="coding_region/@start"/>-<xsl:value-of select="coding_region/@end"/><br/>
	  </p>
	  <table>
		<tr valign="top">
		  <td style="border:0px;"><h4>- cDNA sequence</h4></td>
		  <td style="border:0px;"><a><xsl:attribute name="href">javascript:showhide('cdna');</xsl:attribute>show/hide</a></td>
		</tr>
	  </table>
	  
	  <div id="cdna" style="height:0px; display: none;">
		
		<!--<table border="0" cellpadding="0" cellspacing="0" class="sequence">
		  <loop:for name="i" from="1" to="string-length(cdna/sequence)" step="100">
			<tr>
			  <td class="coord"><xsl:value-of select="$i"/></td>
			  <td class="sequence"><xsl:value-of select="substring(cdna/sequence,$i,100)"/></td>
			</tr>
		  </loop:for>
		  
		  <tr>
			<td colspan="2" class="showhide">
			  <a><xsl:attribute name="href">javascript:showhide('cdna');</xsl:attribute>^^ hide ^^</a>
			</td>
		  </tr>
		</table>-->
		
		<strong>Key: </strong><span class="exonkey">Highlighting</span> indicates alternate exons<br/><br/>
		
		<table>
		  <tr>
			<td width="624" class="sequence">
			  <div class="hardbreak">
				<xsl:variable name="seq" select="cdna/sequence"/>
				<xsl:for-each select="exon">
				  <xsl:choose>
					<xsl:when test="round(position() div 2) = (position() div 2)">
					  <span class="exon">
						<xsl:attribute name="title">Exon <xsl:value-of select="cdna_coords//@start"/>-<xsl:value-of select="cdna_coords/@end"/></xsl:attribute>
						<xsl:value-of select="substring($seq,cdna_coords/@start,(cdna_coords/@end - cdna_coords/@start) + 1)"/>
					  </span>
					</xsl:when>
					<xsl:otherwise>
					  <span class="intron">
						<xsl:attribute name="title">Exon <xsl:value-of select="cdna_coords/@start"/>-<xsl:value-of select="cdna_coords/@end"/></xsl:attribute>
						<xsl:value-of select="substring($seq,cdna_coords/@start,(cdna_coords/@end - cdna_coords/@start) + 1)"/>
					  </span>
					</xsl:otherwise>
				  </xsl:choose>
				</xsl:for-each>
			  </div>
			</td>
		  </tr>
		  <tr>
			<td class="showhide">
			  <a><xsl:attribute name="href">javascript:showhide('cdna');</xsl:attribute>^^ hide ^^</a>
			</td>
		  </tr>
		</table>
	  </div>
	  
	  
	  <xsl:variable name="transname" select="@name"/>
	  <a><xsl:attribute name="name">exons_<xsl:value-of select="$transname"/></xsl:attribute></a>
	  
	  <table>
		<tr valign="top">
		  <td style="border:0px;"><h4>- Exons</h4></td>
		  <td style="border:0px;"><a><xsl:attribute name="href">javascript:showhide('exons');</xsl:attribute>show/hide</a></td>
		</tr>
	  </table>
	  
	  <div id="exons" style="height:0px; display: none;">
		<xsl:variable name="transname" select="@name"/>
		<table>
		  <tr><th colspan="2">LRG</th><th colspan="2">cDNA</th><th colspan="2">Peptide</th><th colspan="100" style="background:#2266BB">Labels</th></tr>
		  <tr>
			<th>Start</th><th>End</th><th>Start</th><th>End</th><th>Start (phase)</th><th>End (phase)</th>
			<xsl:for-each select="/*/updatable_annotation/annotation_set">
			  <xsl:variable name="setnum" select="position()"/>
			  <xsl:for-each select="other_exon_naming/source">
				<xsl:if test="transcript[@name=$transname]"><th style="background:#2266BB"><span><xsl:attribute name="title"><xsl:value-of select="@description"/></xsl:attribute><a><xsl:attribute name="href">#source_<xsl:value-of select="$setnum"/>_<xsl:value-of select="position()"/></xsl:attribute>Source <xsl:value-of select="position()"/></a></span></th></xsl:if>
			  </xsl:for-each>
			</xsl:for-each>
		  </tr>
		  <xsl:for-each select="exon">
			<xsl:variable name="start" select="lrg_coords/@start"/>
			<tr align="right">
			  <td><xsl:value-of select="lrg_coords/@start"/></td>
			  <td><xsl:value-of select="lrg_coords/@end"/></td>
			  <td><xsl:value-of select="cdna_coords/@start"/></td>
			  <td><xsl:value-of select="cdna_coords/@end"/></td>
			  <td><xsl:value-of select="peptide_coords/@start"/> (<xsl:value-of select="peptide_coords/@start_phase"/>)</td>
			  <td><xsl:value-of select="peptide_coords/@end"/> (<xsl:value-of select="peptide_coords/@end_phase"/>)</td>
			  <xsl:for-each select="/*/updatable_annotation/*/other_exon_naming/source">
				<td>
				  <xsl:choose>
					<xsl:when test="transcript[@name=$transname]/exon/lrg_coords[@start=$start]">
					  <xsl:value-of select="transcript[@name=$transname]/exon/lrg_coords[@start=$start]/../label"/>
					</xsl:when>
					<xsl:otherwise>-</xsl:otherwise>
				  </xsl:choose>
				</td>
			  </xsl:for-each>
			</tr>
		  </xsl:for-each>
		</table>
		<p><a><xsl:attribute name="href">javascript:showhide('exons');</xsl:attribute>^^ hide ^^</a></p>
	  </div>

	  <table>
		<tr valign="top">
		  <td style="border:0px;"><h4>- Translated sequence</h4></td>
		  <td style="border:0px;"><a><xsl:attribute name="href">javascript:showhide('translated');</xsl:attribute>show/hide</a></td>
		</tr>
	  </table>
	  <div id="translated" style="height:0px; display: none;">
		<!--<table border="0" cellpadding="0" cellspacing="0" class="sequence">
		  <loop:for name="i" from="1" to="string-length(coding_region/translation/sequence)" step="100">
			<tr>
			  <td class="coord"><xsl:value-of select="$i "/></td>
			  <td class="sequence"><xsl:value-of select="substring(coding_region/translation/sequence,$i,100)"/></td>
			</tr>
		  </loop:for>
		  <tr>
			<td colspan="2" class="showhide">
			  <a><xsl:attribute name="href">javascript:showhide('translated');</xsl:attribute>^^ hide ^^</a>
			</td>
		  </tr>
		</table>-->
		
		<table>
		  <tr>
			<td style="border:0px; padding:0px;"><strong>Key: </strong></td>
			<td style="border:0px;"><span class="exonkey">Highlighting</span> indicates alternate exons</td>
		  </tr>
		  <tr>
			<td style="border:0px;"> </td>
			<td style="border:0px;"><span class="outphasekey">Shading</span> indicates exon boundary is within the codon for this amino acid</td>
		  </tr>
		</table>
		<br/>
		<table>
		  <tr>
			<td width="624" class="sequence">
			  <div class="hardbreak">
				<xsl:variable name="trans_seq" select="coding_region/translation/sequence"/>
				<xsl:for-each select="exon">
				  <xsl:choose>
					<xsl:when test="round(position() div 2) = (position() div 2)">
					  <span class="exon">
						<xsl:attribute name="title">Exon <xsl:value-of select="peptide_coords//@start"/>(<xsl:value-of select="peptide_coords/@start_phase"/>)-<xsl:value-of select="peptide_coords/@end"/>(<xsl:value-of select="peptide_coords/@end_phase"/>)</xsl:attribute>
						<xsl:choose>
						  <xsl:when test="peptide_coords/@start_phase=0">
							<xsl:choose>
							  <xsl:when test="peptide_coords/@end_phase=2">
								<xsl:value-of select="substring($trans_seq,peptide_coords/@start,(peptide_coords/@end - peptide_coords/@start) + 1)"/>
							  </xsl:when>
							  <xsl:otherwise>
								<xsl:value-of select="substring($trans_seq,peptide_coords/@start,(peptide_coords/@end - peptide_coords/@start))"/>
							  </xsl:otherwise>
							</xsl:choose>
						  </xsl:when>
						  <xsl:otherwise>
							<xsl:choose>
							  <xsl:when test="peptide_coords/@end_phase=2">
								<xsl:value-of select="substring($trans_seq,peptide_coords/@start + 1,(peptide_coords/@end - peptide_coords/@start))"/>
							  </xsl:when>
							  <xsl:otherwise>
								<xsl:value-of select="substring($trans_seq,peptide_coords/@start + 1,(peptide_coords/@end - peptide_coords/@start) - 1)"/>
							  </xsl:otherwise>
							</xsl:choose>
						  </xsl:otherwise>
						</xsl:choose>
					  </span>
					</xsl:when>
					<xsl:otherwise>
					  <span class="intron">
						<xsl:attribute name="title">Exon <xsl:value-of select="peptide_coords//@start"/>(<xsl:value-of select="peptide_coords/@start_phase"/>)-<xsl:value-of select="peptide_coords/@end"/>(<xsl:value-of select="peptide_coords/@end_phase"/>)</xsl:attribute>
						<xsl:choose>
						  <xsl:when test="peptide_coords/@start_phase=0">
							<xsl:choose>
							  <xsl:when test="peptide_coords/@end_phase=2">
								<xsl:value-of select="substring($trans_seq,peptide_coords/@start,(peptide_coords/@end - peptide_coords/@start) + 1)"/>
							  </xsl:when>
							  <xsl:otherwise>
								<xsl:value-of select="substring($trans_seq,peptide_coords/@start,(peptide_coords/@end - peptide_coords/@start))"/>
							  </xsl:otherwise>
							</xsl:choose>
						  </xsl:when>
						  <xsl:otherwise>
							<xsl:choose>
							  <xsl:when test="peptide_coords/@end_phase=2">
								<xsl:value-of select="substring($trans_seq,peptide_coords/@start + 1,(peptide_coords/@end - peptide_coords/@start))"/>
							  </xsl:when>
							  <xsl:otherwise>
								<xsl:value-of select="substring($trans_seq,peptide_coords/@start + 1,(peptide_coords/@end - peptide_coords/@start) - 1)"/>
							  </xsl:otherwise>
							</xsl:choose>
						  </xsl:otherwise>
						</xsl:choose>
					  </span>
					</xsl:otherwise>
				  </xsl:choose>
				  
				  <xsl:if test="peptide_coords/@end_phase!=2">
					<span class="outphase">
					<xsl:attribute name="title">Exon/intron boundary at <xsl:value-of select="peptide_coords/@end"/> phase <xsl:value-of select="peptide_coords/@end_phase + 1"/></xsl:attribute>
					  <xsl:value-of select="substring($trans_seq,peptide_coords/@end,1)"/>
					</span>
				  </xsl:if>  
				</xsl:for-each>
			  </div>
			</td>
		  </tr>
		  <tr>
			<td class="showhide">
			  <a><xsl:attribute name="href">javascript:showhide('translated');</xsl:attribute>^^ hide ^^</a>
			</td>
		  </tr>
		</table>
	  </div>

	</xsl:for-each>

	<hr/>

	<h2>UPDATABLE ANNOTATION</h2>

	<xsl:for-each select="updatable_annotation/annotation_set">
	  
	  <xsl:variable name="setnum" select="position()"/>

	  <p><strong>Modification date: </strong><xsl:value-of select="modification_date"/></p>
  
	  <p>
		<strong>Source: </strong><xsl:value-of select="source/name"/><br/>
		<xsl:for-each select="url">
		  <strong>URL: </strong>
		  <a>
			<xsl:attribute name="href">
				<xsl:if test="not(substring(source/url, 0, 4)='http')">http://</xsl:if>
				<xsl:value-of select="source/url"/>
			</xsl:attribute>
			<xsl:attribute name="target">_blank</xsl:attribute>
			<xsl:if test="not(substring(source/url, 0, 4)='http')">http://</xsl:if>
			<xsl:value-of select="source/url"/>
		  </a><br/>
		</xsl:for-each>
		<xsl:for-each select="source/contact">
		  <xsl:if test="name">
			<strong>Name: </strong> <xsl:value-of select="name"/> <br/>
		  </xsl:if>
		  <xsl:if test="address">
			<strong>Address: </strong><xsl:value-of select="address"/><br/>
		  </xsl:if>
		  <xsl:if test="email">
			<strong>Email: </strong><xsl:value-of select="email"/>
		  </xsl:if>
		</xsl:for-each>
	  </p>
  
	  <xsl:if test="other_exon_naming/*">
		<h3>Alternate exon naming</h3>
		<ul>
		  <xsl:for-each select="other_exon_naming/source">
			<li>
			  <a><xsl:attribute name="name">source_<xsl:value-of select="$setnum"/>_<xsl:value-of select="position()"/></xsl:attribute></a>
			  <strong>Source: </strong>
			  <xsl:choose>
				<xsl:when test="contains(@description,'www') or contains(@description,'http') or contains(@description,'/') or contains(@description,'.com')">
				  <a>
					<xsl:attribute name="href">
					  <xsl:if test="not(contains(@description,'http'))">http://</xsl:if>
					  <xsl:value-of select="@description"/>
					</xsl:attribute>
					<xsl:attribute name="target">_blank</xsl:attribute>
					<xsl:if test="not(contains(@description,'http'))">http://</xsl:if>
					<xsl:value-of select="@description"/>
				  </a>
				</xsl:when>
				<xsl:otherwise><xsl:value-of select="@description"/></xsl:otherwise>
			  </xsl:choose>
			
			  <xsl:for-each select="transcript">
				
				<xsl:variable name="transname" select="@name"/>
				
				<xsl:choose>
				  <xsl:when test="/*/fixed_annotation/transcript[@name=$transname]">
					<p>Exon labels for transcript <xsl:value-of select="$transname"/> listed <a><xsl:attribute name="href">#exons_<xsl:value-of select="$transname"/></xsl:attribute>above</a></p>
				  </xsl:when>
				  <xsl:otherwise>
					<p><strong>Transcript: </strong><xsl:value-of select="@name"/></p>
					<table>
					  <tr><th>Exon start</th><th>Exon end</th><th>Label</th></tr>
					  <xsl:for-each select="exon">
						<tr>
						  <td><xsl:value-of select="lrg_coords/@start"/></td>
						  <td><xsl:value-of select="lrg_coords/@end"/></td>
						  <td><xsl:value-of select="label"/></td>
						</tr>
					  </xsl:for-each>
					</table>
				  </xsl:otherwise>
				</xsl:choose>
			  </xsl:for-each>
			</li>
		  </xsl:for-each>
		</ul>
	  </xsl:if>
  
	  <xsl:for-each select="mapping">
		<h3>Mapping (assembly <xsl:value-of select="@assembly"/>)</h3>
		<p>
		  <strong>Region covered: </strong>
		  <xsl:value-of select="@chr_name"/>:<xsl:value-of select="@chr_start"/>-<xsl:value-of select="@chr_end"/>
		</p>
		
		<table>
		  <tr><th>Strand</th><th>LRG start</th><th>LRG end</th><th>Start</th><th>End</th><th>Differences</th></tr>
		  <xsl:for-each select="mapping_span">
			<tr>
			  <td><xsl:value-of select="@strand"/></td>
			  <td><xsl:value-of select="@lrg_start"/></td>
			  <td><xsl:value-of select="@lrg_end"/></td>
			  <td><xsl:value-of select="@start"/></td>
			  <td><xsl:value-of select="@end"/></td>
			  <td>
				<xsl:for-each select="diff">
				  <strong><xsl:value-of select="@type"/>: </strong>
				  (Ref:<xsl:value-of select="@start"/><xsl:if test="@start != @end">-<xsl:value-of select="@end"/></xsl:if>)
				  <xsl:choose>
					<xsl:when test="@genomic_sequence"><xsl:value-of select="@genomic_sequence"/></xsl:when>
					<xsl:otherwise>-</xsl:otherwise>
				  </xsl:choose>
				  ->
				  <xsl:choose>
					<xsl:when test="@lrg_sequence"><xsl:value-of select="@lrg_sequence"/></xsl:when>
					<xsl:otherwise>-</xsl:otherwise>
				  </xsl:choose>
				  (LRG:<xsl:value-of select="@lrg_start"/><xsl:if test="@lrg_start != @lrg_end">-<xsl:value-of select="@lrg_end"/></xsl:if>)
				  <br/>
				</xsl:for-each>
			  </td>
			</tr>
		  </xsl:for-each>
		</table>
	  </xsl:for-each>
  
	  <xsl:if test="alternate_amino_acid_numbering/*">
		<h3>Amino acid mapping</h3>
		<ul>
		  <xsl:for-each select="alternate_amino_acid_numbering/source">
			<li>
			  <strong>Source: </strong>
			  <xsl:choose>
				<xsl:when test="contains(@description,'www') or contains(@description,'http') or contains(@description,'/') or contains(@description,'.com')">
				  <a>
					<xsl:attribute name="href">
					  <xsl:if test="not(contains(@description,'http'))">http://</xsl:if>
					  <xsl:value-of select="@description"/>
					</xsl:attribute>
					<xsl:attribute name="target">_blank</xsl:attribute>
					<xsl:if test="not(contains(@description,'http'))">http://</xsl:if>
					<xsl:value-of select="@description"/>
				  </a>
				</xsl:when>
				<xsl:otherwise><xsl:value-of select="@description"/></xsl:otherwise>
			  </xsl:choose>
			  <xsl:for-each select="transcript">
				<p><strong>Transcript: </strong><xsl:value-of select="@name"/></p>
				<table>
				  <tr><th>LRG start</th><th>LRG end</th><th>Start</th><th>End</th></tr>
				  <xsl:for-each select="align">
					<tr>
					  <td><xsl:value-of select="@lrg_start"/></td>
					  <td><xsl:value-of select="@lrg_end"/></td>
					  <td><xsl:value-of select="@start"/></td>
					  <td><xsl:value-of select="@end"/></td>
					</tr>
				  </xsl:for-each>
				</table>
			  </xsl:for-each>
			</li>
		  </xsl:for-each>
		</ul>
	  </xsl:if>
	  
	  
	  <xsl:if test="features/*">
		  <h3>Features</h3>
  
		  <xsl:if test="features/gene/*">
			  <h4>Genes</h4> 
		  
			  <table>
		  
			  <xsl:for-each select="features/gene">
				<tr><th><xsl:value-of select="@symbol"/></th><th> </th></tr>
				<tr valign="top">
				  <td style="border:0px;">
					<p>
					  <xsl:for-each select="long_name">
						<xsl:value-of select="."/><br/>
					  </xsl:for-each>
					  <xsl:if test="partial">
						Feature partially overlaps LRG (<xsl:for-each select="partial"><xsl:value-of select="."/><xsl:if test="position()!=last()">, </xsl:if></xsl:for-each>)<br/>
					  </xsl:if>
					</p>
					
					<p>
					  <strong>Synonyms: </strong>
					  <xsl:choose>
						<xsl:when test="synonym">
						  <xsl:for-each select="synonym">
							<xsl:value-of select="."/>
							<xsl:if test="position()!=last()">, </xsl:if>
						  </xsl:for-each>
						</xsl:when>
						<xsl:otherwise>-</xsl:otherwise>
					  </xsl:choose>
					  <br/>
					  
					  <strong>LRG coords: </strong>
					  <xsl:value-of select="@start"/>-<xsl:value-of select="@end"/>
					  <br/>
					  
					  <strong>External identifiers:</strong>
					  <xsl:choose>
						<xsl:when test="db_xref">
							<xsl:for-each select="db_xref">
							  <br/>-
								<strong><xsl:value-of select="@source"/>: </strong>
							  
								<xsl:choose>
								  <xsl:when test="@source='RefSeq'">
									<a>
									  <xsl:attribute name="href">
										<xsl:choose>
										  <xsl:when test="contains(@accession,'NP')">http://www.ncbi.nlm.nih.gov/protein/<xsl:value-of select="@accession"/></xsl:when>
										  <xsl:otherwise>http://www.ncbi.nlm.nih.gov/nuccore/<xsl:value-of select="@accession"/></xsl:otherwise>
										</xsl:choose>
									  </xsl:attribute>
									  <xsl:attribute name="target">_blank</xsl:attribute>
									  <xsl:value-of select="@accession"/>
									</a>
								  </xsl:when>
								  <xsl:when test="@source='Ensembl'">
									<a>
									  <xsl:attribute name="href">
										<xsl:choose>
										  <xsl:when test="contains(@accession,'ENST')">http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=<xsl:value-of select="@accession"/></xsl:when>
										  <xsl:when test="contains(@accession,'ENSG')">http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=<xsl:value-of select="@accession"/></xsl:when>
										  <xsl:when test="contains(@accession,'ENSP')">http://www.ensembl.org/Homo_sapiens/Transcript/Idhistory/Protein?protein=<xsl:value-of select="@accession"/></xsl:when>
										  <xsl:otherwise>http://www.ensembl.org/Homo_sapiens/<xsl:value-of select="@accession"/></xsl:otherwise>
										</xsl:choose>
									  </xsl:attribute>
									  <xsl:attribute name="target">_blank</xsl:attribute>
									  <xsl:value-of select="@accession"/>
									</a>
								  </xsl:when>
								  <xsl:when test="@source='UniProtKB'">
									<a>
									  <xsl:attribute name="href">http://www.uniprot.org/uniprot/<xsl:value-of select="@accession"/></xsl:attribute>
									  <xsl:attribute name="target">_blank</xsl:attribute>
									  <xsl:value-of select="@accession"/>
									</a>
								  </xsl:when>
								  <xsl:when test="@source='CCDS'">
									<a>
									  <xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi?REQUEST=ALLFIELDS&amp;DATA=<xsl:value-of select="@accession"/></xsl:attribute>
									  <xsl:attribute name="target">_blank</xsl:attribute>
									  <xsl:value-of select="@accession"/>
									</a>
								  </xsl:when>
								  <xsl:when test="@source='GeneID'">
									<a>
									  <xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&amp;cmd=Retrieve&amp;dopt=Graphics&amp;list_uids=<xsl:value-of select="@accession"/></xsl:attribute>
									  <xsl:attribute name="target">_blank</xsl:attribute>
									  <xsl:value-of select="@accession"/>
									</a>
								  </xsl:when>
								  <xsl:when test="@source='HGNC'">
									<a>
									  <xsl:attribute name="href">http://www.genenames.org/data/hgnc_data.php?hgnc_id=<xsl:value-of select="@accession"/></xsl:attribute>
									  <xsl:attribute name="target">_blank</xsl:attribute>
									  <xsl:value-of select="@accession"/>
									</a>
								  </xsl:when>
								  <xsl:when test="@source='MIM'">
									<a>
									  <xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=<xsl:value-of select="@accession"/></xsl:attribute>
									  <xsl:attribute name="target">_blank</xsl:attribute>
									  <xsl:value-of select="@accession"/>
									</a>
								  </xsl:when>
								  
								</xsl:choose>
							</xsl:for-each>
						</xsl:when>
						<xsl:otherwise>-</xsl:otherwise>
					  </xsl:choose>
					  <br/>
					  
					  <xsl:if test="comment">
						<strong>Comments: </strong>
						<xsl:for-each select="comment">
						  <xsl:value-of select="."/>
						  <xsl:if test="position()!=last()"><br/></xsl:if>
						</xsl:for-each>
					  </xsl:if>
					</p>
				  </td>
				  
				  <!--Transcripts-->
				  <td style="padding:0px; border:0px;">
					<table width="100%">
					  <tr><th>Transcript ID</th><th>Source</th><th>Start</th><th>End</th><th>External identifiers</th><th>Other</th></tr>
					  <xsl:for-each select="transcript">
						<tr valign="top">
						  <td><xsl:value-of select="@transcript_id"/></td>
						  <td><xsl:value-of select="@source"/></td>
						  <td><xsl:value-of select="@start"/></td>
						  <td><xsl:value-of select="@end"/></td>
						  <td>
							<xsl:choose>
							  <xsl:when test="db_xref">
								
								<xsl:for-each select="db_xref">
								  <strong><xsl:value-of select="@source"/>: </strong>
								  
								  <xsl:choose>
									<xsl:when test="@source='RefSeq'">
									  <a>
										<xsl:attribute name="href">
										  <xsl:choose>
											<xsl:when test="contains(@accession,'NP')">http://www.ncbi.nlm.nih.gov/protein/<xsl:value-of select="@accession"/></xsl:when>
											<xsl:otherwise>http://www.ncbi.nlm.nih.gov/nuccore/<xsl:value-of select="@accession"/></xsl:otherwise>
										  </xsl:choose>
										</xsl:attribute>
										<xsl:attribute name="target">_blank</xsl:attribute>
										<xsl:value-of select="@accession"/>
									  </a>
									</xsl:when>
									<xsl:when test="@source='Ensembl'">
									  <a>
										<xsl:attribute name="href">
										  <xsl:choose>
											<xsl:when test="contains(@accession,'ENST')">http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=<xsl:value-of select="@accession"/></xsl:when>
											<xsl:when test="contains(@accession,'ENSG')">http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=<xsl:value-of select="@accession"/></xsl:when>
											<xsl:when test="contains(@accession,'ENSP')">http://www.ensembl.org/Homo_sapiens/Transcript/Idhistory/Protein?protein=<xsl:value-of select="@accession"/></xsl:when>
											<xsl:otherwise>http://www.ensembl.org/Homo_sapiens/<xsl:value-of select="@accession"/></xsl:otherwise>
										  </xsl:choose>
										</xsl:attribute>
										<xsl:attribute name="target">_blank</xsl:attribute>
										<xsl:value-of select="@accession"/>
									  </a>
									</xsl:when>
									<xsl:when test="@source='UniProtKB'">
									  <a>
										<xsl:attribute name="href">http://www.uniprot.org/uniprot/<xsl:value-of select="@accession"/></xsl:attribute>
										<xsl:attribute name="target">_blank</xsl:attribute>
										<xsl:value-of select="@accession"/>
									  </a>
									</xsl:when>
									<xsl:when test="@source='CCDS'">
									  <a>
										<xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi?REQUEST=ALLFIELDS&amp;DATA=<xsl:value-of select="@accession"/></xsl:attribute>
										<xsl:attribute name="target">_blank</xsl:attribute>
										<xsl:value-of select="@accession"/>
									  </a>
									</xsl:when>
									<xsl:when test="@source='GeneID'">
									  <a>
										<xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&amp;cmd=Retrieve&amp;dopt=Graphics&amp;list_uids=<xsl:value-of select="@accession"/></xsl:attribute>
										<xsl:attribute name="target">_blank</xsl:attribute>
										<xsl:value-of select="@accession"/>
									  </a>
									</xsl:when>
									<xsl:when test="@source='HGNC'">
									  <a>
										<xsl:attribute name="href">http://www.genenames.org/data/hgnc_data.php?hgnc_id=<xsl:value-of select="@accession"/></xsl:attribute>
										<xsl:attribute name="target">_blank</xsl:attribute>
										<xsl:value-of select="@accession"/>
									  </a>
									</xsl:when>
									<xsl:when test="@source='MIM'">
									  <a>
										<xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=<xsl:value-of select="@accession"/></xsl:attribute>
										<xsl:attribute name="target">_blank</xsl:attribute>
										<xsl:value-of select="@accession"/>
									  </a>
									</xsl:when>
									
								  </xsl:choose>
								  
								  <xsl:if test="position()!=last()"><br/></xsl:if>
								</xsl:for-each>
								
							  </xsl:when>
							  <xsl:otherwise>-</xsl:otherwise>
							</xsl:choose>
						  </td>
						  <td>
							<xsl:if test="partial">
							  Feature partially overlaps LRG (<xsl:for-each select="partial"><xsl:value-of select="."/><xsl:if test="position()!=last()">, </xsl:if></xsl:for-each>)<br/>
							</xsl:if>
							<xsl:if test="long_name">
							  <strong>Name: </strong><xsl:value-of select="long_name"/><br/>
							</xsl:if>
							<xsl:for-each select="comment">
							  <strong>Comment: </strong><xsl:value-of select="."/><br/>
							</xsl:for-each>
						  </td>
						</tr>
					  </xsl:for-each>
					  
					  <tr><th>Protein ID</th><th>Source</th><th>CDS start</th><th>CDS end</th><th>External identifiers</th><th>Other</th></tr>
					  <xsl:for-each select="transcript">
						<xsl:for-each select="protein_product">
						  <tr valign="top">
							<td><xsl:value-of select="@accession"/></td>
							<td><xsl:value-of select="@source"/></td>
							<td><xsl:value-of select="@cds_start"/></td>
							<td><xsl:value-of select="@cds_end"/></td>
							<td>
							  <xsl:choose>
								<xsl:when test="db_xref">
								  
								  <xsl:for-each select="db_xref">
									<strong><xsl:value-of select="@source"/>: </strong>
									
									<xsl:choose>
									  <xsl:when test="@source='RefSeq'">
										<a>
										  <xsl:attribute name="href">
											<xsl:choose>
											  <xsl:when test="contains(@accession,'NP')">http://www.ncbi.nlm.nih.gov/protein/<xsl:value-of select="@accession"/></xsl:when>
											  <xsl:otherwise>http://www.ncbi.nlm.nih.gov/nuccore/<xsl:value-of select="@accession"/></xsl:otherwise>
											</xsl:choose>
										  </xsl:attribute>
										  <xsl:attribute name="target">_blank</xsl:attribute>
										  <xsl:value-of select="@accession"/>
										</a>
									  </xsl:when>
									  <xsl:when test="@source='Ensembl'">
										<a>
										  <xsl:attribute name="href">
											<xsl:choose>
											  <xsl:when test="contains(@accession,'ENST')">http://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=<xsl:value-of select="@accession"/></xsl:when>
											  <xsl:when test="contains(@accession,'ENSG')">http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=<xsl:value-of select="@accession"/></xsl:when>
											  <xsl:when test="contains(@accession,'ENSP')">http://www.ensembl.org/Homo_sapiens/Transcript/Idhistory/Protein?protein=<xsl:value-of select="@accession"/></xsl:when>
											  <xsl:otherwise>http://www.ensembl.org/Homo_sapiens/<xsl:value-of select="@accession"/></xsl:otherwise>
											</xsl:choose>
										  </xsl:attribute>
										  <xsl:attribute name="target">_blank</xsl:attribute>
										  <xsl:value-of select="@accession"/>
										</a>
									  </xsl:when>
									  <xsl:when test="@source='UniProtKB'">
										<a>
										  <xsl:attribute name="href">http://www.uniprot.org/uniprot/<xsl:value-of select="@accession"/></xsl:attribute>
										  <xsl:attribute name="target">_blank</xsl:attribute>
										  <xsl:value-of select="@accession"/>
										</a>
									  </xsl:when>
									  <xsl:when test="@source='CCDS'">
										<a>
										  <xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/projects/CCDS/CcdsBrowse.cgi?REQUEST=ALLFIELDS&amp;DATA=<xsl:value-of select="@accession"/></xsl:attribute>
										  <xsl:attribute name="target">_blank</xsl:attribute>
										  <xsl:value-of select="@accession"/>
										</a>
									  </xsl:when>
									  <xsl:when test="@source='GeneID'">
										<a>
										  <xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&amp;cmd=Retrieve&amp;dopt=Graphics&amp;list_uids=<xsl:value-of select="@accession"/></xsl:attribute>
										  <xsl:attribute name="target">_blank</xsl:attribute>
										  <xsl:value-of select="@accession"/>
										</a>
									  </xsl:when>
									  <xsl:when test="@source='HGNC'">
										<a>
										  <xsl:attribute name="href">http://www.genenames.org/data/hgnc_data.php?hgnc_id=<xsl:value-of select="@accession"/></xsl:attribute>
										  <xsl:attribute name="target">_blank</xsl:attribute>
										  <xsl:value-of select="@accession"/>
										</a>
									  </xsl:when>
									  <xsl:when test="@source='MIM'">
										<a>
										  <xsl:attribute name="href">http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=<xsl:value-of select="@accession"/></xsl:attribute>
										  <xsl:attribute name="target">_blank</xsl:attribute>
										  <xsl:value-of select="@accession"/>
										</a>
									  </xsl:when>
									  
									</xsl:choose>
									
									<xsl:if test="position()!=last()"><br/></xsl:if>
								  </xsl:for-each>
								  
								</xsl:when>
								<xsl:otherwise>-</xsl:otherwise>
							  </xsl:choose>
							</td>
							<td>
							  <xsl:if test="partial">
								Feature partially overlaps LRG (<xsl:for-each select="partial"><xsl:value-of select="."/><xsl:if test="position()!=last()">, </xsl:if></xsl:for-each>)<br/>
							  </xsl:if>
							  <xsl:if test="long_name">
								<strong>Name: </strong><xsl:value-of select="long_name"/><br/>
							  </xsl:if>
							  <xsl:for-each select="comment">
								<strong>Comment: </strong><xsl:value-of select="."/><br/>
							  </xsl:for-each>
							</td>
						  </tr>
						</xsl:for-each> 
					  </xsl:for-each>
					  
					  <tr><td colspan="6" style="border:0px;"> </td></tr>
					</table>
				  </td>
				</tr>
			  </xsl:for-each>
			</table>
		  </xsl:if>
	  </xsl:if>
	</xsl:for-each>
    </body>
    </html>
  </xsl:template>

</xsl:transform>