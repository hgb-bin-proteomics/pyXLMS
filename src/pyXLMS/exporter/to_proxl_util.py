__local_schema = r"""<?xml version="1.0" encoding="UTF-8"?>
<!-- edited with XMLSpy v2019 (http://www.altova.com) by Michael Riffle (University of Washington) -->
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:vc="http://www.w3.org/2007/XMLSchema-versioning" elementFormDefault="qualified" attributeFormDefault="unqualified" version="1.0" vc:minVersion="1.1">
	<xs:element name="proxl_input">
		<xs:annotation>
			<xs:documentation>Root element for all data</xs:documentation>
		</xs:annotation>
		<!-- Element name must match string in Java class XSD_ElementAttributeNamesConstants  -->
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_program_info"/>
				<xs:element ref="linkers">
					<xs:annotation>
						<xs:documentation>The chemical linkers used in this experiment. There must be at least one.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="reported_peptides">
					<xs:annotation>
						<xs:documentation>All the reported peptides from a PSM search. A reported peptide is a string that uniquely identifies an identified peptide (or linked pair of peptides) IN THIS SEARCH. (This string is not used by ProXL for this purpose.) A unique peptide i this context is one or two peptide sequences linked together at specific positions in those sequences, and with modifications at specific sites. E.g. PEP[16.01]TIDE would uniquely ID a peptide with this sequence with the given mass shift on the given position. or KEKTIDE(3)--PEPT[16.1]IDE(2) would uniquely identify the crosslink between two peptides with these ssequences at positions 3 and 2 (respectively) and where the second peptide has the shown modification. There is no syntax expected for reported peptide, but the reported peptide string must be unique, and no two reported peptides must describe the same peptide or pair of peptides linked at the same positions with the same modifications.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="matched_proteins" minOccurs="0">
					<xs:annotation>
						<xs:documentation>Contains the protein sequences matched by peptides found in the experiment, and associated names, descriptions, and taxonomy IDs for those proteins. This will typically come from the FASTA file used to search the data. Only the proteins matched by the data should be included, not all the proteins in the FASTA file.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="static_modifications" minOccurs="0"/>
				<xs:element ref="decoy_labels" minOccurs="0">
					<xs:annotation>
						<xs:documentation>An optional list of decoy labels (strings that are prefixes on protein names that indicate the given protein is a decoy (e.g., "DECOY_") </xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="configuration_files" minOccurs="0">
					<xs:annotation>
						<xs:documentation>(Optional) The configuration files associated with the search.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
			<xs:attribute name="fasta_filename" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The filename (not including the path) of the FASTA file used to search the data.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="name" type="xs:string">
				<xs:annotation>
					<xs:documentation>An optional name to be used for this search in ProXL. If not present, a default will be used. This name may be changed within ProXL.</xs:documentation>
				</xs:annotation>
				<!-- Attribute name must match string in Java class XSD_ElementAttributeNamesConstants  -->
			</xs:attribute>
			<xs:attribute name="comment" type="xs:string">
				<xs:annotation>
					<xs:documentation>An optional comment to be associated with this search when it is imported into ProXL.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="search_program">
		<xs:annotation>
			<xs:documentation>A search program, such as Kojak or XQuest</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element name="psm_annotation_types">
					<xs:annotation>
						<xs:documentation>The types of annotations this search program has for PSMs--such as quality scores.</xs:documentation>
					</xs:annotation>
					<xs:complexType>
						<xs:sequence>
							<xs:element ref="filterable_psm_annotation_types"/>
							<xs:element ref="descriptive_psm_annotation_types" minOccurs="0"/>
						</xs:sequence>
					</xs:complexType>
				</xs:element>
				<xs:element name="psm_per_peptide_annotation_types" minOccurs="0">
					<xs:annotation>
						<xs:documentation>Optional. If this scoring algorithm reports a separate score for each of two cross-linked peptides in a PSM, this element should include the types of scores/annotations to be associated with each peptide. Note that if per-peptide annotations are reported for cross-linked peptides, it is expected that the same annotation is present for non cross-linked peptides as well.</xs:documentation>
					</xs:annotation>
					<xs:complexType>
						<xs:sequence>
							<xs:element ref="filterable_psm_per_peptide_annotation_types">
								<xs:annotation>
									<xs:documentation>If the search program supports it, the collection of numeric annotation types (such as scores) that may be associated with a each separate peptide that is part that are part of a single PSM (e.g. two cross-linked PSMs)</xs:documentation>
								</xs:annotation>
							</xs:element>
							<xs:element ref="descriptive_psm_per_peptide_annotation_types" minOccurs="0">
								<xs:annotation>
									<xs:documentation>Optional. The collection of annotation types associated with each separate peptide that is part of a single PSM that are descriptive or not meant to be used to filter based on quality</xs:documentation>
								</xs:annotation>
							</xs:element>
						</xs:sequence>
					</xs:complexType>
				</xs:element>
				<xs:element name="reported_peptide_annotation_types" minOccurs="0">
					<xs:annotation>
						<xs:documentation>If the search program reports a peptide-level score for each distinct peptide (e.g., peptide-level scores in Percolator), these are the types of scores/annotations that are reported at that level.</xs:documentation>
					</xs:annotation>
					<xs:complexType>
						<xs:sequence>
							<xs:element ref="filterable_peptide_annotation_types" minOccurs="0"/>
							<xs:element ref="descriptive_peptide_annotation_types" minOccurs="0"/>
						</xs:sequence>
					</xs:complexType>
				</xs:element>
			</xs:sequence>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the search program, e.g., "Kojak" or "XQuest" Used to reference this search program in this document.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="display_name">
				<xs:annotation>
					<xs:documentation>If present, the value for this attribute will be used to display the name of this program in the web site. If abscent the value for "name" will be used.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="version" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The version of the above program that was run to generate these results.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="description" type="xs:string">
				<xs:annotation>
					<xs:documentation>(Optional) A description of the search program.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="reported_peptides">
		<xs:annotation>
			<xs:documentation>All the reported peptides from a PSM search</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="reported_peptide" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>A reported peptide. See the annotation for "reported_peptides" for a detailed explanation.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="reported_peptide">
		<xs:annotation>
			<xs:documentation>A peptide from a peptide spectrum match, may be unlinked, looplinked, or a crosslinked pair of peptides.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="peptides">
					<xs:annotation>
						<xs:documentation>A list of peptides associated with a reported peptide in a PSM search. For example, a reported peptide that is a crosslinked pair of peptides will have two peptides associated with it. A looplink or unlink will have one.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="psms">
					<xs:annotation>
						<xs:documentation>The PSMs found for this reported peptide in this search.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element name="reported_peptide_annotations" minOccurs="0">
					<xs:annotation>
						<xs:documentation>Any annotations generated at the reported-peptide-level for this reported peptide.</xs:documentation>
					</xs:annotation>
					<xs:complexType>
						<xs:sequence>
							<xs:element ref="filterable_reported_peptide_annotations" minOccurs="0"/>
							<xs:element ref="descriptive_reported_peptide_annotations" minOccurs="0"/>
						</xs:sequence>
					</xs:complexType>
				</xs:element>
			</xs:sequence>
			<xs:attribute name="reported_peptide_string" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The string that uniquely identifies this peptide (or peptides) linked at these specific positions (if any), with the given positions and respective masses of modifications (if any). See the annotation of "reported_peptides" for more details.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="type" type="link_type" use="required">
				<xs:annotation>
					<xs:documentation>The type of reported peptide, either "crosslink" for two peptides linked by a crosslinker, "looplink" for a single peptide connected by two ends of a crosslinker, or "unlinked" for peptides that are neither of these.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_reported_peptide_annotations">
		<xs:annotation>
			<xs:documentation>The numeric annotations (usually scores of some type) assigned to a reported peptide in the experiment that are meant to describe the likelihood of the peptide being present given all PSMs. Intended to be used to filter out poor matches. Examples include p-value, q-value, PEP</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="filterable_reported_peptide_annotation" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>An annotation or score assigned to a reported peptide</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_reported_peptide_annotations">
		<xs:annotation>
			<xs:documentation>The descriptive annotations (string or numeric) assigned to a reported peptide in the experiment. Not intended to be used to filter out poor matches.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="descriptive_reported_peptide_annotation" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>An annotation or score assigned to a reported peptide</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_reported_peptide_annotation">
		<xs:annotation>
			<xs:documentation>A numeric annotation (usually a score of some type) assigned to a reported peptide in the experiment that are meant to describe the likelihood of the peptide being present given all PSMs. Intended to be used to filter out poor matches. Examples include p-value, q-value, PEP</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="search_program" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the search program, must be found in the search_programs element</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="annotation_name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the peptide annotation type (e.g., "q-value"). Must be found in the possible score types for the above search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="value" type="xs:decimal" use="required">
				<xs:annotation>
					<xs:documentation>The value or score assigned for the given annotation type from the given search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_reported_peptide_annotation">
		<xs:annotation>
			<xs:documentation>A descriptive annotations (string or numeric) assigned to a reported peptide in the experiment. Not intended to be used to filter out poor matches.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="search_program" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the search program, must be found in the search_programs element</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="annotation_name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the peptide annotation type (e.g., "q-value"). Must be found in the possible score types for the above search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="value" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The value or score assigned for the given annotation type from the given search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="peptides">
		<xs:annotation>
			<xs:documentation>A list of peptides associated with a reported peptide in a PSM search.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="peptide" maxOccurs="2">
					<xs:annotation>
						<xs:documentation>A peptide parsed from the reported peptide. Unlinked and looplinked reported peptides will result in one peptide. Crosslinked reported peptides will result in two peptides.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="peptide">
		<xs:annotation>
			<xs:documentation>A peptide parsed from the reported peptide.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="modifications" minOccurs="0"/>
				<xs:element ref="linked_positions" minOccurs="0">
					<xs:annotation>
						<xs:documentation>The list of linked positions in a peptide. Do not include monolinks here, they are listed as modifications in the modifications element. Only list looplinked and crosslinked positions.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element name="peptide_isotope_labels" minOccurs="0">
					<xs:annotation>
						<xs:documentation>If this peptide is isotopically labelled, this element contains the label present on the peptide.</xs:documentation>
					</xs:annotation>
					<xs:complexType>
						<xs:sequence>
							<xs:element name="peptide_isotope_label" minOccurs="0">
								<xs:annotation>
									<xs:documentation>The isotope label for for this peptide. </xs:documentation>
								</xs:annotation>
								<xs:complexType>
									<xs:attribute name="label" use="required">
										<xs:annotation>
											<xs:documentation>The isotope label. Should be one of: 13C, 15N, 18O, or 2H</xs:documentation>
										</xs:annotation>
									</xs:attribute>
								</xs:complexType>
							</xs:element>
						</xs:sequence>
					</xs:complexType>
				</xs:element>
			</xs:sequence>
			<xs:attribute name="sequence" type="sequence_type" use="required">
				<xs:annotation>
					<xs:documentation>The sequence of the peptide, without any special characters.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="unique_id">
				<xs:annotation>
					<xs:documentation>Optional. This is a string that uniquely identifies this peptide in this reported peptide (e.g. "1" and "2" for cross-linked peptides.) IMPORTANT: If this file contains separate PSM-level annotations for each peptide in a single cross-linked PSM (e.g. a XCorr for each linked peptide), then this attribute is required, and is referenced when reporting the PSM-level score/annotation for each peptide. </xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="matched_proteins">
		<xs:annotation>
			<xs:documentation>Contains the protein sequences matched by peptides found in the experiment, and associated names, descriptions, and taxonomy IDs for those proteins. This will typically come from the FASTA file used to search the data. Only the proteins matched by the data should be included, not all the proteins in the FASTA file.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="protein" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>A protein found in the experiment.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="protein">
		<xs:annotation>
			<xs:documentation>A protein found in the experiment.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="protein_annotation" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>Annotations for this sequence from the sequence database used to search the data. Since there may be multiple annotations for a single sequence (as separated by control-A characters in a FASTA file), each annotation should appear hear associated with a single protein.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element name="protein_isotope_labels" minOccurs="0">
					<xs:annotation>
						<xs:documentation>Optional. If this is an isotopically labelled protein sequence, the isotope label for the protein should be listed here. Note that non-labeled and labeld proteins should be separate protein entries eventhough their protein sequences are the same.</xs:documentation>
					</xs:annotation>
					<xs:complexType>
						<xs:sequence>
							<xs:element name="protein_isotope_label" minOccurs="0">
								<xs:annotation>
									<xs:documentation>The isotope label for this protein.</xs:documentation>
								</xs:annotation>
								<xs:complexType>
									<xs:attribute name="label" use="required">
										<xs:annotation>
											<xs:documentation>The isotope label. Should be one of: 13C, 15N, 18O, or 2H</xs:documentation>
										</xs:annotation>
									</xs:attribute>
								</xs:complexType>
							</xs:element>
						</xs:sequence>
					</xs:complexType>
				</xs:element>
			</xs:sequence>
			<xs:attribute name="sequence" type="sequence_type" use="required">
				<xs:annotation>
					<xs:documentation>The amino acid sequence for this protein.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="protein_annotation">
		<xs:annotation>
			<xs:documentation>Annotations for this sequence from the sequence database used to search the data. Since there may be multiple annotations for a single sequence (as separated by control-A characters in a FASTA file), each annotation should appear hear associated with a single protein.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of this protein, used to report protein names in the ProXL web application. Typically, this is the name for a FASTA entry in a FASTA file in the form of: >NAME DESCRIPTION</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="description" type="xs:string">
				<xs:annotation>
					<xs:documentation>The description of the protein.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="ncbi-taxonomy-id" type="xs:integer">
				<xs:annotation>
					<xs:documentation>The NCBI taxonomy ID for this protein. E.g. 9606 for human, or 4932 for S. cerevisiae.  Used for taxonomy-level filtering in the ProXL web application. If not present, or set to 0, the species name will appear as "Species Unknown" in the web application.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="modifications">
		<xs:annotation>
			<xs:documentation>The modifications identified for a peptide.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="modification" minOccurs="0" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="modification">
		<xs:annotation>
			<xs:documentation>A mass modification of a residue. E.g., a PTM</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="mass" type="mass" use="required">
				<xs:annotation>
					<xs:documentation>The minoisotopic mass, in daltons of the modification. E.g. 79.966331 for phosphorylation</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="position" type="peptide_position" use="optional">
				<xs:annotation>
					<xs:documentation>The position in the peptide (starting at 1) of the modification.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="is_n_terminal" type="xs:boolean">
				<xs:annotation>
					<xs:documentation>Set this to true if this mod is at the n-terminus of the peptide. Omit if this is not an n-terminal mod.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="is_c_terminal" type="xs:boolean">
				<xs:annotation>
					<xs:documentation>Set this to true if this mod is at the c-terminus of the peptide. Omit if this is not a c-terminal mod.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="isMonolink" type="xs:boolean">
				<xs:annotation>
					<xs:documentation>If this modification is a monolink, this is set to "true"</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="linked_position">
		<xs:annotation>
			<xs:documentation>A linked position in a peptide.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="position" type="peptide_position" use="required"/>
		</xs:complexType>
	</xs:element>
	<xs:element name="linked_positions">
		<xs:annotation>
			<xs:documentation>The list of linked positions in a peptide.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="linked_position" minOccurs="0" maxOccurs="2">
					<xs:annotation>
						<xs:documentation>A linked position in a peptide. A looplinked peptide will have two linked positions. A crosslinked peptide will have one (on each linked peptide). The first residue is position 1.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_psm_annotation_type">
		<xs:annotation>
			<xs:documentation>A numeric type of psm annotation (such as a score like a q-value for XCorr) that is meant to be used to filter PSMs based on quality of matches. Try to limit indexed annotation types to a minimum.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name for this score, such as "q-value". Displayed on web site.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="description" type="xs:string">
				<xs:annotation>
					<xs:documentation>A description for what this score is. Will be shown in ProXL as a tooltip for this score.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="filter_direction" type="filter_direction_type" use="required">
				<xs:annotation>
					<xs:documentation>The direction a filterable annotation type is sorted in.  If set to "below", attributes with lower values are considered more significant (such as in the case of p-values). If set to "above", attributes with higher values are considered more significant (such as in the case of XCorr).</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="default_filter" type="xs:boolean">
				<xs:annotation>
					<xs:documentation>If present, and if set to true, this score will be used as a default PSM-level filter when viewing data. Use sparingly.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="default_filter_value" type="xs:decimal">
				<xs:annotation>
					<xs:documentation>This is the default cutoff value for this score when filtering data. If default_filter is set to true, this value must be set.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_psm_per_peptide_annotation_type">
		<xs:annotation>
			<xs:documentation>A numeric type of psm annotation (such as a score like a q-value for XCorr) that is meant to be used to filter PSMs based on quality of matches. Try to limit indexed annotation types to a minimum.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name for this score, such as "q-value". Displayed on web site.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="description" type="xs:string">
				<xs:annotation>
					<xs:documentation>A description for what this score is. Will be shown in ProXL as a tooltip for this score.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="filter_direction" type="filter_direction_type" use="required">
				<xs:annotation>
					<xs:documentation>The direction a filterable annotation type is sorted in.  If set to "below", attributes with lower values are considered more significant (such as in the case of p-values). If set to "above", attributes with higher values are considered more significant (such as in the case of XCorr).</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_peptide_annotation_type">
		<xs:annotation>
			<xs:documentation>A numeric type of peptide annotation (such as a score like a q-value for XCorr) that is meant to be used to filter peptides based on quality of matches. Try to limit indexed annotation types to a minimum.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name for this score, such as "q-value". Displayed on web site.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="description" type="xs:string">
				<xs:annotation>
					<xs:documentation>A description for what this score is. Will be shown in ProXL as a tooltip for this score.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="filter_direction" type="filter_direction_type" use="required">
				<xs:annotation>
					<xs:documentation>The direction a filterable annotation type is sorted in.  If set to "below", attributes with lower values are considered more significant (such as in the case of p-values). If set to "above", attributes with higher values are considered more significant (such as in the case of XCorr).</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="default_filter" type="xs:boolean">
				<xs:annotation>
					<xs:documentation>If present, and if set to true, this score will be used as a default peptide-level filter when viewing data. Use sparingly.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="default_filter_value" type="xs:decimal">
				<xs:annotation>
					<xs:documentation>This is the default cutoff value for this score when filtering data. If default_filter is set to true, this value must be set.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_psm_annotation_type">
		<xs:annotation>
			<xs:documentation>An annotation of any type (numeric or string) that is not meant to be used as a possible filter for filtering out PSMs based on quality of matches. Example could include "scan number" or "charge"</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name for this score, such as "q-value". Displayed on web site.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="description" type="xs:string">
				<xs:annotation>
					<xs:documentation>A description for what this score is. Will be shown in ProXL as a tooltip for this score.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_psm_per_peptide_annotation_type">
		<xs:annotation>
			<xs:documentation>An annotation of any type (numeric or string) that is not meant to be used as a possible filter for filtering out PSMs based on quality of matches. Example could include "scan number" or "charge"</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name for this score, such as "q-value". Displayed on web site.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="description" type="xs:string">
				<xs:annotation>
					<xs:documentation>A description for what this score is. Will be shown in ProXL as a tooltip for this score.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_peptide_annotation_type">
		<xs:annotation>
			<xs:documentation>An annotation of any type (numeric or string) that is not meant to be used as a possible filter for filtering out peptides based on quality of matches.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name for this score, such as "q-value". Displayed on web site.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="description" type="xs:string">
				<xs:annotation>
					<xs:documentation>A description for what this score is. Will be shown in ProXL as a tooltip for this score.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:simpleType name="filter_direction_type">
		<xs:annotation>
			<xs:documentation>The direction a filterable annotation type is sorted in.  If set to "below", attributes with lower values are considered more significant (such as in the case of p-values). If set to "above", attributes with higher values are considered more significant (such as in the case of XCorr).</xs:documentation>
		</xs:annotation>
		<xs:restriction base="xs:string">
			<xs:enumeration value="above"/>
			<xs:enumeration value="below"/>
		</xs:restriction>
	</xs:simpleType>
	<xs:simpleType name="sequence_type">
		<xs:annotation>
			<xs:documentation>Reported sequence strings must be of this type, includes only amino acid sequence</xs:documentation>
		</xs:annotation>
		<xs:restriction base="xs:string">
		</xs:restriction>
	</xs:simpleType>
	<xs:simpleType name="amino_acid_residue">
		<xs:annotation>
			<xs:documentation>The single letter code for an amino acid. E.g., "K" for lysine.</xs:documentation>
		</xs:annotation>
		<xs:restriction base="xs:string">
		</xs:restriction>
	</xs:simpleType>
	<xs:simpleType name="link_type">
		<xs:annotation>
			<xs:documentation>A peptide link type must be this type, limits to 'crosslink' or 'looplink' or 'unlinked' 'unlinked means neither a crosslink nor looplink</xs:documentation>
		</xs:annotation>
		<xs:restriction base="xs:string">
			<xs:enumeration value="crosslink"/>
			<xs:enumeration value="looplink"/>
			<xs:enumeration value="unlinked"/>
		</xs:restriction>
	</xs:simpleType>
	<xs:simpleType name="protein_terminus_designation">
		<xs:annotation>
			<xs:documentation>The terminus of a protein, either "n" or "c".</xs:documentation>
		</xs:annotation>
		<xs:restriction base="xs:string">
			<xs:enumeration value="n"/>
			<xs:enumeration value="c"/>
		</xs:restriction>
	</xs:simpleType>
	<xs:simpleType name="peptide_position">
		<xs:annotation>
			<xs:documentation>A position in a peptide, where 1 is the first position.</xs:documentation>
		</xs:annotation>
		<xs:restriction base="xs:positiveInteger">
			<xs:minInclusive value="1"/>
		</xs:restriction>
	</xs:simpleType>
	<xs:simpleType name="mass">
		<xs:annotation>
			<xs:documentation>Mass represented in daltons</xs:documentation>
		</xs:annotation>
		<xs:restriction base="xs:decimal"/>
	</xs:simpleType>
	<xs:element name="filterable_psm_annotation_types">
		<xs:annotation>
			<xs:documentation>The collection of numeric annotation types (such as scores) that may be associated with a PSM that are meant to be used to filter based on quality</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="filterable_psm_annotation_type" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_psm_per_peptide_annotation_types">
		<xs:annotation>
			<xs:documentation>The collection of numeric annotation types (such as scores) that may be associated with a PSM that are meant to be used to filter based on quality</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="filterable_psm_per_peptide_annotation_type" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_peptide_annotation_types">
		<xs:annotation>
			<xs:documentation>The collection of numeric annotation types (such as scores) that may be associated with a reported peptide that are meant to be used to filter based on quality</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="filterable_peptide_annotation_type" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_psm_annotation_types">
		<xs:annotation>
			<xs:documentation>The collection of annotation types (string or numeric) associated with PSMs that are descriptive or not meant to be used to filter based on quality</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="descriptive_psm_annotation_type" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_psm_per_peptide_annotation_types">
		<xs:annotation>
			<xs:documentation>The collection of annotation types (string or numeric) associated with PSMs that are descriptive or not meant to be used to filter based on quality</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="descriptive_psm_per_peptide_annotation_type" minOccurs="0" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_peptide_annotation_types">
		<xs:annotation>
			<xs:documentation>The collection of annotation types (string or numeric) associated with PSMs that are descriptive or not meant to be used to filter based on quality</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="descriptive_peptide_annotation_type" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>An annotation of any type (numeric or string) that is not meant to be used as a possible filter for filtering out peptides based on quality of matches. Examples include calculated mass or validation status.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_psm_annotation">
		<xs:annotation>
			<xs:documentation>A numeric annotation (usually a score) assigned to a peptide spectrum match (PSM) that is meant to describe the quality of the match and may be used to filter out poor matches. Examples include p-value, q-value, XCorr</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="search_program" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the search program, must be found in the search_programs element</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="annotation_name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the psm annotation type (e.g., "q-value"). Must be found in the possible score types for the above search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="value" type="xs:decimal" use="required">
				<xs:annotation>
					<xs:documentation>The value (ie, score) assigned for the given annotationtype from the given search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_psm_per_peptide_annotation">
		<xs:annotation>
			<xs:documentation>A numeric annotation (usually a score) assigned to a peptide spectrum match (PSM) that is meant to describe the quality of the match and may be used to filter out poor matches. Examples include p-value, q-value, XCorr</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="search_program" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the search program, must be found in the search_programs element</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="annotation_name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the psm annotation type (e.g., "q-value"). Must be found in the possible score types for the above search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="value" type="xs:decimal" use="required">
				<xs:annotation>
					<xs:documentation>The value (ie, score) assigned for the given annotationtype from the given search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_psm_annotation">
		<xs:annotation>
			<xs:documentation>A descriptive annotation (numeric or string) assigned to a peptide spectrum match (PSM) that is not meant to describe the quality of the match or filter out poor matches. Examples include scan number of precursor charge state.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="search_program" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the search program, must be found in the search_programs element</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="annotation_name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the psm annotation type (e.g., "q-value"). Must be found in the possible score types for the above search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="value" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The value (ie, score) assigned for the given annotationtype from the given search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_psm_per_peptide_annotation">
		<xs:annotation>
			<xs:documentation>A descriptive annotation (numeric or string) assigned to a peptide spectrum match (PSM) that is not meant to describe the quality of the match or filter out poor matches. Examples include scan number of precursor charge state.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="search_program" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the search program, must be found in the search_programs element</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="annotation_name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of the psm annotation type (e.g., "q-value"). Must be found in the possible score types for the above search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="value" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The value (ie, score) assigned for the given annotationtype from the given search program</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_psm_annotations">
		<xs:annotation>
			<xs:documentation>The numeric annotations (usually scores of some type) assigned to the peptide spectrum match (PSM) that are meant to describe the quality of the match and may be used to filter out poor matches. Examples include p-value, q-value, XCorr</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="filterable_psm_annotation" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="filterable_psm_per_peptide_annotations">
		<xs:annotation>
			<xs:documentation>The numeric annotations (usually scores of some type) assigned to the peptide spectrum match (PSM) that are meant to describe the quality of the match and may be used to filter out poor matches. Examples include p-value, q-value, XCorr</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="filterable_psm_per_peptide_annotation" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_psm_annotations">
		<xs:annotation>
			<xs:documentation>The descriptive annotations (numeric or string) assigned to the peptide spectrum match (PSM) that are not meant to describe the quality of the match or filter out poor matches. Examples include scan number of precursor charge state.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="descriptive_psm_annotation" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>A descriptive annotation (numeric or string) assigned to a peptide spectrum match (PSM) that is not meant to describe the quality of the match or filter out poor matches.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="descriptive_psm_per_peptide_annotations">
		<xs:annotation>
			<xs:documentation>The descriptive annotations (numeric or string) assigned to the peptide spectrum match (PSM) that are not meant to describe the quality of the match or filter out poor matches. Examples include scan number of precursor charge state.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="descriptive_psm_per_peptide_annotation" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>A descriptive annotation (numeric or string) assigned to a peptide spectrum match (PSM) that is not meant to describe the quality of the match or filter out poor matches.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="psm">
		<xs:annotation>
			<xs:documentation>A PSM</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="filterable_psm_annotations"/>
				<xs:element ref="descriptive_psm_annotations" minOccurs="0">
					<xs:annotation>
						<xs:documentation>The descriptive annotations (numeric or string) assigned to the peptide spectrum match (PSM) that are not meant to describe the quality of the match or filter out poor matches.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element name="per_peptide_annotations" minOccurs="0">
					<xs:annotation>
						<xs:documentation>Optional. Must use when there are separate data to report for each peptide in a cross-link PSM. Examples include cleavable cross-linker experiments or whenever a separate score score is provided for each peptide in a cross-link (e.g. an XCorr for each separate peptide ID).</xs:documentation>
					</xs:annotation>
					<xs:complexType>
						<xs:sequence>
							<xs:element name="psm_peptide" maxOccurs="2">
								<xs:annotation>
									<xs:documentation>Each distinct peptide found in this PSM. For cross-linked PSMs, there will be 2 psm_peptide elements.</xs:documentation>
								</xs:annotation>
								<xs:complexType>
									<xs:sequence>
										<xs:element ref="filterable_psm_per_peptide_annotations">
											<xs:annotation>
												<xs:documentation>The numeric annotations (usually scores of some type) assigned to a specific peptide that is part of a peptide spectrum match (PSM) that are meant to describe the quality of the match and may be used to filter out poor matches. Examples include p-value, q-value, XCorr</xs:documentation>
											</xs:annotation>
										</xs:element>
										<xs:element ref="descriptive_psm_per_peptide_annotations" minOccurs="0">
											<xs:annotation>
												<xs:documentation>The descriptive annotations (numeric or string) assigned to a specific peptide that is part of the peptide spectrum match (PSM) that are not meant to describe the quality of the match or filter out poor matches. Examples include scan number of precursor charge state.</xs:documentation>
											</xs:annotation>
										</xs:element>
									</xs:sequence>
									<xs:attribute name="unique_id" use="required">
										<xs:annotation>
											<xs:documentation>The unique id given to this peptide for this reported_peptide (e.g. "1" or "2". Recall that a reported_peptide is the identified ion, so will contain two separate peptides in the case of a cross-linked peptide. See the peptide element under the peptides element for a reported_peptide.</xs:documentation>
										</xs:annotation>
									</xs:attribute>
									<xs:attribute name="scan_number" type="xs:positiveInteger">
										<xs:annotation>
											<xs:documentation>(Optional) The scan number used to look up the scan in the given file. For cleavable cross-linker data, this will likely be the scan number of the MS3 scan for this peptide.</xs:documentation>
										</xs:annotation>
									</xs:attribute>
									<xs:attribute name="scan_file_name" type="xs:string">
										<xs:annotation>
											<xs:documentation>(Optional) The filename of the file in which this scan may be found. </xs:documentation>
										</xs:annotation>
									</xs:attribute>
									<xs:attribute name="linker_mass" type="mass">
										<xs:annotation>
											<xs:documentation>(Optional) Must be included for a cleavable cross-linker experiment. This is the mass of the side of the cross-linker that is bound to this peptide.</xs:documentation>
										</xs:annotation>
									</xs:attribute>
								</xs:complexType>
							</xs:element>
						</xs:sequence>
					</xs:complexType>
				</xs:element>
			</xs:sequence>
			<xs:attribute name="scan_file_name" type="xs:string">
				<xs:annotation>
					<xs:documentation>(Optional) The filename of the file in which this scan may be found. Will be used to locate the spectrum associated with this PSM when data are uploaded.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="scan_number" type="xs:positiveInteger">
				<xs:annotation>
					<xs:documentation>(Optional) The scan number used to look up the scan in the given file. Will be used to locate the spectrum associated with this PSM when data are uploaded.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="precursor_charge" type="xs:positiveInteger" use="required">
				<xs:annotation>
					<xs:documentation>The predicted charge of the precursor ion. Used for accurately annotating spectra in web interface.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="precursor_retention_time" type="xs:decimal">
				<xs:annotation>
					<xs:documentation>(Optional) The retention time (in seconds) of the precursor ion. If not present, the retention times from the optionally uploaded mz data files (e.g. mzML will be used).</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="precursor_m_z" type="xs:decimal">
				<xs:annotation>
					<xs:documentation>(Optional) The observed m/z of the precursor ion. If not present, the measured m/z from the optionally uploaded mz data files (e.g. mzML will be used).</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="linker_mass" type="mass">
				<xs:annotation>
					<xs:documentation>(Optional) If this PSM identified a looplinked peptide or crosslinked peptidie pair, this is the mass used for the crosslinker in this identification. Used for accurately annotating spectra in web interface.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="psms">
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="psm" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="search_programs">
		<xs:annotation>
			<xs:documentation>The search program(s) used to generate the PSMs and associated scores in this search. More than one search program may be used, as in the case of using Percolator to assign q-values for data search with Kojak. Percolator and Kojak would each have an entry here. However, all PSMs must be annotated by all search programs and all search programs must refer to reported peptides using the same string.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_program" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="linker">
		<xs:annotation>
			<xs:documentation>A chemical crosslinker used in this experiment.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="monolink_masses" minOccurs="0">
					<xs:annotation>
						<xs:documentation>The list of possible masses that monolinks may be reported as having in this search. </xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="crosslink_masses" minOccurs="0">
					<xs:annotation>
						<xs:documentation>The list of possible masses of this linker in this experiment. Required if cross-links or loop-links were searched for and if spectra are being imported.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="linked_ends" minOccurs="0">
					<xs:annotation>
						<xs:documentation>Define where each end of the cross-linker may bind with proteins.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
			<xs:attribute name="name" type="xs:string" use="required">
				<xs:annotation>
					<xs:documentation>The name of this crosslinker (e.g., "dss" or "edc").</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="spacer_arm_length" type="xs:decimal">
				<xs:annotation>
					<xs:documentation>Optional. The length, in angstroms, of the spacer arm of the cross-linker. E.g., 10.3 for DSSO</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="monolink_mass">
		<xs:complexType>
			<xs:attribute name="mass" type="mass" use="required">
				<xs:annotation>
					<xs:documentation>The minoisotopic mass, in daltons.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="crosslink_mass">
		<xs:complexType>
			<xs:attribute name="mass" type="mass" use="required">
				<xs:annotation>
					<xs:documentation>The minoisotopic mass, in daltons.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="chemical_formula" type="xs:string">
				<xs:annotation>
					<xs:documentation>Optional. The chemical formula corresponding to this mass for this cross-linker. E.g., C6H6O3S for DSSO. Required for some reports in proxl.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="monolink_masses">
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="monolink_mass" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>A mass used for this linker in this search when only one end of the linker is bound to a peptide. Use the same precision or better as when monolinks are reported as dynamic mods on peptides.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="crosslink_masses">
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="crosslink_mass" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>A full-length mass used for this linker in this search when both ends are bound to peptides--such as in crosslinks or looplinks</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="cleaved_crosslink_mass" minOccurs="0" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="linkers">
		<xs:annotation>
			<xs:documentation>The chemical linkers used in this experiment. There must be at least one.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="linker" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="configuration_file">
		<xs:annotation>
			<xs:documentation>A configuration or properties file associated with the search. E.g. Kojak.conf</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="file_content"/>
			</xs:sequence>
			<xs:attribute name="search_program" use="required"/>
			<xs:attribute name="file_name" use="required"/>
		</xs:complexType>
	</xs:element>
	<xs:element name="file_content" type="xs:base64Binary">
		<xs:annotation>
			<xs:documentation>The contents of a file, stored as base64</xs:documentation>
		</xs:annotation>
	</xs:element>
	<xs:element name="configuration_files">
		<xs:annotation>
			<xs:documentation>The configuration files associated with the search.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="configuration_file" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="static_modification">
		<xs:annotation>
			<xs:documentation>A mass shift that is always applied to the given amino acid when calculating theoretical peptide or ion masses</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="amino_acid" use="required">
				<xs:annotation>
					<xs:documentation>The one-letter code for the amino acid with the static modification.</xs:documentation>
				</xs:annotation>
				<xs:simpleType>
					<xs:restriction base="xs:string">
						<xs:length value="1"/>
					</xs:restriction>
				</xs:simpleType>
			</xs:attribute>
			<xs:attribute name="mass_change" type="mass" use="required">
				<xs:annotation>
					<xs:documentation>The mass shift applied to all instances of the above amino acid.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="static_modifications">
		<xs:annotation>
			<xs:documentation>A collection of mass shifts that are always applied to the respective amino acids when calculating theoretical peptide or ion masses</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="static_modification" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="decoy_labels">
		<xs:annotation>
			<xs:documentation>A list of decoy labels (strings that are prefixes on protein names that indicate the given protein is a decoy (e.g., "DECOY_")</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="decoy_label" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="decoy_label">
		<xs:annotation>
			<xs:documentation>A decoy label is used to identify a protein to which a peptide is matched during protein inference as a decoy protein (used to calculate FDR statistics). Matches to decoy proteins will be ignored and not inserted into the database as a matched protein for peptides. Example: <decoy_label prefix="DECOY_"/>
			</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="prefix" type="xs:string">
				<xs:annotation>
					<xs:documentation>Any protein name beginning with the supplied string will be considered a decoy protein and ignored.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="search_program_info">
		<xs:annotation>
			<xs:documentation>Information relating to the search program(s) used to generate the peptide spectrum matches and scores/annotations for PSMs and peptides.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_programs"/>
				<xs:element ref="annotation_cutoffs_on_import" minOccurs="0">
					<xs:annotation>
						<xs:documentation>DEPRECATED - Will be ignored by proxl importer. May be removed in future schema release.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="default_visible_annotations"/>
				<xs:element ref="annotation_sort_order" minOccurs="0">
					<xs:annotation>
						<xs:documentation>The sort order applied to filterable annotations when presenting the default listing of PSMs or peptides. The order will be in the order the annotations appear here. E.g. if "q-value" is listed first and "p-value" is listed next, the items will be sorted on q-value first, then, sorted on p-value in the case of ties. Note: all annotation types listed here must be filterable. If this element is not present, no default sorting will be applied to listed PSMs or peptides in ProXL.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="annotation_cutoffs_on_import">
		<xs:annotation>
			<xs:documentation>The cutoffs used to exclude data while importing.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="psm_annotation_cutoffs_on_import">
					<xs:annotation>
						<xs:documentation>DEPRECATED - Will be ignored by proxl importer. May be removed in future schema release.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="reported_peptide_annotation_cutoffs_on_import" minOccurs="0">
					<xs:annotation>
						<xs:documentation>DEPRECATED - Will be ignored by proxl importer. May be removed in future schema release.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="psm_annotation_cutoffs_on_import">
		<xs:annotation>
			<xs:documentation>The PSM-level annotation cutoffs used to exclude data while importing.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_annotation_cutoff" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="reported_peptide_annotation_cutoffs_on_import">
		<xs:annotation>
			<xs:documentation>The peptide-level annotation cutoffs used to exclude data while importing.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_annotation_cutoff" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="search_annotation_cutoff">
		<xs:annotation>
			<xs:documentation>A reference to an annotation in a search program.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="search_program" use="required">
				<xs:annotation>
					<xs:documentation>The name of a search program, must match a name for a search program in this document. E.g., "kojak"</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="annotation_name" use="required">
				<xs:annotation>
					<xs:documentation>The name of the annotation, must match a name of an annotation for the search_program. E.g., "q-value"</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="cutoff_value" type="xs:decimal" use="required">
				<xs:annotation>
					<xs:documentation>This is the cutoff value for this score for excluding data while importing into Proxl.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="default_visible_annotations">
		<xs:annotation>
			<xs:documentation>The annotations from the search program(s) that are visible by default in the ProXL web app.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="visible_psm_annotations"/>
				<xs:element ref="visible_psm_per_peptide_annotations" minOccurs="0">
					<xs:annotation>
						<xs:documentation>Optional. The PSM-level annotations that are associated individual with each peptide (for example, each separate peptide in a cross-link) that should be visible by default when viewing PSM-level data.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="visible_reported_peptide_annotations" minOccurs="0">
					<xs:annotation>
						<xs:documentation>Optional. The peptide-level annotations that are visible by default.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="visible_psm_annotations">
		<xs:annotation>
			<xs:documentation>The PSM-level annotations that are visible by default when viewing PSM-level data.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_annotation" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="visible_psm_per_peptide_annotations">
		<xs:annotation>
			<xs:documentation>The PSM-level annotations that are visible by default when viewing PSM-level data.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_annotation" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="visible_reported_peptide_annotations">
		<xs:annotation>
			<xs:documentation>The peptide-level annotations that are visible by default.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_annotation" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="search_annotation">
		<xs:annotation>
			<xs:documentation>A reference to an annotation in a search program.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="search_program" use="required">
				<xs:annotation>
					<xs:documentation>The name of a search program, must match a name for a search program in this document. E.g., "kojak"</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="annotation_name" use="required">
				<xs:annotation>
					<xs:documentation>The name of the annotation, must match a name of an annotation for the search_program. E.g., "q-value"</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="annotation_sort_order">
		<xs:annotation>
			<xs:documentation>The sort order applied to filterable annotation when presenting the default listing of PSMs or peptides. The order will be in the order the annotations appear here. E.g. if "q-value" is listed first and "p-value" is listed next, the items will be sorted on q-value first, then, sorted on p-value in the case of ties. Note: all annotation types listed here must be filterable. If this element is not present, no default sorting will be applied to listed PSMs or peptides in ProXL.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="psm_annotation_sort_order" minOccurs="0"/>
				<xs:element ref="reported_peptide_annotation_sort_order" minOccurs="0"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="psm_annotation_sort_order">
		<xs:annotation>
			<xs:documentation>The sort order applied to filterable PSM annotations when presenting the default listing of PSMs. The order will be in the order the annotations appear here. E.g. if "q-value" is listed first and "p-value" is listed next, the PSMs will be sorted on q-value first, then, sorted on p-value in the case of ties. Note: all annotation types listed here must be filterable. If this element is not present, no default sorting will be applied to listed PSMs.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_annotation" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="reported_peptide_annotation_sort_order">
		<xs:annotation>
			<xs:documentation>The sort order applied to filterable peptide annotations when presenting the default listing of peptides. The order will be in the order the annotations appear here. E.g. if "q-value" is listed first and "p-value" is listed next, the peptides will be sorted on q-value first, then, sorted on p-value in the case of ties. Note: all annotation types listed here must be filterable. If this element is not present, no default sorting will be applied to listed peptides.</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="search_annotation" maxOccurs="unbounded"/>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="cleaved_crosslink_mass">
		<xs:annotation>
			<xs:documentation>If the cross-linker cleavable, a mass to be used as a cleaved product of this cross-linker. Ie, the mass of each side of the cross-linker bound to a peptide after cleavage. A cleavable linkers will usually have two of these (but not always).</xs:documentation>
		</xs:annotation>
		<xs:complexType>
			<xs:attribute name="mass" type="mass" use="required">
				<xs:annotation>
					<xs:documentation>The minoisotopic mass, in daltons.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="chemical_formula" type="xs:string">
				<xs:annotation>
					<xs:documentation>Optional. The chemical formula corresponding to this mass for this cross-linker. E.g., C3H2O for one of the cleaved ends of DSSO.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
	<xs:element name="linked_ends">
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="linked_end" minOccurs="2" maxOccurs="2">
					<xs:annotation>
						<xs:documentation>Define where one end of the cross-linker may react with a protein. There must be two of these elements, even if both sides are the same. This must contain at least one residues element or protein_termini element.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="linked_end">
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="residues" minOccurs="0">
					<xs:annotation>
						<xs:documentation>A list of amino acid residues with which this end of the cross-linker may bind.</xs:documentation>
					</xs:annotation>
				</xs:element>
				<xs:element ref="protein_termini" minOccurs="0">
					<xs:annotation>
						<xs:documentation>A list of protein terminal positions to which this cross-linker can bind. E.g. 0 or 1 residues away from n-terminus.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="residue" type="amino_acid_residue"/>
	<xs:element name="residues">
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="residue" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>The single letter code for the residue with which this end of the crosslinker may bind. E.g., K</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="protein_termini">
		<xs:complexType>
			<xs:sequence>
				<xs:element ref="protein_terminus" maxOccurs="unbounded">
					<xs:annotation>
						<xs:documentation>A protein terminus, either n or c, and the number of residues away from that terminus the linker can bind.</xs:documentation>
					</xs:annotation>
				</xs:element>
			</xs:sequence>
		</xs:complexType>
	</xs:element>
	<xs:element name="protein_terminus">
		<xs:complexType>
			<xs:attribute name="terminus_end" type="protein_terminus_designation" use="required">
				<xs:annotation>
					<xs:documentation>The terminal end, either "n" or "c".</xs:documentation>
				</xs:annotation>
			</xs:attribute>
			<xs:attribute name="distance_from_terminus" type="xs:nonNegativeInteger" use="required">
				<xs:annotation>
					<xs:documentation>The distance from the given terminus. 0 indicates the terminus itself. 1 indicates one position removed from the terminus. E.g., in the protein PEPTIDE, a distance of 0 from the c terminus would be E, a distance of 1 would be D.</xs:documentation>
				</xs:annotation>
			</xs:attribute>
		</xs:complexType>
	</xs:element>
</xs:schema>
"""
