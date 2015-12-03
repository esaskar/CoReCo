#!/usr/bin/env python

import os
import sys
import re

#sys.path.append("/usr/local/lib/")

import libsbml

sys.path.append("../model_training_scripts/")
import common

SUM_STOIC_COEFFS = 1  # sum coefficients of both sides: 2A + 2B <=> 3A + 1C ---> 2B <=> 1A + 1C

NOTES = """
<body xmlns="http://www.w3.org/1999/xhtml">
<p>
This model was reconstructed with CoReCo method from protein sequence and phylogeny data. CoReCo is described in Pitkanen, E., Jouhten, P., Hou, J., Syed, M. F., Blomberg, P., Kludas, J., Oja, M., Holm, L., Penttila, M., Rousu, J. and Arvas, M. (2014). Comparative Genome-Scale Reconstruction of Gapless Metabolic Networks for Present and Ancestral Species. PLoS Computational Biology, 10(2), e1003465. doi:10.1371/journal.pcbi.1003465.

CoReCo annotates each reaction with the following attributes: 
<ul>
<li>reco:balanced: True, if reaction is in stoichiometric balance.</li>
<li>reco:posterior: CoReCo phase II probability</li>
<li>reco:cost: CoReCo phase II log-score</li>
<li>reco:level: High, if cost below reconstruction threshold, and Low, if cost above reconstruction threshold. A reaction with level=Low only appears in reconstruction to gapfill a reaction with level=High.</li>
<li>reco:naivep: CoReCo phase I probability.</li>
<li>reco:btscore: CoReCo phase II probability considering only BLAST as input.</li>
<li>reco:gtscore: CoReCo phase II probability considering only GTG as input.</li>
<li>reco:bscore: CoReCo BLAST score</li>
<li>reco:bseq1: Protein sequence in target species with highest BLAST score</li>
<li>reco:bseq2: Protein sequence in UniProt with highest BLAST score</li>
<li>reco:gscore: CoReCo GTG score</li>
<li>reco:gseq1: Protein sequence in target species</li>
<li>reco:gseq2: Protein sequence in GTG database with highest GTG score</li>
</ul>
</p>
</body>
"""

ANNOTATION = """
      <rdf:RDF>
        <rdf:Description rdf:about="#metaid_ymn4_0">
           <dc:description>     
              <p>Metabolic model of %s reconstructed with CoReCo.</p>
          </dc:description>
          <dc:creator rdf:parseType="Resource">
             <rdf:Bag>
              <rdf:li rdf:parseType="Resource">
                <vCard:N rdf:parseType="Resource">
                  <vCard:Family>Pitkanen</vCard:Family>
                  <vCard:Given>Esa</vCard:Given>
                </vCard:N>
                <vCard:EMAIL>esa.pitkanen@helsinki.fi</vCard:EMAIL>
                <vCard:ORG>
                  <vCard:Orgname>University of Helsinki</vCard:Orgname>
                </vCard:ORG>
              </rdf:li>
            </rdf:Bag>
          </dc:creator>
          <dcterms:created rdf:parseType="Resource">
            <dcterms:W3CDTF>2013-01-20T12:00:00Z</dcterms:W3CDTF>
          </dcterms:created>
          <dcterms:modified rdf:parseType="Resource">
            <dcterms:W3CDTF>2013-01-20T12:00:00Z</dcterms:W3CDTF>
          </dcterms:modified>
          <bqbiol:is>
            <rdf:Bag>
              <rdf:li rdf:resource="urn:miriam:taxonomy:%s"/>
            </rdf:Bag>
          </bqbiol:is>
        </rdf:Description>
      </rdf:RDF>
"""

def get_reaction_id(r):
    """R00001_1_rev -> R00001"""
    return r.split("_")[0].replace("#","_")
#    return r.split("_")[0]

def convert_to_SBML(reco, eqns, mol2name, taxon, modelid, modelname, species, sbmlversion, bounds,  rxnnames=None, pathways=None, pathwayasnames=None):

    DEFAULT_COMPARTMENT = "cytosol"
    reKEGGRID = re.compile("R\d\d\d\d\d")
    if sbmlversion == 2:	
        d = libsbml.SBMLDocument(2, 4)
    else:
        d = libsbml.SBMLDocument(3, 1)

    ns = d.getNamespaces()
    ns.add("http://www.w3.org/1999/02/22-rdf-syntax-ns", "rdf")
    ns.add("http://www.w3.org/2001/vcard-rdf/3.0#", "vCard")
    ns.add("http://purl.org/dc/elements/1.1/", "dc")
    ns.add("http://purl.org/dc/terms/", "dcterms")
    ns.add("http://biomodels.net/biology-qualifiers/", "bqbiol")
    ns.add("http://biomodels.net/model-qualifiers/", "bqmodel")
    ns.add("http://www.cs.helsinki.fi/group/sysfys/coreco", "coreco")

# biomodels.net instructions
#Enter all the relevant information you believe is necessary for the curation (reference publication, modifications or clarifications of the model, etc.) either directly into the model file if allowed (for example using the notes elements if your model is under one of the SBML formats), or into the Curation comment text field provided by the form below.
#If you created the model, or collaborated to its creation, and you are not an author of the reference publication, add to the model element a dc:creator annotation containing your data (first and last name, organisation, email), so that your contribution can be acknowledged. Click here to view an example of a dc:creator annotation which you can re-use (skip blue part if already present).
#Choose a meaningful value for the attribute name of the model element. Examples of good model names are NameAuthorYear_Topic_Method, Levchenko2000_MAPK_noScaffold or Edelstein1996_EPSP_AChEvent.
#Check the validity of the model (for example by using this online validator if your model is encoded in the SBML format). All the models undergo a primary XML validity check upon submission anyway, and, as mentioned before, a more thorough testing during the curation phase, but an already valid model is of great help nevertheless!
#If the model was not created directly in SBML, or if it requires a specific software to be simulated adequately, please enter in the Original Model form a URL pointing to the model in the original repository. Refrain from entering a generic URL to the repository itself.

    m = d.createModel()
    m.setName(modelname)
    m.setId(modelid)
    m.setNotes(NOTES)
    m.setAnnotation(ANNOTATION % (species, taxon))

    comp = m.createCompartment()
    comp.setId(DEFAULT_COMPARTMENT)
    #comp.setCompartmentType()
    comp.setConstant(False) # Sets the value of the 'constant' attribute of this Compartment.
    comp.setName(DEFAULT_COMPARTMENT)
    #comp.setOutside() # sid the identifier of a compartment that encloses this one.
    comp.setSBOTerm("SBO:0000290")
    comp.setSize(1)
    comp.setAnnotation("""
<rdf:RDF>
<rdf:Description rdf:about="">
<bqbiol:is>
<rdf:Bag>
<rdf:li rdf:resource="urn:miriam:obo.go:GO%3A0005737" />  
</rdf:Bag>
</bqbiol:is>
</rdf:Description>
</rdf:RDF>
""")

    reactions = set()
    rid2re = {}
    for r in reco:
        rid = get_reaction_id(r)	
        if not rid in eqns:
            #print "Warning: %s not in balanced equations" % (rid)
            continue
        if rid not in reactions:
            reactions.add(rid)
            rid2re[rid] = r

    mols = set()
    for rid in reactions:
        isbalanced, ostatus, nstatus, eqn = eqns[rid]
        lhs, rhs = eqn
        for mol, qty in lhs:
            mols.add(mol)
        for mol, qty in rhs:
            mols.add(mol)

    for mol in mols:
        if mol in mol2name:
            name = mol2name[mol]
        else:
            name = mol

#      <speciesType metaid="metaid_t_0001" id="t_0001" name="(1-&gt;3)-beta-D-glucan" sboTerm="SBO:0000247">
        st = m.createSpecies()
        metaId = "meta_%s" % (mol)
        st.setId(mol)
        st.setMetaId(metaId)
        st.setSBOTerm("SBO:0000247")
        st.setName("%s [%s]" % (name, DEFAULT_COMPARTMENT))
	st.setInitialConcentration(0)
	st.setInitialAmount(0)
	st.setBoundaryCondition(False)
	st.setHasOnlySubstanceUnits(False)
	st.setConstant(True)
	st.setCompartment(DEFAULT_COMPARTMENT)

        if mol.find("CHEBI") > -1:
            
            chebiNO= mol[5:]
    # CHEBI goes into Bag: <rdf:li rdf:resource="http://identifiers.org/obo.chebi/CHEBI:36464"/>
            st.setAnnotation("""
             <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
               <rdf:Description rdf:about="#%s">
                 <bqbiol:is>
                   <rdf:Bag>
                     <rdf:li rdf:resource="http://identifiers.org/obo.chebi/CHEBI:%s"/> 
                   </rdf:Bag>
                 </bqbiol:is>
               </rdf:Description>
             </rdf:RDF>
    """ % (metaId, chebiNO))
        
	elif mol.find("Cluster") > -1:
    		print ("No link created for metabolite: %s" % mol)
        else:
            st.setAnnotation("""
             <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
               <rdf:Description rdf:about="#%s">
                 <bqbiol:is>
                   <rdf:Bag>
                     <rdf:li rdf:resource="http://identifiers.org/kegg.compound/%s" />
                   </rdf:Bag>
                 </bqbiol:is>
               </rdf:Description>
             </rdf:RDF>
    """ % (metaId, mol))        
        
    for rid in reactions:
        rec = reco[rid2re[rid]]

        isbalanced, ostatus, nstatus, eqn = eqns[rid]
        lhs, rhs = eqn

        create = 0
        if SUM_STOIC_COEFFS == False:
            r = m.createReaction()
            for mol, qty in lhs:
                rea = r.createReactant()
                rea.setSpecies(mol)
                rea.setStoichiometry(qty)
		rea.setConstant(True)
            for mol, qty in rhs:
                rea = r.createProduct()
                rea.setSpecies(mol)
                rea.setStoichiometry(qty)
		rea.setConstant(True)		
            create = 1
        else:

            total = {}
            for mol, qty in lhs:
                if mol not in total:
                    total[mol] = -qty
                else:
                    total[mol] -= qty
            for mol, qty in rhs:
                if mol not in total:
                    total[mol] = qty
                else:
                    total[mol] += qty
            r = None
            for mol in total:
                qty = total[mol]
                if qty < 0:
                    if r is None:
                        r = m.createReaction()
                    rea = r.createReactant()
                    rea.setSpecies(mol)
                    rea.setStoichiometry(-qty)
		    rea.setConstant(True)		
                    create = 1
                elif qty > 0:
                    if r is None:
                        r = m.createReaction()
                    rea = r.createProduct()
                    rea.setSpecies(mol)
                    rea.setStoichiometry(qty)
		    rea.setConstant(True)	
                    create = 1
                else:
                    pass # coeff sum zero -> leave this metabolite out of the model

        if create == 0:
            continue
	
        sbmlrid=rid.replace("-","")
        r.setId(str(sbmlrid))
        parts = rid.split("-")
        keggrid = parts[len(parts)-1]
        if not reKEGGRID.match(keggrid):
            keggrid = None
         
       
        try:
            r_bounds = bounds[rid]
            lb, ub = r_bounds.split(",")
            kinLaw = r.createKineticLaw()
            lbParameter= kinLaw.createLocalParameter()
            lbParameter.setId("LOWER_BOUND")        
            lbParameter.setValue(float(lb))
            ubParameter= kinLaw.createLocalParameter()
            ubParameter.setId("UPPER_BOUND")
            ubParameter.setValue(float(ub))       
        except Exception as ex:
            print "No bounds for %s",  rid
            kinLaw = r.createKineticLaw()
            lbParameter= kinLaw.createLocalParameter()
            lbParameter.setId("LOWER_BOUND")        
            lbParameter.setValue(float(-1000))
            ubParameter= kinLaw.createLocalParameter()
            ubParameter.setId("UPPER_BOUND")
            ubParameter.setValue(float(1000))       


	if not r.isSetId():
            ridfix = rid.replace("-","")
            ridfix = "N"+ridfix
            r.setId(str(ridfix))
            #print ("Adding N to reaction name: %s -> %s" % (rid,ridfix))
        r.setCompartment(DEFAULT_COMPARTMENT)
        r.setReversible(True)
        r.setMetaId("meta_%s" % (rid))	
        r.setSBOTerm("SBO:0000176")
	r.setFast(False)
        if rxnnames is not None:
            if rid in rxnnames:
                rn=rxnnames[rid]
                rn=rn.replace("<i>","")
                rn=rn.replace("</i>","")
                rn=rn.replace("<I>","")
                rn=rn.replace("</I>","")
                rn=rn.replace("<sup>","^")
                rn=rn.replace("</sup>","")
                rn=rn.replace("<sub>","_")
                rn=rn.replace("</sub>","")
                rn=rn.replace("<SUP>","^")
                rn=rn.replace("</SUP>","")
                rn=rn.replace("<SUB>","_")
                rn=rn.replace("</SUB>","")
                rn=rn.replace("<em>","")
                rn=rn.replace("</em>","")
                rn=rn.replace("<small>","_")
                rn=rn.replace("</small>","")
                rn=rn.replace("&alpha;","alpha")
                rn=rn.replace("&beta;","beta")
                rn=rn.replace("&gamma;","gamma")
                rn=rn.replace("&Delta;","Delta")
                rn=rn.replace("&omega;","omega")
                r.setName(rn)
            else:
                #print ("rid: %s not found in rxnnames" % rid) 
                r.setName("")
   
# <reaction metaid="meta_r_0001" sboTerm="SBO:0000176" id="r_0001" name="(R)-lactate:ferricytochrome-c 2-oxidoreductase"> - <notes> - <body xmlns="http://www.w3.org/1999/xhtml">
#  <p>GENE_ASSOCIATION:((YDL174C and YEL039C) or (YDL174C and YJR048W) or (YEL039C and YEL071W) or (YEL071W and YJR048W))</p>
#  </body>
#  </notes>
# - <annotation>
# - <rdf:RDF>
# - <rdf:Description rdf:about="#meta_r_0001">
# - <bqbiol:is>
# - <rdf:Bag>
#  <rdf:li rdf:resource="http://identifiers.org/kegg.reaction/R00197" />  </rdf:Bag>
#  </bqbiol:is>
#  </rdf:Description>
#  </rdf:RDF>
#  </annotation> 


        #kl = libsbml.KineticLaw(SBML_LEVEL, SBML_VERSION)
        #r.setKineticLaw(kl)
        #kl = r.createKineticLaw()

        #ast = libsbml.ASTNode()
        #ast.setName("dummy")
        #kl.setMath(ast)

        #para = kl.createParameter()
        #para.setId("dummy")
        #para.setValue(1)
        #para.setUnits("time")

        annotation = """
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
"""
        if keggrid:
            annotation += """
<rdf:Description rdf:about="#%s">
<bqbiol:is>
<rdf:Bag>
<rdf:li rdf:resource="http://identifiers.org/kegg.reaction/%s" />  
</rdf:Bag>
</bqbiol:is>
""" % ("meta_%s" % (rid), keggrid)
            if pathways is not None:
                if keggrid in pathways:
                    paths = pathways[keggrid].split(",")
                    annotation += """
<bqbiol:isVersionOf>
  <rdf:Bag>
"""
                    for p in paths:
                        annotation += """
     <rdf:li rdf:resource="http://identifiers.org/kegg.pathway/%s"/>
""" % p
                    
                    annotation += """
  </rdf:Bag>
</bqbiol:isVersionOf>
"""
                    if pathwayasnames is not None:
                        pathns = pathwayasnames[keggrid].split(",")
                        note = """
<notes>
  <body xmlns="http://www.w3.org/1999/xhtml">
"""
                        for p in pathns:
                            note += """
     <p>Pathway: %s</p>
""" % p
                        note += """
  </body>
</notes>
"""                     
                        r.setNotes(note)

###
            annotation += """
</rdf:Description>
""" 

        annotation += """
<coreco:ec>%s</coreco:ec>
<coreco:balanced>%s</coreco:balanced>
<coreco:posterior>%s</coreco:posterior>
<coreco:cost>%s</coreco:cost>
<coreco:level>%s</coreco:level>
""" % (rec.ec, isbalanced, rec.p, rec.cost, rec.level)

        if rec.score != None:
            annotation += """
<coreco:naivep>%s</coreco:naivep>
<coreco:btscore>%s</coreco:btscore>
<coreco:gtscore>%s</coreco:gtscore>
<coreco:bscore>%s</coreco:bscore>
<coreco:bseq1>%s</coreco:bseq1>
<coreco:bseq2>%s</coreco:bseq2>
<coreco:gscore>%s</coreco:gscore>
<coreco:gseq1>%s</coreco:gseq1>
<coreco:gseq2>%s</coreco:gseq2>
""" % (rec.score.npscore,
       rec.score.btscore,
       rec.score.gtscore,
       rec.score.bscore,
       rec.score.bseq1,
       rec.score.bseq2,
       rec.score.gscore,
       rec.score.gseq1,
       rec.score.gseq2)

        annotation += "</rdf:RDF>"
        r.setAnnotation(annotation)
        #print r.toSBML()
        #exit()

    return d

def convert_eqn_side(side):
    new = []
    for mol in side:
        vals = mol.split(" ")
        if len(vals) == 2:
            coeff, mol = vals
            try:
                coeff = float(coeff)
            except ValueError:
                coeff = float("NaN")
        else:
            ##mol = re.findall("\w\d{5}", vals[0])[0]
            mol = vals[0].strip()
           # print vals[0] + "......" + mol

            #mol = vals[0]
            coeff = 1
        new.append((mol, coeff))
    return new

def convert_eqn(e):
    lhs, rhs = e.split(" <=> ")
    lhs = lhs.split(" + ")
    rhs = rhs.split(" + ")
    return (convert_eqn_side(lhs), convert_eqn_side(rhs))

def read_balanced_reactions(f):
    reactions = {}
    for s in f:
        if s.startswith("#"):
            continue
        rid, isbalanced, ostatus, nstatus, eqn = s.strip().split("\t")
        if isbalanced == "True":
            eqn = convert_eqn(eqn)
            reactions[rid] = (isbalanced, ostatus, nstatus, eqn)
    return reactions

def main(rdir, eqnfn, molfn, taxon, modelid, modelname, species, outfn, sbmlversion, boundsfile,  rxnnamefile=None, pathwayfile=None):
    print("Loading molecule names... ")
    # dictionary where keys are mol ids (C00001) and
    # items are names from second column of kegg-compounds file
    mol2name = common.parse_molecule_names(molfn)
    print "Loading reconstruction: %s/%s" % (rdir, common.NETWORK_REACTION_FILE)
    f = open("%s/%s" % (rdir, common.NETWORK_REACTION_FILE))
    bf = open(eqnfn)
    reco = common.read_reconstruction(f)
    eqns = read_balanced_reactions(bf)
    if pathwayfile is not None:
        fp = open(pathwayfile)
        pathways={};
        pathwayasnames={};
        for s in fp: 
            #print s
            sisalto = s.strip().split("\t") 
            #print len(sisalto)
            if len(sisalto)>1:
                pathways[sisalto[0]]=sisalto[1]
            if len(sisalto)>2:
                pathwayasnames[sisalto[0]]=sisalto[2]
        if len(pathwayasnames) == 0:
            pathwayasnames = None
    else:
        pathways = None
    if rxnnamefile is not None:
        fp = open(rxnnamefile)
        rxnnames={};
        for s in fp: 
            sisalto = s.strip().split("\t") 
            if len(sisalto)>1:
                rxnnames[sisalto[0]]=sisalto[1]
                #print ("rxnnames[%s] = %s" % (sisalto[0],sisalto[1]))                
    else:
        rxnnames = None
        
    bounds={};
    if boundsfile is not None:
        fp = open(boundsfile)        
        for s in fp: 
            sisalto = s.strip().split("\t") 
            if len(sisalto)>1:
                bounds[sisalto[0].strip()]=sisalto[3] +","+ sisalto[4]
  #  print "%d reactions" % (len(reco))
    sbml = convert_to_SBML(reco, eqns, mol2name, taxon, modelid, modelname, species, sbmlversion, bounds,  rxnnames, pathways, pathwayasnames)
    libsbml.writeSBMLToFile(sbml, outfn)

if __name__ == "__main__":
    rdir = sys.argv[1]  # dir with network.reactions
    eqnfn = sys.argv[2] # balanced reactions, "balances.eqn"
    molfn = sys.argv[3]  # metabolite id -> name
    taxon = sys.argv[4] # NCBI taxonomy id
    modelid = sys.argv[5]
    modelname = sys.argv[6]
    species = sys.argv[7]
    outfn = sys.argv[8] # sbml output
    sbmlversion = sys.argv[9] # sbml version (2 or 3)
    boundsfile = sys.argv[10] # sbml version (2 or 3)
    pathwayfile = None
    rxnnamefile = None
    if len(sys.argv) > 11:
        if(os.path.isfile(sys.argv[11])):
            #print "Reading reaction paths"
            rxnnamefile=sys.argv[11]
    if len(sys.argv) > 12:
        if(os.path.isfile(sys.argv[12])):
            #print "Reading reaction paths"
            pathwayfile=sys.argv[12]

    main(rdir, eqnfn, molfn, taxon, modelid, modelname, species, outfn, sbmlversion, boundsfile,  rxnnamefile, pathwayfile)

