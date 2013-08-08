#!/usr/bin/env python

import sys

#sys.path.append("/usr/local/lib/")

import libsbml

sys.path.append("../model_training_scripts/")
import common

SUM_STOIC_COEFFS = 1  # sum coefficients of both sides: 2A + 2B <=> 3A + 1C ---> 2B <=> 1A + 1C

SBML_LEVEL = 2
SBML_VERSION = 4

NOTES = """
<body xmlns="http://www.w3.org/1999/xhtml">
<p>
This model was reconstructed with CoReCo method from protein sequence and phylogeny data. CoReCo is described in the manuscript "Comparative genome-scale reconstruction of gapless metabolic networks for present and ancestral species" (submitted) by Esa Pitkanen, Paula Jouhten, Peter Blomberg, Liisa Holm, Merja Penttila, Juho Rousu and Mikko Arvas.

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
    return r.split("_")[0]

def convert_to_SBML(reco, eqns, mol2name, taxon, modelid, modelname, species):

    DEFAULT_COMPARTMENT = "cytosol"

    d = libsbml.SBMLDocument(SBML_LEVEL, SBML_VERSION)

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
    comp.setConstant(True) # Sets the value of the 'constant' attribute of this Compartment.
    comp.setName(DEFAULT_COMPARTMENT)
    #comp.setOutside() # sid the identifier of a compartment that encloses this one.
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
            print "Warning: %s not in balanced equations" % (rid)
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
        st = m.createSpeciesType()
        metaId = "meta_%s" % (mol)
        st.setId(metaId)
        st.setMetaId(metaId)
        st.setName(name)
        st.appendNotes("SBO:0000247")
# CHEBI goes into Bag: <rdf:li rdf:resource="http://identifiers.org/obo.chebi/CHEBI:36464"/>
        st.setAnnotation("""
         <rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">
           <rdf:Description rdf:about="#%s">
             <bqbiol:is>
               <rdf:Bag>
                 <rdf:li rdf:resource="http://identifiers.org/obo.chebi/CHEBI:37671" />  <rdf:li rdf:resource="http://identifiers.org/kegg.compound/%s" />
               </rdf:Bag>
             </bqbiol:is>
           </rdf:Description>
         </rdf:RDF>
""" % (metaId, mol))

# <speciesType metaid="meta_t_0001" sboTerm="SBO:0000247" id="t_0001" name="(1->3)-beta-D-glucan"> - <annotation> - <rdf:RDF> - <rdf:Description rdf:about="#meta_t_0001"> - <bqbiol:is> - <rdf:Bag>
#  <rdf:li rdf:resource="http://identifiers.org/obo.chebi/CHEBI:37671" />  <rdf:li rdf:resource="http://identifiers.org/kegg.compound/C00965" />  </rdf:Bag>
#  </bqbiol:is>
#  </rdf:Description>
#  </rdf:RDF>
#  </annotation>
#  </speciesType>


#      <species id="s_0001" name="(1-&gt;3)-beta-D-glucan [cell envelope]" speciesType="t_0001" compartment="c_01"/>
        sp = m.createSpecies()
        sp.setCompartment(DEFAULT_COMPARTMENT)
        sp.setId(mol)
        sp.setSpeciesType(metaId)
        sp.setName("%s [%s]" % (name, DEFAULT_COMPARTMENT))
        
    for rid in reactions:
        rec = reco[rid2re[rid]]

        isbalanced, ostatus, nstatus, eqn = eqns[rid]
        lhs, rhs = eqn

        create = 0
        if SUM_STOIC_COEFFS == False:
            r = m.createReaction()
            for mol, qty in lhs:
                re = r.createReactant()
                re.setSpecies(mol)
                re.setStoichiometry(qty)
            for mol, qty in rhs:
                re = r.createProduct()
                re.setSpecies(mol)
                re.setStoichiometry(qty)
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
                    if r == None:
                        r = m.createReaction()
                    re = r.createReactant()
                    re.setSpecies(mol)
                    re.setStoichiometry(-qty)
                    create = 1
                elif qty > 0:
                    if r == None:
                        r = m.createReaction()
                    re = r.createProduct()
                    re.setSpecies(mol)
                    re.setStoichiometry(qty)
                    create = 1
                else:
                    pass # coeff sum zero -> leave this metabolite out of the model

        if create == 0:
            continue

        r.setId(rid)
        r.setCompartment(DEFAULT_COMPARTMENT)
        r.setReversible(True)
        r.setMetaId("meta_%s" % (rid))
        r.appendNotes("SBO:0000176")

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
        kl = r.createKineticLaw()

        ast = libsbml.ASTNode()
        ast.setName("dummy")
        kl.setMath(ast)

        para = kl.createParameter()
        para.setId("dummy")
        para.setValue(1)
        para.setUnits("time")

        annotation = """
<rdf:RDF>
<rdf:Description rdf:about="#%s">
<bqbiol:is>
<rdf:Bag>
<rdf:li rdf:resource="http://identifiers.org/kegg.reaction/%s" />  
</rdf:Bag>
</bqbiol:is>
</rdf:Description>
""" % ("meta_%s" % (rid), rid)

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

import re

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
            mol = re.findall("\w\d{5}", vals[0])[0]
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

def main(rdir, eqnfn, molfn, taxon, modelid, modelname, species, outfn):
    mol2name = {}
    f = open(molfn)
    for s in f:
        molid, name, name2 = s.strip().split("\t")
        mol2name[molid] = name
    print "Loading reconstruction: %s/%s" % (rdir, common.NETWORK_REACTION_FILE)
    f = open("%s/%s" % (rdir, common.NETWORK_REACTION_FILE))
    bf = open(eqnfn)
    reco = common.read_reconstruction(f)
    eqns = read_balanced_reactions(bf)
    print "%d reactions" % (len(reco))
    sbml = convert_to_SBML(reco, eqns, mol2name, taxon, modelid, modelname, species)
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
    main(rdir, eqnfn, molfn, taxon, modelid, modelname, species, outfn)
