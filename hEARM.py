#!/usr/bin/python


###########################################################################
# Description of the 'hEARM' model of TRAIL-induced apoptosis, proposed in
# Bertaux et al., PLoS Computational Biology, 2014.
#
# It extends the EARM model used in Spencer et al., Nature, 2009,
# by accounting for protein fluctuations.
###########################################################################

# cellpop module
import cellpop

### define native (i.e. fluctuating) proteins
def defineNativeProteins (model) :

	# direct path to death
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="R" , EP=1000. , dilutionHalfLife=27. )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="pC8" , EP=10000. , dilutionHalfLife=27. )
	model.addFluctuatingProtein ( cell_type="HeLa" , name="pC3" , EP=10000. , HLP=27. , EM=17. , HLM=9. , Ton=0.099 , Toff=3.46 )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="pC6" , EP=10000. , dilutionHalfLife=27. )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="PARP" , EP=1000000. , dilutionHalfLife=27. )

	# direct path to death : inhibitors
	model.addFluctuatingProtein ( cell_type="HeLa" , name="flip" , EP=2000. , HLP=0.4 , EM=17. , HLM=1. , Ton=16. , Toff=24. )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="Bar" , EP=1000. , dilutionHalfLife=27. )
	model.addFluctuatingProtein ( cell_type="HeLa" , name="XIAP" , EP=100000. , HLP=27. , EM=17. , HLM=9. , Ton=0.098 , Toff=3.64 )

	# indirect path to death : mitchondrial pathway
	model.addFluctuatingProtein ( cell_type="HeLa" , name="Bid" , EP=60000. , HLP=27. , EM=17. , HLM=9. , Ton=0.098 , Toff=3.64 )
	model.addFluctuatingProtein ( cell_type="HeLa" , name="Bax" , EP=80000. , HLP=27. , EM=17. , HLM=9. , Ton=0.099 , Toff=3.14 )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="Pore" , EP=500000. , dilutionHalfLife=27. )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="smac_m" , EP=100000. , dilutionHalfLife=27. )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="CyC_m" , EP=500000. , dilutionHalfLife=27. )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="pC9" , EP=100000. , dilutionHalfLife=27. )
	model.addFluctuatingProtein ( cell_type="HeLa" , standard=True , name="Apaf" , EP=100000. , dilutionHalfLife=27. )

	# indirect path to death : inhibitors
	model.addFluctuatingProtein ( cell_type="HeLa" , name="Mcl1" , EP=20000. , HLP=0.4 , EM=17. , HLM=1. , Ton=16. , Toff=24.  )
	model.addFluctuatingProtein ( cell_type="HeLa" , name="Bcl2" , EP=30000. , HLP=27. , EM=17. , HLM=9. , Ton=0.098 , Toff=3.73 )


### Signaling reactions : receptor activation and caspase cascade
def defineCaspaseCascadeReactions (model,ref_kb,ref_ku,ref_kc) :
	model.addCellularReaction (cell_type="HeLa",name="TrailBindsReceptor",reactants=["Env.TRAIL","R"],products=["Env.TRAIL","TRAIL_R"],rate=("kb",4.*ref_kb))
	model.addCellularReaction (cell_type="HeLa",name="TrailUnbindsReceptor",reactants=["TRAIL_R"],products=["R"],rate=("ku",1.e-3*ref_ku))
	model.addCellularReaction (cell_type="HeLa",name="ReceptorActivation",reactants=["TRAIL_R"], products=["R_act"],rate=("kc",0.01*ref_kc))
	model.addCellularReaction (cell_type="HeLa",name="ActiveReceptorGivesTrailBack",reactants=["R_act"], products=["R"],rate=("k",10.*ref_ku))
	model.addCellularCatalyticReaction (cell_type="HeLa",name="Caspase8Activation",substrate="pC8",catalyst="R_act",product="C8",rates=[("kb",ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="Caspase3Activation",substrate="pC3",catalyst="C8",product="C3",rates=[("kb",ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="Caspase6Activation",substrate="pC6",catalyst="C3",product="C6",rates=[("kb",0.),("ku",0.),("kc",0.)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="Caspase8ActivationByC6",substrate="pC8",catalyst="C6",product="C8",rates=[("kb",0.),("ku",0.),("kc",0.)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="PARPCleavage",substrate="PARP",catalyst="C3",product="cPARP",rates=[("kb",10.*ref_kb),("ku",ref_ku),("kc",20.*ref_kc)])

### Signaling reactions : inhibitors of caspase cascade
def defineCaspaseCascadeInhibitionReactions (model,ref_kb,ref_ku,ref_kc) :
	model.addCellularReversibleReaction (cell_type="HeLa",name="flipBindsActivatedReceptor",reactants=["flip","R_act"],products=["flip_R_act"],rates=[("kb",10.*ref_kb),("ku",ref_ku)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="BarBindsCaspase8",reactants=["Bar","C8"],products=["Bar_C8"],rates=[("kb",10.*ref_kb),("ku",ref_ku)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="XIAPdegradesCaspase3",substrate="C3",catalyst="XIAP",product="C3_ub",rates=[("kb",20.*ref_kb),("ku",ref_ku),("kc",0.1*ref_kc)])

### Signaling reactions : from Bid cleavage to mitochondrial pore formation
def defineUpstreamMompReactions (model,ref_kb,ref_ku,ref_kc,mito_cyto_vol_ratio) :
	model.addCellularCatalyticReaction (cell_type="HeLa",name="BidCleavage",substrate="Bid",catalyst="C8",product="tBid",rates=[("kb",ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="BaxActivation",substrate="Bax",catalyst="tBid",product="Bax_act",rates=[("kb",ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="BaxTranslocation",reactants=["Bax_act"],products=["Bax_act_m"],rates=[("ktransf",3600.*0.01),("ktransb",3600.*1.)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="BaxDimerization",reactants=["Bax_act_m","Bax_act_m"],products=["Bax_act_m_2"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="BaxTetramerization",reactants=["Bax_act_m_2","Bax_act_m_2"],products=["Bax_act_m_4"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="BaxTetramerBindsPore",reactants=["Pore","Bax_act_m_4"],products=["Pore_Bax_act_m_4"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addCellularReaction (cell_type="HeLa",name="PoreFormation",reactants=["Pore_Bax_act_m_4"],products=["Pore_act"],rate=("kc",ref_kc))


### Signaling reactions : inhibition of MOMP
def defineUpstreamMompInhibitionReactions (model,ref_kb,ref_ku,mito_cyto_vol_ratio) :
	model.addCellularReversibleReaction (cell_type="HeLa",name="Mcl1BindstBid",reactants=["Mcl1","tBid"],products=["Mcl1_tBid"],rates=[("kb",10.*ref_kb),("ku",ref_ku)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="Bcl2BindsBax",reactants=["Bcl2","Bax_act_m"],products=["Bcl2_Bax_act_m"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="Bcl2BindsBax2",reactants=["Bcl2","Bax_act_m_2"],products=["Bcl2_Bax_act_m_2"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="Bcl2BindsBax4",reactants=["Bcl2","Bax_act_m_4"],products=["Bcl2_Bax_act_m_4"],rates=[("kb",10.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="XIAPBindsApop",reactants=["XIAP","Apop"],products=["XIAP_Apop"],rates=[("kb",20.*ref_kb),("ku",ref_ku)])

### Signaling reactions : MOMP reactions
def defineMompReactions (model,ref_kb,ref_ku,ref_kc,mito_cyto_vol_ratio) :
	model.addCellularCatalyticReaction (cell_type="HeLa",name="CyCRelease1",substrate="CyC_m",catalyst="Pore_act",product="CyC_r",rates=[("kb",20.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku),("kc",10.*ref_kc)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="SmacRelease1",substrate="smac_m",catalyst="Pore_act",product="smac_r",rates=[("kb",20.*ref_kb/mito_cyto_vol_ratio),("ku",ref_ku),("kc",10.*ref_kc)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="CyCRelease2",reactants=["CyC_r"],products=["CyC"],rates=[("ktransf",3600.*1.),("ktransb",3600.*0.01)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="SmacRelease2",reactants=["smac_r"],products=["smac"],rates=[("ktransf",3600.*1.),("ktransb",3600.*0.01)])

### Signaling reactions : Downstream MOMP
def defineDownstreamMompReactions (model,ref_kb,ref_ku,ref_kc) :
	model.addCellularReversibleReaction (cell_type="HeLa",name="SmacBindsXIAP",reactants=["smac","XIAP"],products=["smac_XIAP"],rates=[("kb",70.*ref_kb),("ku",ref_ku)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="ApafActivation",substrate="Apaf",catalyst="CyC",product="Apaf_act",rates=[("kb",5.*ref_kb),("ku",ref_ku),("kc",ref_kc)])
	model.addCellularReversibleReaction (cell_type="HeLa",name="ApoptosomeFormation",reactants=["Apaf_act","pC9"],products=["Apop"],rates=[("kb",0.5*ref_kb),("ku",ref_ku)])
	model.addCellularCatalyticReaction (cell_type="HeLa",name="ApopCleavesCaspase3",substrate="pC3",catalyst="Apop",product="C3",rates=[("kb",0.05*ref_kb),("ku",ref_ku),("kc",ref_kc)])



### set the degradation of non native protein species
def setDegradationRates (model) :
	model.setAllModifiedProteinDegradation (cell_type="HeLa" , halfLife=5.)
	model.setModifiedProteinDegradation (cell_type="HeLa" , name="Pore_act" , halfLife=1.9)
	model.setModifiedProteinDegradation (cell_type="HeLa" , name="flip_R_act" , halfLife=0.4)
	model.setModifiedProteinDegradation (cell_type="HeLa" , name="Mcl1_tBid" , halfLife=0.4)
	model.setModifiedProteinDegradation (cell_type="HeLa" , name="cPARP" , halfLife=27)
	model.setEnvironmentSpeciesDegradation (name="TRAIL",halfLife=9.)


### MAIN

## the model object
model = cellpop.Model ("hEARM")

# create HeLa cell type
model.addCellType ("HeLa")

## define fluctuation models
defineNativeProteins (model=model)

## references rates parameters
ref_kb = 1.e-7*3600.
ref_ku = 0.001*3600.
ref_kc = 1.*3600.
mito_cyto_vol_ratio = 0.07

## caspase cascade
defineCaspaseCascadeReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc)
defineCaspaseCascadeInhibitionReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc)

## Upstream MOMP
defineUpstreamMompReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc,mito_cyto_vol_ratio=mito_cyto_vol_ratio)
defineUpstreamMompInhibitionReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,mito_cyto_vol_ratio=mito_cyto_vol_ratio)

## MOMP and downstream MOMP
defineMompReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc,mito_cyto_vol_ratio=mito_cyto_vol_ratio)
defineDownstreamMompReactions (model=model,ref_kb=ref_kb,ref_ku=ref_ku,ref_kc=ref_kc)

## degradation of modified forms
setDegradationRates (model=model)

## define cell death rule
model.setCellDeathRule (cell_type="HeLa",params={"death_cPARP_threshold":100000.},death_bool_expression="cPARP > death_cPARP_threshold")

## define cell division rule
model.defineDivision (cell_type="HeLa",cell_cycle_length_avg=27.,cell_cycle_length_std=3.)

# test translation
transmodel = cellpop.translateCellPopModel (model)
# transmodel.display ()
cellpop.simulating.writeAllCode ( model=transmodel , target_folder="hEARM" )



