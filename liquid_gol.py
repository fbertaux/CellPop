#!/usr/bin/python

import modeling
import simulating


### model construction

model = modeling.Model ("liquid_gol")

# agents
model.addAgent ( name="Cell" , unique=False , properties=[ "age" , "CC_length" ] )
model.addAgent ( name="Medium" , unique=True , properties=[ "IP" ] )

# event: cell division
model.addDeterministicEvent ( name="cell_division" , 
							  agent="Cell" ,
							  kind="creation" , 
							  parameters={"CC_avg":3.0,"CC_std":0.25} ,
							  trigger_expression="age > CC_length" ,
							  realization_function=
							  		["new.age = age - CC_length" ,
							  		 "age = new.age" ,
							  		 "new.CC_length = ran.norm ( CC_avg , CC_std ) " ,
							  		 "CC_length = ran.norm ( CC_avg , CC_std ) " ] )

# event: cell death
model.addStochasticEvent ( name="cell_death" ,
						   agent="Cell" , 
						   kind="destruction" , 
						   parameters={ "IP_threshold":1.0 , "IP_death_rate_slope":0.01 } ,
						   propensity_expression="IP > IP_threshold ? IP * IP_death_rate_slope : 0." )

# continuous dynamics
model.addContinuousChange ( name="cell_aging" , agent="Cell" , property="age" , rate_expression="1." )

model.addContinuousChange ( name="IP_degradation" , agent="Medium" , property="IP" ,
							parameters={ "IP_degradation_rate":10. } ,
							rate_expression = "- IP * IP_degradation_rate" )

model.addContinuousChange ( name="IP_production" , agent="Medium" , property="IP" , agent_source="Cell" , 
							parameters={ "IP_prod_rate":0.01 } ,
							rate_expression="IP_prod_rate" )


### model simulation (i.e. generation of cpp code for simulation)

simulating.writeAllCode ( model=model , target_folder="liquid_gol" )
