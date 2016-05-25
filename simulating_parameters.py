############# PARAMETERS ###################

def giveParametersDeclaration ( model ) :
	param_names = []
	params_def = "struct Params\n{\n"
	event_types = [ "det_events" , "stoch_events" , "ode_dynamics" ]
	for event_type in event_types :
		for event in eval ( "model." + event_type ) :
			for param in event.parameters.keys () :
				if param not in param_names :
					param_names.append ( param )
					params_def += "\tstatic Doub " + param + " ; \n"
	params_def += "} ;\n\n"
	return params_def

def giveParametersDefinition ( model ) :
	param_names = []
	params_decl = ""
	event_types = [ "det_events" , "stoch_events" , "ode_dynamics" ]
	for event_type in event_types :
		for event in eval ( "model." + event_type ) :
			for param in event.parameters.keys () :
				if param not in param_names :
					param_names.append ( param )
					params_decl += "Doub Params::" + param + " = " + str ( event.parameters[param] ) + " ; \n"
	params_decl += "\n\n"					
	return params_decl



def writeParameterCode ( model , target_folder ) :
	writer = open ( target_folder + "/parameters.hpp" , "w" )
	writer.write ( "\n#ifndef PARAMETERS_HPP\n#define PARAMETERS_HPP\n\n" )
	writer.write ( "\n#include \"nr3.hpp\"\n\n" )
	writer.write ( giveParametersDeclaration ( model ) )
	writer.write ( "#endif\n" )
	writer = open ( target_folder + "/parameters.cpp" , "w" )
	writer.write ( "\n#include \"parameters.hpp\"\n\n" )
	writer.write ( giveParametersDefinition ( model ) )	


