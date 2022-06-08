# Thermodynamic
Miscellaneous functions to compute thermodynamic properties of mixture and pure component

1) cp_perry: 	
	The function returns the value of cp @ a given temperature T according to
	the hyperbolic function
				
2) cpdT_perry:
	The function returns the value of the analitic integral of cp in dT with T from
	T0 to Tf, "cp_param" can be a matrix (integral for each species are evaluated @
	the same function call) or "cp_param" can be a vector (integral for the specifiec
	species are evaluated @ function call)				
	
3) reactionHentalpy:
	The function returns the reaction hentalpy given the temperature and the
	stoichiometric coefficients of the reaction using cpdT_perry for cp
	integration
	
4) dG0R_VantHoff:
	The function returns the value of Keq = exp(-dG0R(T)/R/T) known all the input
	variables using the analitically integrated formulation of the Van't Hoff equation
	
5) viscosity_perry:
	The function returns the value of dynamic viscosity @ a given temperature T
	according to the hyperbolic function (Perry's)
	
6) thermalConductivity_perry:
	The function returns the value of thermal conductivity @ a given temperature T
	according to the hyperbolic function (Perry's)
	
7) mixtureThermalConductivity:
	The function returns the value of thermal conductivity of a mixture @ a given temperature T
	given its composition x. The function need the function "thermalConductivity_perry".
	The mixing rule (3) of the following article is used.
	https://www.researchgate.net/publication/311803060_On_Thermal_Conductivity_of_Gas_Mixtures_Containing_Hydrogen

8) viscosityCorrectionWithPressure:
	The function returns the value of dynamic viscosity @ a given temperature T and
	pressure according to the function "viscosity_perry" and correlation from Perry's text book
	
9) mixtureViscosity:
	The function returns the value of dynamic viscosity of a mixture @ a given temperature T
	given its composition x. The function needs the function "viscosity_perry" to compute
	components dynamic viscosity.
	The mixing rule (3) of the following article is used.
	https://stacks.cdc.gov/view/cdc/10045/cdc_10045_DS1.pdf
	
10) NcpdT:
	The function returns the value of the sum of the product between the integral of
	cp from T0 to Tf and the molar flow f(i) for every species
	The function needs the function "cpdT_perry"
	
Please report any problem or bug to the developer at: eliaferretti@outlook.it
