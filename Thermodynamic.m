%	Politecnico di Milano (2022)
%
%	Bachelor's degree in chemical engineering
%
%	This code was developed and tested by Elia Ferretti
%
%	You can redistribute the code and/or modify it
%	Whenever the code is used to produce any publication or document,
%	reference to this work and author should be reported
%	No warranty of fitness for a particular purpose is offered
%	The user must assume the entire risk of using this code
%
%------------------------------------------------------------------------------------------------------------




function y = cp_perry(T_r,cp_param)
%The function returns the value of cp @ a given temperature T according to
%the hyperbolic function

    [row,column] = size(cp_param);
    if row>1 && column>1
        %return a vector
        a = cp_param(:,1);
        b = cp_param(:,2);
        c = cp_param(:,3);
        d = cp_param(:,4);
        e = cp_param(:,5);
    else
        %return a single value
        a = cp_param(1);
        b = cp_param(2);
        c = cp_param(3);
        d = cp_param(4);
        e = cp_param(5);
    end
        
    perry = @(T) a + b.*(c./T./sinh(c./T)).^2 + d.*(e./T./cosh(e./T)).^2;
    y = perry(T_r);
end

function y = cpdT_perry(T0,Tf,cp_param)
%The function returns the value of the analitic integral of cp in dT with T from
%T0 to Tf, "cp_param" can be a matrix (integral for each species are evaluated @
%the same function call) or "cp_param" can be a vector (integral for the specifiec
%species are evaluated @ function call)

    [row,column] = size(cp_param);
    if row>1 && column>1
        %return a vector
        a = cp_param(:,1);
        b = cp_param(:,2);
        c = cp_param(:,3);
        d = cp_param(:,4);
        e = cp_param(:,5);
    else
        %return a single value
        a = cp_param(1);
        b = cp_param(2);
        c = cp_param(3);
        d = cp_param(4);
        e = cp_param(5);
    end

    
    perry = @(T) a*T+b.*c.*coth(c/T)-d.*e.*tanh(e/T);

    y = perry(Tf)-perry(T0);
end

function y = reactionHentalpy(v,T,dh0f_perry,cp_param,Tr)
%The function returns the reaction hentalpy given the temperature and the
%stoichiometric coefficients of the reaction using cpdT_perry for cp
%integration

    y = v'*dh0f_perry;
    for i=1:length(v)
        y = y + v(i)*cpdT_perry(Tr,T,cp_param(i,:));
    end
end

function y = dG0R_VantHoff(dg0f,dh0f,cp_param,T,Tr,v)
%The function returns the value of Keq = exp(-dG0R(T)/R/T) known all the input variables
%using the analitically integrated formulation of the Van't Hoff equation

    R = 8.31446261815324;   %[J/mol/K]
    NC = length(v);
    dG0R_rif = dg0f'*v;
    alfa = dG0R_rif/R/Tr;
    beta = 0;
    
    for i=1:NC
        costante = dh0f(i)-(cp_param(i,1)*Tr+cp_param(i,2)*cp_param(i,3)*coth(cp_param(i,3)/Tr)-cp_param(i,4)*cp_param(i,5)*tanh(cp_param(i,5)/Tr));
        beta = beta + v(i)*( costante*(1/Tr-1/T) + cp_param(i,1)*log(T/Tr)...
            - cp_param(i,2)*log(tanh(cp_param(i,3)/T)*cosh(cp_param(i,3)/T)/(tanh(cp_param(i,3)/Tr)*cosh(cp_param(i,3)/Tr)))...
            + cp_param(i,4)*log(cosh(cp_param(i,5)/T)/cosh(cp_param(i,5)/Tr)) );
    end
    beta = beta/R;
    adimensionalGibbsOfReaction = alfa-beta;
    y = exp(-adimensionalGibbsOfReaction);
end

function y = viscosity_perry(mu_param,T)
%The function returns the value of dynamic viscosity @ a given temperature T according to
%the hyperbolic function (Perry)

    y = mu_param(:,1).*T.^mu_param(:,2)./(1+mu_param(:,3)./T+mu_param(:,4)./T.^2);
end

function y = thermalConductivity_perry(k_param,T)
%The function returns the value of thermal conductivity @ a given temperature T according to
%the hyperbolic function (Perry)

    y = k_param(:,1).*T.^k_param(:,2)./(1+k_param(:,3)./T+k_param(:,4)./T.^2);
end

function y = mixtureThermalConductivity(x,T,k_param)
%The function returns the value of thermal conductivity of a mixture @ a given temperature T
%given its composition x. The function need the function "thermalConductivity_perry".
%The mixing rule (3) of the following article is used.
%https://www.researchgate.net/publication/311803060_On_Thermal_Conductivity_of_Gas_Mixtures_Containing_Hydrogen

    sum1 = 0;
    sum2 = 0;
    k = thermalConductivity_perry(k_param,T);
    for i=1:length(x)
        sum1 = sum1 + x(i)*k(i);
        sum2 = sum2 + x(i)/k(i);
    end
    y = 0.5*(sum1+1/sum2);
end

function y = viscosityCorrectionWithPressure(T,p,NC,Tc,pc,w,rho_c,mu_param,PM)
%The function returns the value of dynamic viscosity @ a given temperature T and pressure according to
%the function "viscosity_perry" and correlation from Perry's text book

    %p [Pa], T[K], PM [kg/mol]
	
    y = zeros(NC,1);
    rho = zeros(NC,1);
    R = 8.3144621;              %[J/mol/K]

    psi = 2173.4.*(Tc).^(1/6).*(PM*1000).^(-0.5).*(pc.*1e-6).^(-2/3);
    

    for i=1:NC
		%Comprimibility factor must be defined with some EoS
		%for pure copmponent (no mixing rules needed)
		zeta = zetaPR(i,T,p,'V');
        rho(i) = p/R/T/zeta;
    end
	
    rho_r = rho./rho_c;
    mu0 = viscosity_perry(mu_param,T);

    for i=1:NC
        if i==2 || i==3
			%Correlation for POLAR gasses
			%In this case compounds 2 and 3 are POLAR
            if rho_r(i)<=0.1
                %case 1                
                y(i) = (1.656*(rho_r(i))^1.111)/psi(i) + mu0(i);
            elseif rho_r(i)>0.1 && rho_r(i)<=0.9
                %case 2
                y(i) = (0.0607*(9.045*rho_r(i)+0.63)^1.739)/psi(i) + mu0(i);
            elseif rho_r(i)>0.9 && rho_r(i)<=2.2
                %case 3
                polinomiale = 0.6439-0.1005*rho_r(i);
                y(i) = exp(4-exp(polinomiale))/psi(i) + mu0(i);
            elseif rho_r(i)>2.2
                %case 4
                polinomiale = 0.6439-0.1005*rho_r(i)-0.000475*(rho_r(i)^3-10.65)^2;
                y(i) = exp(4-exp(polinomiale))/psi(i) + mu0(i);
            end
        else
			%Correlation for NON-POLAR gasses
            polinomiale = 1.0230+0.23364*rho_r(i)+0.58533*rho_r(i)^2-0.40758*rho_r(i)^3+0.093324*rho_r(i)^4;
            y(i) = (polinomiale^(4)-1)/psi(i) + mu0(i);
        end 
    end
end

function y = mixtureViscosity(T,x,mu_param,PM)
%The function returns the value of dynamic viscosity of a mixture @ a given temperature T
%given its composition x. The function needs the function "viscosity_perry" to compute
%components dynamic viscosity.
%The mixing rule (3) of the following article is used.
%https://stacks.cdc.gov/view/cdc/10045/cdc_10045_DS1.pdf
	
	NC = length(x);
    phi = zeros(NC,NC);
	mu = viscosity_perry(mu_param,T);
	
    for i=1:NC
        for j=1:NC
            phi(i,j) = ( 1+ (mu(i)/mu(j))^0.5*(PM(j)/PM(i))^0.25 )^2/( 4/sqrt(2)*(1+PM(i)/PM(j))^0.5 );
        end
    end
    y = 0;
    somma = 0;
    for i=1:NC
        for j=1:NC
            if not(j==i)
                somma = somma + x(j)*phi(i,j);
            end
        end
        y = y + x(i)*mu(i)/(x(i)+somma);
    end
end

function y = NcpdT(f,T0,Tf,cp_param)
%The function returns the value of the sum of the product between the integral of
%cp from T0 to Tf and the molar flow f(i) for every species
%The function needs the function "cpdT_perry"

	y = 0;
	for i=1:length(f)
		y = y + f(i)*cpdT_perry(T0,Tf,cp_param(i,:));
	end
end