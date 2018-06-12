function irf = SIME_getIRF(K,t_int,model,Vnd)

    switch model
        case '1TCM'
            irf = K(1)*exp(-K(2)*t_int);
            
        case '2TCM'
            s   = K(2) + K(3) + K(4);
            a1  = .5* (s - sqrt(s^2 - 4*K(2)*K(4)) );
            a2  = .5* (s + sqrt(s^2 - 4*K(2)*K(4)) );
            p   = (K(3) + K(4) - a1) / (a2 - a1);
            irf = K(1) * (p*exp(-t_int*a1) + (1-p)*exp(-t_int*a2));
            

        case '2TCMfixedVnd'  % nota bene: [K(1) K(2) K(3)] = [k2 k3 k4]
            s   = K(1) + K(2) + K(3);
            a1  = .5* (s - sqrt(s^2 - 4*K(1)*K(3)) );
            a2  = .5* (s + sqrt(s^2 - 4*K(1)*K(3)) );
            p   = (K(2) + K(3) - a1) / (a2 - a1);
            irf = Vnd*K(1) * (p*exp(-t_int*a1) + (1-p)*exp(-t_int*a2));                                   
    end 