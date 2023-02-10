function [L1, L2, L3, L4] = bh_get_preconfigured_machine(act_str)

act_str = upper(string(act_str));

switch(act_str)
    case "MACHINE_A"
            L1 = 0.4;
            L2 = 1.4;
            L3 = 1.2;
            L4 = 1;        
    case "MACHINE_B"
            L1 = 1.2;
            L2 = 1.4;
            L3 = 0.4;
            L4 = 1;
    case "MACHINE_C"
            L1 = 0.3;
            L2 = 1.4;
            L3 = 0.4;
            L4 = 1;        
    case "MACHINE_D"
            L1 = 0.39;
            L2 = 1;
            L3 = 0.4;
            L4 = 1;       
    case "MACHINE_D2"
            L1 = 0.4;
            L2 = 1;
            L3 = 0.4;
            L4 = 1;                   
    case "MACHINE_E"
            L1 = 0.5;
            L2 = 0.5;
            L3 = 0.5;
            L4 = 1;        
    case "MACHINE_F"
            L1 = 0.4;
            L2 = 0.5;
            L3 = 1;
            L4 = 1;  
    case "MACHINE_G"
           L1        = 1;   % m
           L2        = 1.4;   % m
           L3        = 1.2;   % m
           L4        = 0.5;     % m     
    case "MACHINE_X"
           L1        = 0.4;   % m
           L2        = 0.7;   % m
           L3        = 0.8;   % m
           L4        = 1.0;     % m             
    otherwise
        error("###_ERROR:  UNknown machine type !")
end

