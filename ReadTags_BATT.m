function [Data] = ReadTags_BATT(NameFile)

DIM.S ={'Water', 0.00;
        'Sac1', 0.001;
        'Sac2', 0.002;
        'Sac3', 0.003;
        'Sac4', 0.004;
        'Sac5', 0.005};
    
[Data] = MedPC_Tags(NameFile, DIM);

end %Function