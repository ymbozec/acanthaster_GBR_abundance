%
% BOOTSTRAP RE-SAMPLING OF MANTA TOW COUNTS FOR CONFIDENCE INTERVAL CALCULATIONS
%
% -----------------------------------------------------------------------------
% Yves-Marie Bozec, The University of Queensland, Australia (y.bozec@uq.edu.au)
% March 2023
% -----------------------------------------------------------------------------

clear

load('LTMP_COTS_MTC.mat')
% Reef-level manta count of CoTS from AIMS Long Term Monitoring Program
% 
% Australian Institute of Marine Science (AIMS). (2015). AIMS Long-term Monitoring Program:
% Crown-of-thorns starfish and benthos Manta Tow Data (Great Barrier Reef). https://doi.org/10.25845/5c09b0abf315a.
%
% 1- SampledYear: year of monitoring (from original tow date)
% 2- ReefID: Reef IDs beginning 10- to 24- are generally within the Marine Park, 
%       those beginning 26- and 27- are from the Moreton Bay area and those beginning 99- are from Torres Straight
% 3- MngmtRegion: GBRMPA designation of management areas (Far Northern, Cairns/Cooktown, Townsville/Whitsunday , Mackay/Capricorn 
% 4- N: number of conducted tows
% 5- MEAN_COUNT: mean count per tow (number of counted CoTS divided by the number of tows)

% Create boostrap samples by resampling the reef-level mean counts
nb_simul = 500;
BOOTSTRAP_COUNTS = nan(size(LTMP_MTC,1), nb_simul);

for c=1:nb_simul
    
    c
    n = 1;

    while n < size(LTMP_MTC,1)
        
        I =find(LTMP_MTC.SampledYear == LTMP_MTC.SampledYear(n)); % select all surveys conducted in a given year
        [BOOTSTRAP_COUNTS(I,c), J] = datasample(LTMP_MTC.MEAN_COUNT(I), length(I)); % sample among those surveys (with replacement)
        n = n + length(I);
    end
end

save('BOOTSTRAP_MTC.mat', 'BOOTSTRAP_COUNTS') 
% Each columnn contains a bootstrapped sample of every surveyed year
