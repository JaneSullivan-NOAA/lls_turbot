# AFSC longline survey (LLS) data for Greenland turbot

The AFSC longline survey (LLS) [database documentation](https://akfinbi.psmfc.org/analyticsRes/Documentation/Database_Background_Instructions_AKFIN_20210915.pdf) (hosted on AKFIN Answers and AKFIN's new assessment portal ASAP) is a great place to start for understanding the available data for this survey. The AFSC also publishes an [annual survey report](https://repository.library.noaa.gov/view/noaa/37523), which includes several useful references on survey history and methods. Maps and detailed descriptions of geographic areas and depth strata can be found in [Echave et al. 2013](https://repository.library.noaa.gov/view/noaa/11869). 

![](https://github.com/JaneSullivan-NOAA/lls_turbot/blob/main/data/lls_stations_echave.PNG)

Figure reproduced from [Echave 2017](https://media.fisheries.noaa.gov/dam-migration/first-results-of-the-tagging-of-shortspine-thornyhead-508.pdf)

The indices of abundance available for the LLS include the following: CPUE = numbers per skate (1 skate = 45 hooks); relative population numbers (RPNs) = area-weighted CPUE (area estimates are defined by depth strata and geographic area), relative population weight (RPW) = RPN multiplied by mean fish weight, which is calculated using an allometric relationship and the mean length of fish collected in a given strata and geographic area. The methods for variance estimation of these indices are documented on [p 26 the 2016 sablefish SAFE pdf](https://apps-afsc.fisheries.noaa.gov/REFM/Docs/2016/GOAsablefish.pdf) (p 350 of the GOA SAFE). 

The time series starts for the domestic (i.e. operated by the AFSC) LLS in the BS and AI starts in 1996, and there are historical data available in the BS/AI from the cooperative Japanese/U.S. survey. In the modern/domestic survey, the eastern AI are surveyed in even years, and the BS is surveyed in odd years. The longline survey does not sample the western AI. Estimates in NW and SW AI are based on fixed ratios in the NE and SE AI, respectively, from historical cooperative Japanese/U.S. surveys in 1979-1994. Note that definitions of eastern and western AI in the LLS are not equivalent to the AI bottom trawl survey definitions. It is a work in progress to reconcile these survey areas for purposes of comparison. Survey depths range from ~100 to 1,000 m; depth strata are defined using 100 m increments at depths < 400 m and 200 m increments at depths >= 400 m. These depth strata differ from the AFSC bottom trawl survey definitions, making one-on-one comparisons challenging at depths > 300 m.

![](https://github.com/JaneSullivan-NOAA/lls_turbot/blob/main/data/depth_strata.PNG)



