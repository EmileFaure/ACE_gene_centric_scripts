01-ACE CTD Profiles
================
Lisa-Marie Delpech
2/23/2022

- [Compute absolute salinity and potential
  temperature](#compute-absolute-salinity-and-potential-temperature)
- [T-S plots per event](#t-s-plots-per-event)
- [Oxygen-Salinity and Oxygen-Temperature
  plots](#oxygen-salinity-and-oxygen-temperature-plots)

``` r
require(tidyverse)
require(viridis)
require(gsw)
require(oce)
require(ocedata)
require(ggoce)
require(cowplot)
require(ggrepel)
```

``` r
ctd.profile.data <- readRDS("./R_Data/CTD_profile_data.rds") # CTD data concatenated in long format df
head(ctd.profile.data)
```

    ##   EVENTNBR STNNBR CASTNO        DATETIME_UTC  LATITUDE LONGITUDE CTDPRS DEPTH
    ## 1     3117    103     14 2017-03-16 10:20:01 -43.97961  14.06983      1 0.990
    ## 2     3117    103     14 2017-03-16 10:20:01 -43.97961  14.06983      2 1.976
    ## 3     3117    103     14 2017-03-16 10:20:01 -43.97961  14.06983      3 2.967
    ## 4     3117    103     14 2017-03-16 10:20:01 -43.97961  14.06983      4 3.987
    ## 5     3117    103     14 2017-03-16 10:20:01 -43.97961  14.06983      5 4.949
    ## 6     3117    103     14 2017-03-16 10:20:01 -43.97961  14.06983      6 5.949
    ##    CTDTMP CTDTMP_FLAG_W  CTDSAL CTDSAL_FLAG_W  CTDDENS  CTDOXY CTDOXY_FLAG_W
    ## 1 10.7156             4  9.9693             4 1007.385 273.358             4
    ## 2 10.6822             2 34.3687             3 1026.350 255.778             1
    ## 3 10.6812             2 34.3729             2 1026.358 256.259             1
    ## 4 10.6818             2 34.3727             2 1026.362 256.139             1
    ## 5 10.6809             2 34.3735             2 1026.367 256.070             1
    ## 6 10.6750             2 34.3742             2 1026.373 256.197             1
    ##   CTDOXYSAT CTDOXYSAT_FLAG_W CTDFLUOR1Q CTDFLUOR2Q     PAR CASTTYPE
    ## 1  323.3252                4       -999  0.4282905  363.54   upcast
    ## 2  271.8630                1       -999  0.4282905 -999.00   upcast
    ## 3  271.8602                1       -999  0.4282905 -999.00   upcast
    ## 4  271.8571                1       -999  0.4282905 -999.00   upcast
    ## 5  271.8609                1       -999  0.4282905  363.66   upcast
    ## 6  271.8937                1       -999  0.4282905  363.66   upcast
    ##                                                     FILE
    ## 1 ./Metadata/ACE_CTD_Profile_Data/uACE201603_014_ct1.csv
    ## 2 ./Metadata/ACE_CTD_Profile_Data/uACE201603_014_ct1.csv
    ## 3 ./Metadata/ACE_CTD_Profile_Data/uACE201603_014_ct1.csv
    ## 4 ./Metadata/ACE_CTD_Profile_Data/uACE201603_014_ct1.csv
    ## 5 ./Metadata/ACE_CTD_Profile_Data/uACE201603_014_ct1.csv
    ## 6 ./Metadata/ACE_CTD_Profile_Data/uACE201603_014_ct1.csv

``` r
length(table(ctd.profile.data$EVENTNBR))
```

    ## [1] 61

``` r
ctd.profile.data$EVENTNBR <- as.factor(ctd.profile.data$EVENTNBR)

meta <- read.table(file = "./Metadata/meta_CTD_ess_V4_4.csv", header=TRUE, sep = ";", dec = ".", stringsAsFactors = TRUE) # Metadata to get event numbers
length(table(meta$Event_number)) # 49 event numbers vs. 61 available CTD profiles
```

    ## [1] 49

``` r
meta.sum <- meta %>% 
  select("Event_number","Depth_m","Latitude","Longitude","Temperature_ITS90","Salinity_PSS78","Density_kg.m2","dO2_µmol.kg", "Water_mass_detailed") %>% 
  group_by(Event_number, Depth_m, Water_mass_detailed) %>% # Remove filter information 
  summarize_if(.predicate = is.numeric, .funs = mean)
```

The README file accompanying the metadata in Henry *et al.* (2020)
describes quality flags as follow:

> Quality Flags Each variable was inspected and quality flags
> (WHP-Exchange Format, Barna, Swift and Diggs, 2016;
> <https://exchange-format.readthedocs.io/en/latest/>) applied using the
> following conditions: \* Flag 1 Not Calibrated - Oxygen variables
> (CTD) flagged until comparison with bottle samples or suitable
> climatology dataset e.g. GLODAP, WOCE or WOA  
> \* Flag 2 Acceptable Measurement  
> \* Flag 3 Questionable Measurement - CTDTEMP or CTDSAL failed spike or
> gradient tests (SeaDataNet Data Quality Control Procedures Version
> 2.0). BACKSC where \> 50 % of values flagged as Flag 4, see below.  
> \* Flag 4 Bad Measurement - Bad flag applied during Sea-Bird Data
> Processing Modules, scan variable value is negative number, variable
> value identified as being anomalous at beginning or end of upcast or
> downcast during visual inspection (see ace_ctd_visual_inspection.csv
> for full list of data points identified), values fall outside
> acceptable sensor specific ranges (see information on sensor specific
> flagging below).  
> \* Flag 5 Not Reported - No instances where this was applied. \* Flag
> 6 Interpolated over pressure interval larger than 2dbar. \* Flag 7
> Despiked - No instances where this was applied.  
> \* Flag 9 Not Sampled - Sensor was not connected. Flags can be
> identified from the following header format \[PARAMETER\]\_FLAG_W (WHP
> Exchange Format).

We here remove flags 3 and 4 for salinity, temperature and oxygen.

``` r
CTDTMP_QUAL <- ctd.profile.data$CTDTMP
CTDTMP_QUAL[which(ctd.profile.data$CTDTMP_FLAG_W %in% c(3,4))] <- NA
CTDSAL_QUAL <- ctd.profile.data$CTDSAL
CTDSAL_QUAL[which(ctd.profile.data$CTDSAL_FLAG_W %in% c(3,4))] <- NA
CTDOXY_QUAL <- ctd.profile.data$CTDOXY
CTDOXY_QUAL[which(ctd.profile.data$CTDOXY_FLAG_W %in% c(3,4))] <- NA

ctd.profile.filtered <- ctd.profile.data %>% 
  mutate(CTDTMP = CTDTMP_QUAL) %>% # , .after = CTDTMP
  mutate(CTDSAL = CTDSAL_QUAL) %>% # , .after = CTDSAL
  mutate(CTDOXY = CTDOXY_QUAL) # , .after = CTDOXY
```

Check availability of CTD profiles for metadata and metagenomic samples.
These are the events not available in in CTD profiles.

``` r
meta.sum[which(!meta.sum$Event_number %in% ctd.profile.filtered$EVENTNBR),]
```

    ## # A tibble: 3 × 9
    ## # Groups:   Event_number, Depth_m [3]
    ##   Event_number Depth_m Water_mass_detailed Latitude Longitude Temperature_ITS90
    ##          <int>   <int> <fct>                  <dbl>     <dbl>             <dbl>
    ## 1           75      15 SASW                   -46.9      38.1                NA
    ## 2           75      60 SASW                   -46.9      38.1                NA
    ## 3           75     150 SASW                   -46.9      38.1                NA
    ## # ℹ 3 more variables: Salinity_PSS78 <dbl>, Density_kg.m2 <dbl>,
    ## #   dO2_µmol.kg <dbl>

Event 75 is not available as a CTD profile while it is in the metadata.
Hence we will not be able to classify this sample’s water mass based on
the CTD data. Note that event 317 is not in the metagenomes dataset but
included here with all available CTD profiles for the metadata.

# Compute absolute salinity and potential temperature

``` r
ctd <- ctd.profile.filtered %>% 
  mutate(Absolute_Salinity = gsw_SA_from_SP(CTDSAL, CTDPRS, LONGITUDE, LATITUDE)) %>% 
  mutate(Potential_temperature = gsw_pt_from_t(Absolute_Salinity, CTDTMP, CTDPRS, p_ref = 0),
         Conservative_temperature = gsw_CT_from_t(Absolute_Salinity, CTDTMP, CTDPRS))
```

# T-S plots per event

``` r
# TS plot for all
(plot.TS.all <- ggplot(ctd %>% filter(CASTTYPE == "downcast"),
                       aes(x = CTDSAL, y = Conservative_temperature)) +
   geom_isopycnal(salinity_type = "practical",
                  temperature_type = "conservative",
                  ref_longitude = 0,
                  ref_latitude = 0,
                  breaks = seq(25,28,0.2),
                  linetype = "dotted") +
   geom_path(aes(group = EVENTNBR),
             size = .3,
             col = "#777777") +
   labs(x = "Practical Salinity (psu)", y = "Conservative Temperature (°C)") +
   scale_x_continuous(breaks = seq(32,35,0.2)) +
   scale_y_continuous(breaks = seq(-1,12,2)) +
   theme_bw())
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).

![](01-ACE_CTD_Profiles_files/figure-gfm/Plot%20for%20all%20events-1.png)<!-- -->

``` r
# Depth-T for all
(plot.Depth.T.all <- ggplot(data = ctd %>% filter(CASTTYPE == "downcast",
                                                  DEPTH < 2000),              # Only affects event 1465 for which depth goes over 3000m
                            aes(x = Conservative_temperature, y = DEPTH)) +
    geom_path(aes(group = EVENTNBR),
              size = .3, col = "#777777") +
    scale_x_continuous(breaks = seq(-1,12,2)) +
    scale_y_reverse() +
    labs(x = "Conservative Temperature (°C)", y = "Depth (m)") +
    theme_bw())
```

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).

![](01-ACE_CTD_Profiles_files/figure-gfm/Plot%20for%20all%20events-2.png)<!-- -->

``` r
# Depth-S for all
(plot.Depth.S.all <- ggplot(data = ctd %>% filter(CASTTYPE == "downcast",
                                                  DEPTH < 2000),
                            aes(x = CTDSAL, y = DEPTH)) +
    geom_path(aes(group = EVENTNBR),
              size = .3, col = "#777777") +
    scale_x_continuous(breaks = seq(32,35,0.4)) +
    scale_y_reverse() +
    labs(x = "Practical Salinity (psu)", y = "Depth (m)") +
    theme_bw())
```

    ## Warning: Removed 80 rows containing missing values or values outside the scale range
    ## (`geom_path()`).

![](01-ACE_CTD_Profiles_files/figure-gfm/Plot%20for%20all%20events-3.png)<!-- -->

``` r
for (event in unique(ctd$EVENTNBR)){
  if (event %in% meta.sum$Event_number){
    df.ctd = ctd %>% filter(CASTTYPE == "downcast",
                            ##CTDSAL >= 32 & CTDSAL <= 35,
                            EVENTNBR == event)
    df.meta = meta.sum %>% filter(Event_number == event)
    
    # Add TS plot for one event
    plot.TS <- plot.TS.all +
      geom_path(data = df.ctd,
                aes(x = CTDSAL, y = Conservative_temperature, col = DEPTH),
                size = .7) +
      scale_color_gradient(high = "red", low = "blue") +
      geom_point(data = df.meta,
                 aes(x = Salinity_PSS78, y = Temperature_ITS90, col = Depth_m)) +
      geom_point(data = df.meta,
                 aes(x = Salinity_PSS78, y = Temperature_ITS90),
                 col = "black", fill = "transparent", shape = 21) +
      geom_label_repel(data = df.meta,
                       aes(x = Salinity_PSS78, y = Temperature_ITS90, label = Water_mass_detailed),
                       size = 3,
                       label.padding = 0.15, force_pull = .5,
                       label.r = 0.07) +
      scale_fill_gradient(high = "red", low = "blue", guide = "none") +
      labs(col = "Depth (m)", title = paste("Event number", event)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Add Depth-T plot for one event
    plot.Depth.T <- plot.Depth.T.all +
      geom_path(data = df.ctd %>% filter(DEPTH < 2000),
                aes(x = Conservative_temperature, y = DEPTH, col = DEPTH),
                size = .7) +
      geom_point(data = df.meta,
                 aes(x = Temperature_ITS90, y = Depth_m, col = Depth_m)) +
      scale_color_gradient(high = "red", low = "blue") +
      geom_point(data = df.meta,
                 aes(x = Temperature_ITS90, y = Depth_m),
                 col = "black", fill = "transparent", shape = 21) +
      #scale_fill_gradient(high = "red", low = "blue", guide = "none") +
      geom_label_repel(data = df.meta,
                       aes(x = Temperature_ITS90, y = Depth_m, label = paste(Depth_m,"m")),
                       size = 2) +
      labs(col = "Depth (m)")
    
    # Add Depth-S plot for one event
    plot.Depth.S <- plot.Depth.S.all +
      geom_path(data = df.ctd %>% filter(DEPTH < 2000),
                aes(x = CTDSAL, y = DEPTH, col = DEPTH),
                size = .7) +
      scale_color_gradient(high = "red", low = "blue") +
      geom_point(data = df.meta,
                 aes(x = Salinity_PSS78, y = Depth_m, col = Depth_m)) +
      geom_label_repel(data = df.meta,
                       aes(x = Salinity_PSS78, y = Depth_m, label = paste(Depth_m,"m")),
                       size = 2) +
      geom_point(data = df.meta,
                 aes(x = Salinity_PSS78, y = Depth_m), col = "black", fill = "transparent", shape = 21) +
      #scale_fill_gradient(high = "red", low = "blue", guide = "none") +
      labs(col = "Depth (m)") +
      theme(legend.position = "none")
    
    # All together
    plot.Depth.S.T <- plot_grid(plot.Depth.S, plot.Depth.T, rel_widths = c(1,1.3))
    print(
      plot_grid(plot.TS, plot.Depth.S.T, nrow = 2, ncol = 1)
    )
    
    #ggsave(paste0("./R_Plots/CTD_Profiles/CTD_profile_",event,".pdf"), height = 10, width = 8, paper = "a4")
  }
}
```

![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-1.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-2.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-3.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-4.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-5.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-6.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-7.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-8.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-9.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-10.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-11.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-12.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-13.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-14.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-15.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-16.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-17.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-18.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-19.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-20.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-21.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-22.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-23.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-24.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-25.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-26.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-27.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-28.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-29.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-30.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-31.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-32.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-33.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-34.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-35.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-36.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-37.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-38.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-39.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-40.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-41.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-42.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-43.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-44.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-45.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-46.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-47.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Add%20highlight%20for%20each%20event-48.png)<!-- -->

Note that event 1029 has a limited downcast but the upcast might be
useful

``` r
# Refine salinity between 33.6 and 34.8
# And temperature between -2 and 5

plot.TS.all.zoom <- ggplot(ctd %>% filter(CASTTYPE == "downcast",
                                          CTDSAL >= 33.6, CTDSAL <= 34.8,
                                          Conservative_temperature >=-2, Conservative_temperature <= 5),
                           aes(x = CTDSAL, y = Conservative_temperature)) +
  geom_isopycnal(salinity_type = "practical",
                 temperature_type = "conservative",
                 ref_longitude = 0,
                 ref_latitude = 0,
                 breaks = seq(26.7,28,0.1),
                 linetype = "dotted") +
  geom_path(aes(group = EVENTNBR),
            size = .3,
            col = "#777777") +
  labs(x = "Practical Salinity (psu)", y = "Conservative Temperature (°C)") +
  scale_x_continuous(breaks = seq(33.6,34.8,0.1)) +
  scale_y_continuous(breaks = seq(-2,5,1)) +
  theme_bw()

for (event in unique(ctd %>% filter(CASTTYPE == "downcast", CTDSAL >= 33.6, CTDSAL <= 34.8, Conservative_temperature >=-2, Conservative_temperature <= 5) %>% pull(EVENTNBR))){
  if (event %in% meta.sum$Event_number){
    
    df.ctd = ctd %>% filter(CASTTYPE == "downcast",
                            CTDSAL >= 33.6, CTDSAL <= 34.8,
                            Conservative_temperature >=-2, Conservative_temperature <= 5,
                            EVENTNBR == event)
    df.meta = meta.sum %>% filter(Event_number == event,
                                  Salinity_PSS78 >= 33.6, Salinity_PSS78 <= 34.8,
                                  Temperature_ITS90 >= -2, Temperature_ITS90 <= 5)
    
    # Add TS plot for one event
    plot.TS.path <- plot.TS.all.zoom +
      geom_path(data = df.ctd,
                aes(x = CTDSAL, y = Conservative_temperature, col = DEPTH),
                size = .7) +
      scale_color_gradient(high = "red", low = "blue") +
      labs(col = "Depth (m)", title = paste("Zoom - Event number", event)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    if (nrow(df.meta) != 0){
      print(
      plot.TS.path +
        geom_point(data = df.meta,
                   aes(x = Salinity_PSS78, y = Temperature_ITS90, col = Depth_m)) +
        geom_point(data = df.meta,
                   aes(x = Salinity_PSS78, y = Temperature_ITS90),
                   col = "black", fill = "transparent", shape = 21) +
        scale_fill_gradient(high = "red", low = "blue", guide = "none") +
        geom_label_repel(data = df.meta,
                         aes(x = Salinity_PSS78, y = Temperature_ITS90, label = paste(Depth_m,"m")),
                         size = 2) +
        labs(col = "Depth (m)", title = paste("Zoom - Event number", event)) +
        theme(plot.title = element_text(hjust = 0.5))
      )
    }
  }
  
  #ggsave(paste0("./R_Plots/CTD_Profiles/TS_zoom_",event,".pdf"), height = 5, width = 8, paper = "a4")
}
```

![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-1.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-2.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-3.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-4.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-5.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-6.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-7.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-8.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-9.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-10.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-11.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-12.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-13.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-14.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-15.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-16.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-17.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-18.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-19.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-20.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-21.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-22.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-23.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-24.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-25.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-26.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-27.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-28.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-29.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-30.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-31.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-32.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-33.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-34.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-35.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-36.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-37.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-38.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-39.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Zoom%20in%20on%20TS%20diagram-40.png)<!-- -->

# Oxygen-Salinity and Oxygen-Temperature plots

In order to determine water masses and remove some ambiguity around some
samples, we need to also add salinity oxygen plots to our data.

Notes from the data collection:

> A note on oxygen: During the data conversion process, the ‘tau’ and
> ‘hysteresis’ corrections were applied to the oxygen data. The ‘tau’
> correction improves the response of the measured signal in regions of
> large oxygen gradients, and also may amplify residual noise in the
> signal (CalCOFI CTD Data Algorithms). The ‘tau’ correction computes
> the derivative of the oxygen signal with respect to time with a
> user-input window size), using a linear regression to determine the
> slope. The correction was applied using the Data Conversion step where
> the module uses a window looking backwards in time and a window size
> of two seconds was defined. Although it is recommended that the ‘tau’
> correction be applied during the ‘Derive’ module for more accurate
> results, a subset of files were processed using that processing order
> and there was no difference in the derived oxygen values. The
> hysteresis correction corrects the oxygen voltage values for changes
> in membrane permeability as pressure varies. Correction coefficients
> stored in the XMLCON file are applied with calculations based on the
> current pressure and how long the session spent at previous pressures.
> Two oxygen variables are provided: (1) dissolved oxygen concentration
> in micromols/kg and (2) dissolved oxygen saturation, the theoretical
> saturation limit of the water at the local temperature and salinity
> value, with local pressure reset to zero (1 atmosphere) (See SBE Data
> Processing Software Manual for SBE Data Processing Software v7.26.8 &
> CalCOFI CTD Data Algorithms). Dissolved oxygen saturation is
> calculated using the equations of Garcia and Gordon (1992) which is
> based on the equation of Weiss (1970) but reduced the error in cold
> waters and presents improved performance in waters between -5 to 50
> degrees Celsius.

Regarding the quality of those values, Henry *et al.* report:

> Variables CTDOXY and CTDOXYSAT were flagged with Flag 1 as per the CTD
> cast data (See section Quality Flags).  
> The oxygen variables have not been corrected further, but should be
> compared to oxygen bottle data when made available. An initial
> comparison to World Ocean Atlas 2018 data (WOA18; Figure 3; Garcia *et
> al.*, 2018) suggests a substantial offset which should be corrected
> for when using the oxygen data.

``` r
# Oxygen-Salinity diagram for all samples
(plot.dOS.all <- ggplot(ctd %>% filter(CASTTYPE == "downcast"),
                      aes(x = CTDSAL, y = CTDOXY)) +
  geom_path(aes(group = EVENTNBR),
            size = .3,
            col = "#777777") +
  labs(x = "Practical Salinity (psu)", y = "Dissolved Oxygen (µmol/kg)") +
  scale_x_continuous(breaks = seq(32,35,0.2)) +
  scale_y_continuous(breaks = seq(150,400,20)) +
  theme_bw())
```

![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-1.png)<!-- -->

``` r
#Temperature-Oxygen diagram for all samples
(plot.TdO.all <- ggplot(ctd %>% filter(CASTTYPE == "downcast"),
                        aes(x = CTDOXY, y = Conservative_temperature)) +
  geom_path(aes(group = EVENTNBR),
            size = .3,
            col = "#777777") +
  labs(x = "Dissolved Oxygen (µmol/kg)", y = "Conservative Temperature (°C)") +
  scale_x_continuous(breaks = seq(150,400,20)) +
  scale_y_continuous(breaks = seq(-1,12,2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
```

![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-2.png)<!-- -->

``` r
#Depth-Oxygen for all samples
(plot.Depth.O2.all <- ggplot(data = ctd %>% filter(CASTTYPE == "downcast",
                                                  DEPTH < 2000),
                             aes(x = CTDOXY, y = DEPTH)) +
  geom_path(aes(group = EVENTNBR),
            size = .3, col = "#777777") +
  scale_x_continuous(breaks = seq(150,400,20)) +
  scale_y_reverse() +
  labs(x = "Dissolved Oxygen (µmol/kg)", y = "Depth (m)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
```

![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-3.png)<!-- -->

``` r
for (event in unique(ctd$EVENTNBR)){
  if (event %in% meta.sum$Event_number){
  #print(event)
  df.ctd = ctd %>% filter(CASTTYPE == "downcast",
                          EVENTNBR == event)
  df.meta = meta.sum %>% filter(Event_number == event)
  
  # Add TS plot for one event
  plot.dOS <- plot.dOS.all +
    geom_path(data = df.ctd,
              aes(x = CTDSAL, y = CTDOXY, col = DEPTH),
              size = .7) +
    scale_color_gradient(high = "red", low = "blue") +
    geom_point(data = df.meta,
               aes(x = Salinity_PSS78, y = dO2_µmol.kg, col = Depth_m)) +
    geom_point(data = df.meta,
               aes(x = Salinity_PSS78, y = dO2_µmol.kg),
               col = "black", fill = "transparent", shape = 21) +
    scale_fill_gradient(high = "red", low = "blue", guide = "none") +
    labs(col = "Depth (m)", title = paste("Event number", event)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_label_repel(data = df.meta,
                     aes(x = Salinity_PSS78, y = dO2_µmol.kg, label = paste(Depth_m,"m")),
                     size = 2)
  
  plot.TdO <- plot.TdO.all +
    geom_path(data = df.ctd,
              aes(x = CTDOXY, y = CTDTMP, col = DEPTH),
              size = .7) +
    scale_color_gradient(high = "red", low = "blue") +
    geom_point(data = df.meta,
               aes(x = dO2_µmol.kg, y = Temperature_ITS90, col = Depth_m)) +
    geom_point(data = df.meta,
               aes(x = dO2_µmol.kg, y = Temperature_ITS90),
               col = "black", fill = "transparent", shape = 21) +
    scale_fill_gradient(high = "red", low = "blue", guide = "none") +
    labs(col = "Depth (m)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_label_repel(data = df.meta,
                     aes(x = dO2_µmol.kg, y = Temperature_ITS90, label = paste(Depth_m,"m")),
                     size = 2) +
    theme(legend.position = "none")
  
  plot.Depth.O2 <- plot.Depth.O2.all +
    geom_path(data = df.ctd %>% filter(DEPTH < 2000),
              aes(x = CTDOXY, y = DEPTH, col = DEPTH),
              size = .7) +
    scale_color_gradient(high = "red", low = "blue") +
    geom_point(data = df.meta,
               aes(x = dO2_µmol.kg, y = Depth_m, col = Depth_m)) +
    geom_label_repel(data = df.meta,
                     aes(x = dO2_µmol.kg, y = Depth_m, label = paste(Depth_m,"m")),
                     size = 2) +
    geom_point(data = df.meta,
               aes(x = dO2_µmol.kg, y = Depth_m), col = "black", fill = "transparent", shape = 21) +
    labs(col = "Depth (m)")
  
  
  # All together
  plot2 <- plot_grid(plot.TdO, plot.Depth.O2, rel_widths = c(1.1,1))
  print(
    plot_grid(plot.dOS, plot2, nrow = 2, ncol = 1)
  )
  #ggsave(paste0("./R_Plots/CTD_Profiles/CTD_Oxygen_",event,".pdf"), height = 11, width = 8, paper = "a4")
  }
}
```

![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-4.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-5.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-6.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-7.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-8.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-9.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-10.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-11.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-12.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-13.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-14.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-15.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-16.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-17.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-18.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-19.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-20.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-21.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-22.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-23.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-24.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-25.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-26.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-27.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-28.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-29.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-30.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-31.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-32.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-33.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-34.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-35.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-36.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-37.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-38.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-39.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-40.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-41.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-42.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-43.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-44.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-45.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-46.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-47.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-48.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-49.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-50.png)<!-- -->![](01-ACE_CTD_Profiles_files/figure-gfm/Plots%20dO-S-51.png)<!-- -->
