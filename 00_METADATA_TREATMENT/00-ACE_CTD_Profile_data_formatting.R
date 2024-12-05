# ------------ Libraries ------------
require(tidyverse)
require(lubridate)

# ------------ Load data ------------
files <- list.files(path = "./Metadata/ACE_CTD_Profile_Data", pattern = "*.csv", full.names = TRUE)

# ------------ Process ------------
ctd.profile.data <- tibble()
ctd.headers <- tibble()

for (file in files){
  ctd <- read.table(file, header = TRUE, skip = 17, sep = ",") # Import CTD data
  info <- read.table(file, header = FALSE, skip = 5, nrows = 12, sep = ",") # Import other info (inc event number)
  
  Name <-  str_split(as.vector(info$V1), " = ", simplify = TRUE)[,1]
  Attribute <-  str_split(as.vector(info$V1), " = ", simplify = TRUE)[,2]
  
  ctd.info <- data.frame(Name, Attribute) %>% 
    column_to_rownames(var = "Name")
  
  ctd.profile.data <- ctd[-1,] %>% # Remove column info (units)
    rownames_to_column(var = "ROWNBR") %>% # Get rid of row numbering
    mutate(DATETIME_UTC = ymd_hms(DATETIME_UTC)) %>% 
    select(c("STNNBR","CASTNO","DATETIME_UTC","LATITUDE","LONGITUDE","CTDPRS","DEPTH","CTDTMP","CTDTMP_FLAG_W","CTDSAL","CTDSAL_FLAG_W","CTDDENS","CTDOXY","CTDOXY_FLAG_W","CTDOXYSAT","CTDOXYSAT_FLAG_W","CTDFLUOR1Q","CTDFLUOR2Q","PAR")) %>% # Remove FLAG columns (mainly)
    mutate(EVENTNBR = as.numeric(ctd.info["EVENTNBR",]), .before = STNNBR) %>% # Append event number in long format
    mutate_if(.predicate = is.character, .funs = as.numeric) %>% 
    mutate(CASTTYPE = as.factor(ctd.info["CASTTYPE",])) %>% 
    mutate(FILE = file) %>% 
    bind_rows(ctd.profile.data)
  
  ctd.info <- ctd.info %>% 
    t() %>%
    as.data.frame()
  rownames(ctd.info) <- NULL
  ctd.headers <- ctd.info %>% 
    select(-c("NUMBER_HEADERS")) %>% 
    relocate("EVENTNBR", .before = "EXPOCODE") %>%
    mutate(BOTTOM_DATE = ymd(BOTTOM_DATE),
           BOTTOM_TIME = paste(substr(.$BOTTOM_TIME, 1,2),substr(.$BOTTOM_TIME, 3,4), sep = ":"),
           BOTTOM_TIME = hm(BOTTOM_TIME)) %>%
    bind_rows(ctd.headers)
  
}

ctd.profile.data
ctd.headers

# ------------ Save RDS ------------
saveRDS(ctd.profile.data, file = "./R_Data/CTD_profile_data.rds")
write.table(ctd.headers, file = "./Metadata/ACE_CTD_Profile_Data/CTD_profile_headers.tsv", sep = "\t", row.names= FALSE)


