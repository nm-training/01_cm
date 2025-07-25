# R Trainer: Binod Jung Bogati (linkedin.com/in/bjungbogati)
# copyright: © 2025. Unauthorized distribution or reuse prohibited.


# load
library(tidyverse)


load(url("https://bit.ly/44UDUFG"))

names(raw_cm)


raw_cm2 <- raw_cm |> 
  select(-c("PROJECTID","PROJECT","STUDYID","SITEGROUP","SITENUMBER","SITEID",
            "INSTANCENAME","INSTANCEREPEATNUMBER","FOLDER","FOLDERID","FOLDERSEQ",
            "FOLDERNAME","PAGEREPEATNUMBER","DATAPAGEID","RECORDDATE","RECORDID",
            "MINCREATED","MAXUPDATED","SITEID","ENVIRONMENTNAME","STUDYENVSITENUMBER",
            "SAVETS","SITE","SUBJECTID","INSTANCEID","RECORDPOSITION","TARGETDAYS"), 
         -starts_with("Z_"))

# Exercise 1: Write a program to map proper formats for CMSTDTC and CMENDTC from raw CM dataset.

raw_ex1 <- raw_cm |>
  select(SUBJECT, CMSTDT_RAW, CMENDT_RAW) |>
  mutate(
    CMSTDTC = case_when(
      str_detect(CMSTDT_RAW, fixed("UNK")) ~ str_sub(CMSTDT_RAW, 8, 11),
      TRUE ~ dmy(CMSTDT_RAW) |> as.character()
    ),
    CMENDTC = case_when(
      str_detect(CMENDT_RAW, fixed("UN")) ~ str_replace(CMENDT_RAW, fixed("UN"), "01") |>
        dmy() |> as.character() |> str_sub(1, 7),
      TRUE ~ dmy(CMENDT_RAW) |> as.character()
    ), 
    
    USUBJID = str_c("TRA-025", SUBJECT,sep = "-")
  ) |> select(USUBJID, CMSTDTC, CMENDTC) 



# Exercise 2: Write a program to create Unique Subject Identifier variable and Sequence variable.

raw_ex1 <- raw_cm |>
  USUBJID = str_c("TRA-025", SUBJECT,sep = "-") |>
  arrange(STUDYID, USUBJID, CMTRT, CMSTDTC) |>
  group_by(USUBJID) |>
  mutate(CMSEQ = row_number()) |>
  ungroup()
  

# Exercise 3: Write a program to calculate Study Day variable. (CMSTDTC, CMENDTC)

raw_ex3 <- raw_ex1 |> 
  right_join(dm, by = "USUBJID") |> 
  mutate(
    CMSTDT = as.Date(CMSTDTC), 
    CMENDT = as.Date(CMENDTC),
    RFSTDT = as.Date(RFSTDTC), 
    
    CMSTDY = case_when(
      CMSTDT < RFSTDT ~ as.numeric(CMSTDT - RFSTDT), 
      CMSTDT >= RFSTDT ~ as.numeric(CMSTDT - RFSTDT) + 1
    )
  ) |> select(USUBJID, CMSTDTC, CMENDTC, CMSTDY)




# Exercise 4: Take the concomitant medication data.
# The sample source data also has ATC coded variables for the medication names.
# Plan the mapping of these into cm_mapping.xlsx excel file




# Exercise 5: Compare the mapping you proposed for exercise 4 with the mapping file given by Binod.


#
# Exercise 6. After completing the comparison, create a single cm.R program to create CM and SUPPCM .
# After that, validate the CM & SUPPCM datasets with provided dataset (upon request).
#
