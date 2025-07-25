# *******************************************************************************************************
# Program Name:                 cm.R
# Project:                      TRA-025
# Purpose:                      Create SDTM Dataset of CM
# Original Author:              Binod Jung Bogati (linkedin.com/in/bjungbogati)
# Copyright:                    © 2025. Unauthorized distribution or reuse prohibited.
# Date Created:                 7/23/2025
# Parameters:                   NA
#
# Input:                        raw_cm dm
# Output:                       cm suppcm
# Modifications:
#   Date        By             Changes
# **********  *************  *************************************
#
# ********************************************************************************************************/

library(tidyverse)
library(labelled)
library(diffdf)

load(url("https://bit.ly/44UDUFG"))

# set formats
route_fmt <- c("Intravenous" = "INTRAVENOUS", "Oral" = "ORAL")
unit_fmt <- c("mL" = "mL", "Tablet" = "TABLET")
freq_fmt <- c("QD = 1x per day" = "QD", "X1 = once" = "ONCE")


# study day
cal_days <- function(date, ref_date) {
  case_when(
    is.na(date) | is.na(ref_date) ~ NA_integer_, # Handle missing dates
    as.Date(date) < as.Date(ref_date) ~ as.numeric(as.Date(date) - as.Date(ref_date)),
    TRUE ~ as.numeric(as.Date(date) - as.Date(ref_date)) + 1L
  )
}

# drop variable
raw_cm_drop <- raw_cm |>
  select(
    -c(
      "PROJECTID", "PROJECT", "STUDYID", "SITEGROUP", "SITENUMBER", "SITEID",
      "INSTANCENAME", "INSTANCEREPEATNUMBER", "FOLDER", "FOLDERID", "FOLDERSEQ",
      "FOLDERNAME", "PAGEREPEATNUMBER", "DATAPAGEID", "RECORDDATE", "RECORDID",
      "MINCREATED", "MAXUPDATED", "SITEID", "ENVIRONMENTNAME", "STUDYENVSITENUMBER",
      "SAVETS", "SITE", "SUBJECTID", "INSTANCEID", "RECORDPOSITION", "TARGETDAYS"
    ),
    -starts_with("Z_"), -ends_with("_STD")
  ) |>
  filter(CMTRT != "")

cm_data <- raw_cm_drop |>
  mutate(
    STUDYID = "TRA-025",
    DOMAIN = "CM",
    USUBJID = str_c(STUDYID, SUBJECT, sep = "-"),
    CMSTDTC = case_when(
      str_detect(CMSTDT_RAW, fixed("UNK")) ~ str_sub(CMSTDT_RAW, 8, 11),
      TRUE ~ CMSTDT_RAW |> as.Date(format = "%d %b %Y") |> as.character()
    ),
    CMENDTC = case_when(
      str_detect(CMENDT_RAW, fixed("UN")) ~ str_replace(CMENDT_RAW, "UN", "01") |>
        as.Date(format = "%d %b %Y") |> str_sub(1, 7),
      TRUE ~ CMENDT_RAW |> as.Date(format = "%d %b %Y") |> as.character()
    ),
    CMROUTE = route_fmt[CMROUTE] |> unname(),
    CMDOSU = unit_fmt[CMUNIT] |> unname(),
    CMDOSFRQ = freq_fmt[CMFREQ] |> unname(),
    # CMINDC = CMINDIC,
    CMDOSE = as.numeric(CMDOSE),
    CMENRF = if_else(CMONGO == "Yes", "ONGOING", ""),
    CMCAT = str_to_upper(DATAPAGENAME)
  ) |>
  left_join(dm, by = "USUBJID") |>
  mutate(
    CMSTDY = cal_days(CMSTDTC, RFSTDTC),
    CMENDY = cal_days(CMENDTC, RFSTDTC),
  ) |>
  arrange(STUDYID, USUBJID, CMTRT, CMSTDTC) |>
  group_by(USUBJID) |>
  mutate(CMSEQ = row_number()) |>
  ungroup()


cm <- cm_data |>
  select(
    STUDYID,
    DOMAIN,
    USUBJID,
    CMSEQ,
    CMTRT,
    CMCAT,
    CMDOSE,
    CMDOSU,
    CMDOSFRQ,
    CMROUTE,
    CMSTDTC,
    CMENDTC,
    CMSTDY,
    CMENDY,
    CMENRF
  )

# Create a named vector of labels
variable_labels <- c(
  STUDYID = "Study Identifier",
  DOMAIN = "Domain Abbreviation",
  USUBJID = "Unique Subject Identifier",
  CMSEQ = "Sequence Number",
  CMTRT = "Reported Name of Drug, Med, or Therapy",
  CMCAT = "Category for Medication",
  CMDOSE = "Dose",
  CMDOSU = "Dose Units",
  CMDOSFRQ = "Dose Frequency per Interval",
  CMROUTE = "Route of Administration",
  CMSTDTC = "Start Date/Time of Medication",
  CMENDTC = "End Date/Time of Medication",
  CMSTDY = "Study Day of Start of Medication",
  CMENDY = "Study Day of End of Medication",
  CMENRF = "End Relative to Reference Period"
)

# Apply labels
cm <- cm |>
  set_variable_labels(.labels = variable_labels)


diffdf(base = cm, compare = v_cm)


suppcm <- cm_data |>
  mutate(
    CMDICVER = str_c(CMTRT_CODERDICTNAME, CMTRT_CODERDICTVERSION, sep = " "),
    RDOMAIN = "CM",
    IDVAR = "CMSEQ",
    ATC = CMTRT_ATC,
    ATCCD = CMTRT_ATC_CODE
  ) |>
  pivot_longer(
    cols = c(ATC, ATCCD, CMDICVER),
    names_to = "QNAM",
    values_to = "QVAL"
  ) |>
  mutate(
    QLABEL = case_when(
      QNAM == "ATC" ~ "ATC Level 1 Term",
      QNAM == "ATCCD" ~ "ATC Level 1 Code",
      QNAM == "CMDICVER" ~ "Medical Coding Dictionary & Version",
      # QNAM == "CMYN" ~ "CM Taken?"
    ),
    QORIG = if_else(QNAM %in% "CMYN", "CRF", "eDT"),
    QEVAL = ""
  ) |>
  select(
    STUDYID,
    RDOMAIN,
    USUBJID,
    IDVAR,
    IDVARVAL = CMSEQ,
    QNAM,
    QLABEL,
    QVAL,
    QORIG,
    QEVAL
  ) |>
  mutate(
    IDVARVAL = as.character(IDVARVAL)
  ) |>
  arrange(STUDYID, USUBJID, IDVAR, IDVARVAL, QNAM)



diffdf(base = suppcm, compare = v_suppcm)
