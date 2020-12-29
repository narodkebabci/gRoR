#' @title Build duplets based on mutations
#'
#' @description This function takes the data frame and mutation column to form WT-MT duplets. It also
#' prints the missing duplets within the dataset.
#'
#' @param data A data frame containing mutations that will be classified
#' @param mutations Specify the column name that contains the mutations
#'
#' @examples
#' duplet_formation(data = df, mutations = "Mutation")
#'
#' @export

duplet_formation <- function(data, mutations){

  # extract the wt in a vector
  WT <- vector("character", length = nrow(data))

  for (i in 1:nrow(data)){
    WT[i] <- stri_sub(data[, mutations][i], from=-stri_length(data[, mutations][i]),
                           to=-stri_length(data[, mutations][i]))
  }

  # extract the mt in a vector
  MT <- vector("character", length = nrow(data))

  for (j in 1:nrow(data)){
    MT[j] <- stri_sub(data[, mutations][j], from=stri_length(data[, mutations][j]),
                           to=stri_length(data[, mutations][j]))
  }

  # create a vector of duplets
  Duplets <- gsub(" ","",paste(WT, MT))

  # add Duplets vector to the original dataframe
  data[, "Duplets"] <- Duplets

  # relocate the position of Duplets vector
  data <- data %>% relocate(Duplets, .after = mutations)

  # list duplet labels, return missing
  duplet_labels = c("AC","AD","AE","AF","AG","AH","AI","AK","AL","AM","AN","AP","AQ","AR","AS","AT","AV","AW","AY",
                    "CA","CD","CE","CF","CG","CH","CI","CK","CL","CM","CN","CP","CQ","CR","CS","CT","CV","CW","CY",
                    "DA","DC","DE","DF","DG","DH","DI","DK","DL","DM","DN","DP","DQ","DR","DS","DT","DV","DW","DY",
                    "EA","EC","ED","EF","EG","EH","EI","EK","EL","EM","EN","EP","EQ","ER","ES","ET","EV","EW","EY",
                    "FA","FC","FD","FE","FG","FH","FI","FK","FL","FM","FN","FP","FQ","FR","FS","FT","FV","FW","FY",
                    "GA","GC","GD","GE","GF","GH","GI","GK","GL","GM","GN","GP","GQ","GR","GS","GT","GV","GW","GY",
                    "HA","HC","HD","HE","HF","HG","HI","HK","HL","HM","HN","HP","HQ","HR","HS","HT","HV","HW","HY",
                    "IA","IC","ID","IE","IF","IG","IH","IK","IL","IM","IN","IP","IQ","IR","IS","IT","IV","IW","IY",
                    "KA","KC","KD","KE","KF","KG","KH","KI","KL","KM","KN","KP","KQ","KR","KS","KT","KV","KW","KY",
                    "LA","LC","LD","LE","LF","LG","LH","LI","LK","LM","LN","LP","LQ","LR","LS","LT","LV","LW","LY",
                    "MA","MC","MD","ME","MF","MG","MH","MI","MK","ML","MN","MP","MQ","MR","MS","MT","MV","MW","MY",
                    "NA","NC","ND","NE","NF","NG","NH","NI","NK","NL","NM","NP","NQ","NR","NS","NT","NV","NW","NY",
                    "PA","PC","PD","PE","PF","PG","PH","PI","PK","PL","PN","PM","PQ","PR","PS","PT","PV","PW","PY",
                    "QA","QC","QD","QE","QF","QG","QH","QI","QK","QL","QN","QM","QP","QR","QS","QT","QV","QW","QY",
                    "RA","RC","QD","RE","RF","RG","RH","RI","RK","RL","RN","RM","RP","RQ","RS","RT","RV","RW","RY",
                    "SA","SC","SD","SE","SF","SG","SH","SI","SK","SL","SM","SN","SP","SQ","SR","ST","SV","SW","SY",
                    "TA","TC","TD","TE","TF","TG","TH","TI","TK","TL","TM","TN","TP","TQ","TR","TS","TV","TW","TY",
                    "VA","VC","VD","VE","VF","VG","VH","VI","VK","VL","VM","VN","VP","VQ","VR","VS","VT","VW","VY",
                    "WA","WC","WD","WE","WF","WG","WH","WI","WK","WL","WM","WN","WP","WQ","WR","WS","WT","WV","WY",
                    "YA","YC","YD","YE","YF","YG","YH","YI","YK","YL","YM","YN","YP","YQ","YR","YS","YT","YV","YW")

  absent_labels <- base::setdiff(duplet_labels, unique(data[,"Duplets"]))

  if (length(absent_labels) > 0){
    print("Here is a list of missing duplets:")
    print(absent_labels)
  }else{
    print("All 380 possible duplets are formed")
  }

  return(data)

}
