describe = function(obj) {
  if ('help' %in% names(attributes(obj))) {
    writeLines(attr(obj, 'help'))
  }
}
attr(describe, 'help') = "
This function prints the contents of the 'help' attribute of any R object. 
It is meant to provide help documentation in the same vein as Docstrings in Python. 
"

read_flow_exp_file = function(f, cn_expected) {
  ### Parse the Treg panel - first sheet of the excel workbook
  ## Read the file
  dat1 = read.xls(f, sheet=1, verbose=F, header=F, as.is=T, na.strings=c("NA", "", " ", "#DIV/0!"))
  
  ## Identify the first row containing information
  start_row = min(which(!is.na(dat1[,1])))
  ## Get the column names from the first row
  cn1 = dat1[start_row,]
  ## Subset the dataframe and list of column names to remove empty columns
  dat1 = dat1[,!is.na(cn1)]
  cn1 = cn1[!is.na(cn1)]
  
  ## DEBUGGING / TESTING
  #print(cn1)
  
  ## Parse nested column names
  cn_start = which(!is.na(dat1[1,]))
  lenTreg = cn_start[2]-cn_start[1]-1
  cn1[cn_start[1]:(cn_start[1]+lenTreg-1)] = paste(rep('CD4pos_Foxp3neg_', lenTreg), cn1[cn_start[1]:(cn_start[1]+lenTreg-1)], sep='')
  cn1[cn_start[2]:(cn_start[2]+lenTreg-1)] = paste(rep('Tregs_', lenTreg), cn1[cn_start[2]:(cn_start[2]+lenTreg-1)], sep='')
  
  ## Identify the time column
  start_col = which(cn1=='time')
  ## Standardize the column names
  cn1 = fix_column_names(cn1, 'treg', start_col)
  ## Identify the last row with data
  end_row = which(is.na(dat1[,1]))[2]
  ## Update the column names and subset the dataframe
  colnames(dat1) = cn1
  dat1 = dat1[(start_row+1):(end_row-1),]
  
  ## DEBUGGING / TESTING
  #print(end_row)
  #print(dim(dat1))
  
  ## Parse the T-cell panel - second and third sheets in the excel workbook
  ## Initialize counter to keep track of what sheet is being parsed
  i = 1
  ## Identify the T-cell sheets (since both D7 and D21 panels are not always included)
  CD8_sheets = grep("CD8", sheetNames(f))
  print(sheetNames(f)[CD8_sheets])
  ## Iterate through the sheets containing the T-cell panels
  for (sheet in CD8_sheets) {
    if (i == 1) {
      sheet_name = sheetNames(f)[sheet]
      ## Read the file, identify the first row containing information, and get column names
      dat2 = read.xls(f, sheet=sheet, verbose=F, header=F, as.is=T, na.strings=c("NA", "", " ", "#DIV/0!"))
      start_row = min(which(dat2[,1]=='ID' | dat2[,1]=='UNC strain'))
      cn2 = dat2[start_row,]
      dat2 = dat2[,!is.na(cn2)]
      cn2 = cn2[!is.na(cn2)]
      
      ## DEBUGGING / TESTING
      #print(cn2)
      
      ## Parse nested column names
      cn_start = which(!is.na(dat2[1,]))
      lenCD8 = cn_start[2]-cn_start[1]-1
      lenCD4 = length(cn2)-cn_start[2]+1
      cn2[cn_start[1]:(cn_start[1]+lenCD8-1)] = paste(rep('CD8pos_', lenCD8), cn2[cn_start[1]:(cn_start[1]+lenCD8-1)], sep='')
      cn2[cn_start[2]:(cn_start[2]+lenCD4-1)] = paste(rep('CD4pos_', lenCD4), cn2[cn_start[2]:(cn_start[2]+lenCD4-1)], sep='')
      
      ## Identify the time column
      start_col = which(cn2=='time')
      ## Determine which T-cell panel is being parsed (D7 or D21) and standardize the column names
      if (length(grep("d7", sheet_name)) > 0 | length(grep("d12", sheet_name)) > 0) {
        cn2 = fix_column_names(cn2, 'tcell_d7', start_col)
      } else {
        cn2 = fix_column_names(cn2, 'tcell_d21', start_col)
      }
      ## Identify the last row with data
      end_row = which(is.na(dat2[,1]))[2]
      ## Update the column names and subset the dataframe
      colnames(dat2) = cn2
      dat2 = dat2[(start_row+1):(end_row-1),]
    } else {
      ## Repeat the above process for the second T-cell panel if it exists in the file
      sheet_name = sheetNames(f)[sheet]
      dat2b = read.xls(f, sheet=sheet, verbose=F, header=F, as.is=T, na.strings=c("NA", "", " ", "#DIV/0!"))
      start_row = min(which(dat2b[,1]=='ID' | dat2b[,1]=='UNC strain'))
      cn2 = dat2b[start_row,]
      dat2b = dat2b[,!is.na(cn2)]
      cn2 = cn2[!is.na(cn2)]
      
      ## DEBUGGING / TESTING
      #print(cn2)
      
      cn_start = which(!is.na(dat2b[1,]))
      lenCD8 = cn_start[2]-cn_start[1]-1
      lenCD4 = length(cn2)-cn_start[2]+1
      cn2[cn_start[1]:(cn_start[1]+lenCD8-1)] = paste(rep('CD8pos_', lenCD8), cn2[cn_start[1]:(cn_start[1]+lenCD8-1)], sep='')
      cn2[cn_start[2]:(cn_start[2]+lenCD4-1)] = paste(rep('CD4pos_', lenCD4), cn2[cn_start[2]:(cn_start[2]+lenCD4-1)], sep='')
      start_col = which(cn2=='time')
      if (length(grep("d7", sheet_name)) > 0 | length(grep("d12", sheet_name)) > 0) {
        cn2 = fix_column_names(cn2, 'tcell_d7', start_col)
      } else {
        cn2 = fix_column_names(cn2, 'tcell_d21', start_col)
      }
      end_row = which(is.na(dat2b[,1]))[2]
      colnames(dat2b) = cn2
      dat2b = dat2b[(start_row+1):(end_row-1),]
      
      ## Merge the two T-cell panels, since only some animals (timepoints) will be included in each
      dat2 = merge(dat2, dat2b, by=intersect(colnames(dat2), colnames(dat2b)), all=T)
    }
    i = i + 1
  }
  
  ## DEBUGGING / TESTING
  #print(end_row)
  #print(dim(dat2))
  
  ## Parse the ICS Panel - the third (or fourth) sheet in the excel workbook
  ## Read the file, identify the first row containing information, and get column names
  dat3 = read.xls(f, sheet=(max(CD8_sheets)+1), verbose=F, header=F, as.is=T, na.strings=c("NA", "", " ", "#DIV/0!"))
  start_row = min(which(dat3[,1]=='ID' | dat3[,1]=='UNC strain'))
  cn3 = dat3[start_row,]
  dat3 = dat3[,!is.na(cn3)]
  cn3 = cn3[!is.na(cn3)]
  
  ## Parse the nested column names
  cn_start1 = which(!is.na(dat3[1,]))
  cn_start2 = which(!is.na(dat3[2,]))
  lenICSCD8_1 = cn_start2[2]-cn_start2[1]-1
  lenICSCD4_1 = cn_start2[3]-cn_start2[2]-6
  lenICSCD8_2 = cn_start2[4]-cn_start2[3]-1
  lenICSCD4_2 = cn_start2[5]-cn_start2[4]-6
  lenICSCD8_3 = cn_start2[6]-cn_start2[5]-1
  lenICSCD4_3 = cn_start2[7]-cn_start2[6]-6
  lenICSCD8_4 = cn_start2[8]-cn_start2[7]-1
  lenICSCD4_4 = length(cn3)-cn_start2[8]+1
  cn3[cn_start2[1]:(cn_start2[1]+lenICSCD8_1-1)] = paste(rep('CD8pos_', lenICSCD8_1), cn3[cn_start2[1]:(cn_start2[1]+lenICSCD8_1-1)], sep='')
  cn3[cn_start2[2]:(cn_start2[2]+lenICSCD4_1-1)] = paste(rep('CD4pos_', lenICSCD4_1), cn3[cn_start2[2]:(cn_start2[2]+lenICSCD4_1-1)], sep='')
  cn3[cn_start2[3]:(cn_start2[3]+lenICSCD8_2-1)] = paste(rep('CD8_', lenICSCD8_2), cn3[cn_start2[3]:(cn_start2[3]+lenICSCD8_2-1)], sep='')
  cn3[cn_start2[4]:(cn_start2[4]+lenICSCD4_2-1)] = paste(rep('CD4_', lenICSCD4_2), cn3[cn_start2[4]:(cn_start2[4]+lenICSCD4_2-1)], sep='')
  cn3[cn_start2[5]:(cn_start2[5]+lenICSCD8_3-1)] = paste(rep('CD8_', lenICSCD8_3), cn3[cn_start2[5]:(cn_start2[5]+lenICSCD8_3-1)], sep='')
  cn3[cn_start2[6]:(cn_start2[6]+lenICSCD4_3-1)] = paste(rep('CD4_', lenICSCD4_3), cn3[cn_start2[6]:(cn_start2[6]+lenICSCD4_3-1)], sep='')
  cn3[cn_start2[7]:(cn_start2[7]+lenICSCD8_4-1)] = paste(rep('CD8_', lenICSCD8_4), cn3[cn_start2[7]:(cn_start2[7]+lenICSCD8_4-1)], sep='')
  cn3[cn_start2[8]:(cn_start2[8]+lenICSCD4_4-1)] = paste(rep('CD4_', lenICSCD4_4), cn3[cn_start2[8]:(cn_start2[8]+lenICSCD4_4-1)], sep='')
  
  lenDMSO = cn_start1[2]-cn_start1[1]
  lenNS4B = cn_start1[3]-cn_start1[2]
  lenHIWNV = cn_start1[4]-cn_start1[3]
  lenCD3CD28 = length(cn3)-cn_start1[4]+1
  cn3[cn_start1[1]:(cn_start1[1]+lenDMSO-1)] = paste(rep('DMSO_', lenDMSO), cn3[cn_start1[1]:(cn_start1[1]+lenDMSO-1)], sep='')
  cn3[cn_start1[2]:(cn_start1[2]+lenNS4B-1)] = paste(rep('NS4B_', lenNS4B), cn3[cn_start1[2]:(cn_start1[2]+lenNS4B-1)], sep='')
  cn3[cn_start1[3]:(cn_start1[3]+lenHIWNV-1)] = paste(rep('HIWNV_', lenHIWNV), cn3[cn_start1[3]:(cn_start1[3]+lenHIWNV-1)], sep='')
  cn3[cn_start1[4]:(cn_start1[4]+lenCD3CD28-1)] = paste(rep('CD3CD28_', lenCD3CD28), cn3[cn_start1[4]:(cn_start1[4]+lenCD3CD28-1)], sep='')
  
  ## Identify the first time column (DMSO ICS panel)
  start_col = which(cn3=='DMSO_time')
  ## Standardize the column names
  cn3 = fix_column_names(cn3, 'ics', start_col)
  ## Identify the last column with data
  end_row = which(is.na(dat3[,1]))[3]
  ## Update the column names and subset the dataframe
  colnames(dat3) = cn3
  dat3 = dat3[(start_row+1):(end_row-1),]
  
  ## DEBUGGING / TESTING
  #print(end_row)
  #print(dim(dat3))
  #print(colnames(dat1))
  #print(colnames(dat2))
  #print(colnames(dat3))
  
  ## Merge all the panels
  dat = merge(dat1, dat2 , by=c('UNC_strain', 'UW_strain', 'RIX_ID', 'Timepoint', 'Tissue', 'Total_Cell_Count'), suffixes=c('', ''))
  dat = merge(dat, dat3, by=c('UNC_strain', 'UW_strain', 'RIX_ID', 'Timepoint', 'Tissue', 'Total_Cell_Count'), suffixes=c('',''))
  
  ## Fix final column names
  cn_final = colnames(dat)
  cn_final = gsub(" ", "_", cn_final)
  cn_final = gsub("UNC_strain", "Mating", cn_final)
  cn_final = gsub("UW_strain", "UW_Line", cn_final)
  colnames(dat) = cn_final
  
  ## Update the ID column
  ids = paste(dat$Mating, dat$RIX_ID, sep='_')
  dat$ID = ids
  
  ## Fix timepoints
  dat$Timepoint[dat$Timepoint=='d7'] = '7'
  dat$Timepoint[dat$Timepoint=='d12'] = '12'
  dat$Timepoint[dat$Timepoint=='d21'] = '21'
  dat$Timepoint[dat$Timepoint=='d28'] = '28'
  dat$Timepoint[dat$Timepoint=='d12m'] = '12m'
  dat$Timepoint[dat$Timepoint=='d28m'] = '28m'
  
  ## Create Virus column
  dat$Virus = 'WNV'
  dat$Virus[dat$Timepoint=='12m'] = 'Mock'
  dat$Virus[dat$Timepoint=='28m'] = 'Mock'
  dat$Timepoint[dat$Timepoint=='12m'] = '12'
  dat$Timepoint[dat$Timepoint=='28m'] = '28'
  
  ## Print unexpected column names, and add columns with NAs if necessary
  if (length(setdiff(colnames(dat), cn_expected)) > 0) {
    print("Unexpected columns: ")
    print(setdiff(colnames(dat), cn_expected))
  }
  for (col in cn_expected) {
    if (!(col %in% colnames(dat))) {
      dat[,col] = NA
    }
  }
  return(dat)
}
attr(read_flow_exp_file, 'help') = "
This function parses flow cytometry data from an Excel workbook.

Parameters
f: The Excel file name.
cn_expected: A character vector containing the expected column names.

Returns
A dataframe containing the processed flow cytometry data.
"

fix_column_names = function(cn, panel, start_col=8) {
  ## DEBUGGING / TESTING
  #print(start_col)
  #print(cn)
  
  ## Standardize the column names for all panels
  cn = trim(gsub("\\+([^[:space:]])", "+ \\1", cn))
  cn = gsub("\\+", "pos", cn)
  cn = gsub("IL-17", "IL_17", cn)
  cn = gsub("CTLA-4", "CTLA_4", cn)
  cn = gsub("PSGL-1", "PSGL_1", cn)
  cn = trim(gsub("-([^[:space:]])", "- \\1", cn))
  cn = gsub(" ", "_", cn)
  cn = gsub("-", "neg", cn)
  cn = gsub("cell", "Cell", cn)
  cn = gsub("count", "Count", cn)
  cn = gsub("live", "Live", cn)
  cn = gsub("singlets", "Singlets", cn)
  cn = gsub("lymphocytes", "Lymphocytes", cn)
  cn = gsub("total", "Total", cn)
  cn = gsub("time", "Time", cn)
  cn = gsub("Organ", "Tissue", cn)
  cn = gsub("Mouse_#", "RIX_ID", cn)
  cn = gsub("NS4B_CD4pos_Tbetpos", "NS4B_CD4_Tbetpos", cn)
  cn = gsub("HIWNV_CD4pos_Tbetpos", "HIWNV_CD4_Tbetpos", cn)
  cn = gsub("CD3CD28_CD4pos_Tbetpos", "CD3CD28_CD4_Tbetpos", cn)
  cn = gsub("%_CTLA_4pos", "CTLA_4pos", cn)
  cn = gsub("NS4b", "NS4B", cn)
  cn = gsub("TNFa", "TNFA", cn)
  cn = gsub("IFNg", "IFNG", cn)
  
  panel = paste(panel, '_', sep='')
  cn[start_col:length(cn)] = gsub("^", panel, cn[start_col:length(cn)])
  
  ## DEBUGGING / TESTING
  #print(cn)
  
  return(cn)
}
attr(fix_column_names, 'help') = "
This function standardizes the column names of each flow cytometry panel.

Parameters
cn: A character vector containing the original column names.
panel: A string which will be appended to column names indicating the panel being processed (e.g. 'treg').
start_column: A number indicating in which column the flow data starts (preceding columns are sample identifiers) 

Returns
A character vector containing the fixed column names.
"