
## Load functions for parsing the flow cytometry spreadsheets
## The gdata library is necessary for reading Excel spreadsheets; it will be loaded as well.
source('flow_data_cleaning_functions.r')

## View help documentation on the functions listed above
describe(read_flow_exp_file)

## Load all expected flow variables (expected column names)
flow_cn = read.delim('./data/final_flow.txt', sep='\t', as.is=T, header=F)
flow_cn = flow_cn[,1]

## Move to the directory holding the data
flow_dir = "~/Documents/MyDocuments/SystemsImmunogenetics/WNV/Lund_Flow_fixed_Nov_13"
setwd(flow_dir)

## Get a list of data files to read (in this case all flow spreadsheets begin with the prefix 'Expt')
flow_files = list.files('.', pattern="Expt.*\\.xls")

print(flow_files)

## Iterate through all the files, parse each, and merge all data into a single dataframe
i = 1
for (file in flow_files) {
    print(file)
    flow_dat = read_flow_exp_file(file, flow_cn)
    
    ## Check if there are any unexpected columns
    new_columns = setdiff(colnames(flow_dat), flow_cn)
    if (length(new_columns) > 0) {
        flow_cn = c(flow_cn, new_columns)
    }
    if (i > 1) {
        ## Fill extra columns with NAs
        for (col in new_columns) {
            flow_all[,col] = NA
        }
        ## Merge data
        flow_all = rbind(flow_all[,flow_cn], flow_dat[,flow_cn])
    } else {
        flow_all = flow_dat
    }
    i = i + 1
}

## Check the dimensions of the dataframe
dim(flow_all)

## Check that all expected columns are present
setdiff(flow_cn, colnames(flow_all))

setdiff(colnames(flow_all), flow_cn)

## Order columns, add Lab column and fix formatting
flow_all = flow_all[, flow_cn]
flow_all$Lab = "Lund"

flow_all$ID = gsub(" ", "", flow_all$ID)
flow_all$ID = gsub("X", "x", flow_all$ID)
flow_all$Mating = gsub(" ", "", flow_all$Mating)
flow_all$Mating = gsub("X", "x", flow_all$Mating) 
flow_all$UW_Line = as.numeric(flow_all$UW_Line)

## Check for duplicate IDs
new_flow_ids = paste(flow_all$ID, flow_all$Tissue, sep='_')
sum(duplicated(new_flow_ids))

## Change all data columns to numeric
for (i in 11:277) {
    flow_all[,i] = as.numeric(flow_all[,i])
}

## Calculate cell counts and ratios
flow_full = flow_all
flow_full = calc_treg_counts(flow_full)
flow_full = calc_tcell_counts(flow_full)
flow_full = calc_ics_counts(flow_full)
flow_full = calc_ics_percent_ratios(flow_full)
flow_full = calc_ics_count_ratios(flow_full)
flow_full = clean_inf_nan(flow_full)

## Save R data file
write.table(flow_full, file='Lund_Flow_Full_29-Jan-2016_final.txt', col.names=T, row.names=F, quote=F, sep='\t', na='')
save(flow_full, file='lund_flow_full_29-jan-2016_final.rda')
