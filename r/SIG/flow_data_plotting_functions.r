library(colorRamps)
library(RColorBrewer)
library(pheatmap)
#options(warn=-1)

## Help Documentation
describe = function(obj) {
  if ('help' %in% names(attributes(obj))) {
    writeLines(attr(obj, 'help'))
  }
}
attr(describe, 'help') = "
This function prints the contents of the 'help' attribute of any R object. 
It is meant to provide help documentation in the same vein as Docstrings in Python. 
"

## Format the data for boxplots
flow_boxplot_data = function(flow_df, lines, tissue, flow_vars, tp, line_colors=NULL, mocks_only=FALSE, data_type=1) {
	## Initialize lists
	bdat = list()
	pdat = list()
	
	## Get Mock data
	n_mocks = 0
	mock_lines = c()
	for (l in lines) {
		label = paste(l, 'MOCK', sep='_')
		d = flow_df[flow_df$UW_Line==l & tolower(flow_df$Virus)=='mock' & flow_df$Tissue==tissue & flow_df$Timepoint %in% as.numeric(tp), flow_vars]
		if (!all(is.na(d))) {
			bdat[[label]] = d
			n_mocks = n_mocks + 1
			mock_lines = c(mock_lines, l)
		}
		tmp_list = list()
		for (t in tp) {
			d = flow_df[flow_df$UW_Line==l & tolower(flow_df$Virus)=='mock' & flow_df$Tissue==tissue & flow_df$Timepoint == as.numeric(t), flow_vars]
			if (!all(is.na(d))) {
				tmp_list[[t]] = d
			}
		}
		if (length(tmp_list) > 0) {
			pdat[[label]] = tmp_list
		}
	}
	
	## Get WNV data
	n_wnv = 0
	wnv_lines = c()
	if (!mocks_only) {
		for (l in lines) {
			label = paste(l, 'WNV', sep='_')
			d = flow_df[flow_df$UW_Line==l & tolower(flow_df$Virus)=='wnv' & flow_df$Tissue==tissue & flow_df$Timepoint %in% as.numeric(tp), flow_vars]
			if (!all(is.na(d))) {
				bdat[[label]] = d
				n_wnv = n_wnv + 1
				wnv_lines = c(wnv_lines, l)
			}
			tmp_list = list()
			for (t in tp) {
				d = flow_df[flow_df$UW_Line==l & tolower(flow_df$Virus)=='wnv' & flow_df$Tissue==tissue & flow_df$Timepoint == as.numeric(t), flow_vars]
				if (!all(is.na(d))) {
					tmp_list[[t]] = d
				}
			}
			if (length(tmp_list) > 0) {
				pdat[[label]] = tmp_list
			}
		}
	}
	
	## Create boxplot colors
	union_lines = union(mock_lines, wnv_lines)
	n_colors = length(union_lines)
	l_colors = colorRampPalette(brewer.pal(9,'Set1'))(n_colors)
	if (!mocks_only) {
		mock_colors = c()
		for (l in mock_lines) {
			mock_colors = c(mock_colors, l_colors[which(union_lines==l)])
		}
		wnv_colors = c()
		for (l in wnv_lines) {
			wnv_colors = c(wnv_colors, l_colors[which(union_lines==l)])
		}
		line_colors = c(mock_colors, wnv_colors)
	} else {
		line_colors = c(l_colors[1:n_mocks])
	}
	
	## DEBUGGING / TESTING
	#print(line_colors)
	
	## Determine data type (e.g. percentages vs. absolute counts)
	data_type = as.numeric(data_type)
	if (length(grep("^ics_.*", flow_vars)) > 0) {
		data_type = data_type + 2
	}
	
	## Create plot title with heritability annotation
	#icc_var1 = paste(tissue, '_mock_icc_gelman_hill', sep="")
	#icc_var2 = paste(tissue, '_icc_gelman_hill', sep="")
	#title = paste(flow_vars, " (ICC - Mock = ", round(flow_heritability[flow_vars, icc_var1], digits=3), "; ICC - All = ", round(flow_heritability[flow_vars, icc_var2], digits=3), ")", sep="")
	title = flow_vars
	
	## Return heatmap matrix
	return(list(bdat_list=bdat, pdat_list=pdat, colors=line_colors, num_mocks=n_mocks, time_points=tp, dtype=data_type, title=title))	
}
attr(flow_boxplot_data, 'help') = "
This function aggregates the data needed for creating boxplots of the flow cytometry data.

Parameters:
flow_df: The dataframe containing the flow data.
lines: A numeric vector containing the mouse lines that should be plotted.
tissue: A string indicating the tissue (e.g. 'brain' or 'spleen')
flow_vars: A string indicating the flow variable to be plotted.
tp: A character vector containing the time points to be included.
line_colors: A vector of colors (default=NULL; colors will be determined automatically).
mocks_only: Logical value indicating whether mocks only will be plotted.
data_type: A number indicating whether percentages (1) or absolute cell counts (2) will be plotted.

Returns:
A list containing the data to be plotted; input for the flow_boxplots() function.
"

## Boxplots
flow_boxplots = function(boxplot_data, ...) {
	tp_symbols = c("7"=21, "12"=22, "21"=23, "28"=24)
	if (boxplot_data$dtype == 1) {
		ylab="Percent"
	} else if (boxplot_data$dtype == 2) {
		ylab="Cell Count"
	} else if (boxplot_data$dtype == 3) {
		ylab="Percent Ratio"
	} else {
		ylab="Cell Count Ratio"
	}
	
	if (is.na(boxplot_data$y_min) | is.na(boxplot_data$y_max)) {
		if (boxplot_data$dtype == 1) {
			ylim = c(0,105)
		} else {
			ylim = NULL
		}
	} else {
		ylim = c(boxplot_data$y_min, boxplot_data$y_max)
	}
	
	if (boxplot_data$rm_outliers) {
		if (boxplot_data$dtype == 1) {
			boxplot(boxplot_data$bdat_list, col=boxplot_data$colors, xaxt = "n",  xlab = "", ylab=ylab, main=boxplot_data$title, outline=FALSE, ylim=ylim)
			if (boxplot_data$show_data) {
				for (t in boxplot_data$time_points) {
					d = lapply(boxplot_data$pdat_list, function(x){x[[t]]})
					stripchart(d, method='jitter', pch=tp_symbols[t], lwd=2, cex=1.25, vertical=T, add=T)
				}
			}
		} else {
			boxplot(boxplot_data$bdat_list, col=boxplot_data$colors, xaxt = "n",  xlab = "", ylab=ylab, main=boxplot_data$title, outline=FALSE, ylim=ylim)
			if (boxplot_data$show_data) {
				for (t in boxplot_data$time_points) {
					d = lapply(boxplot_data$pdat_list, function(x){x[[t]]})
					stripchart(d, method='jitter', pch=tp_symbols[t], lwd=2, cex=1.25, vertical=T, add=T)
				}
			}
		}
	} else {
		if (boxplot_data$dtype == 1) {
			boxplot(boxplot_data$bdat_list, col=boxplot_data$colors, xaxt = "n",  xlab = "", ylab=ylab, main=boxplot_data$title, pch=8, cex=1.5, ylim=ylim)
			if (boxplot_data$show_data) {
				for (t in boxplot_data$time_points) {
					d = lapply(boxplot_data$pdat_list, function(x){x[[t]]})
					stripchart(d, method='jitter', pch=tp_symbols[t], lwd=2, cex=1.25, vertical=T, add=T)
				}
			}
		} else {
			boxplot(boxplot_data$bdat_list, col=boxplot_data$colors, xaxt = "n",  xlab = "", ylab=ylab, main=boxplot_data$title, pch=8, cex=1.5, ylim=ylim)
			if (boxplot_data$show_data) {
				for (t in boxplot_data$time_points) {
					d = lapply(boxplot_data$pdat_list, function(x){x[[t]]})
					stripchart(d, method='jitter', pch=tp_symbols[t], lwd=2, cex=1.25, vertical=T, add=T)
				}
			}
		}
	}
	if (boxplot_data$show_data) {
		legend('topleft', legend=c('D7', 'D12', 'D21', 'D28'), pch=tp_symbols, pt.lwd=2, pt.cex=1.25, ncol=4)
	}
	labels = names(boxplot_data$bdat_list)
	axis(1, at=1:length(labels), labels=FALSE)
	y_interval = (par("usr")[4] - par("usr")[3])/675
	y_label_pos = par("usr")[3]-(80*y_interval)
	text(x=seq_along(labels), y=y_label_pos, srt=270, adj=1, labels=labels, xpd=TRUE, ...)
	abline(v=boxplot_data$num_mocks+0.5, lty=3)
}
attr(flow_boxplots, 'help') = "
This function plots the data returned by flow_boxplot_data().

Parameters:
boxplot_data: This should be a list combining the value returned by flow_boxplot_data() and additional options.

Returns:
NULL

Example:
boxplot_data = flow_boxplot_data(flow_full, c(7,8,9), 'brain', 'treg_T_regs', c('7','12','21','28'))
opts = list(rm_outliers=F, show_data=F, y_min=0, y_max=100)
flow_boxplots(c(boxplot_data, opts))
"