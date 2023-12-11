#!/usr/bin/env	python3

"""
Obtain metadata summary reports
By Veronica Mixao
@INSA
"""


import os
import sys
import argparse
import textwrap
import pandas
from pandas.api.types import is_datetime64_any_dtype as is_datetime
import datetime

version = "1.1.0"
last_updated = "2023-12-11"

# functions	----------

def partitions2metadata(partitions_name, partitions, mx_metadata, partitions2report, filters, log):
	""" 
	create a combined matrix with metadata and partitions information
	
	input:
	partitions = partitions *.tsv file
	metadata = metadata *.tsv file
	partitions2report = list of partition columns to include in the final matrix
	filter = list of filters
	log = log file already opened 
	
	returns:
	dataframe with the required metadata and partitions
	"""

	columns_names = [col.replace(" ", "_") for col in mx_metadata.columns]
	#mx_metadata.replace(" ", "_", regex=True, inplace=True)
	mx_metadata.columns = columns_names
	sample_column = mx_metadata.columns[0]
	possible_subset = False

	# checking for 'date' column
	if "date" in mx_metadata.columns:
		index_no = mx_metadata.columns.get_loc("date")
		mx_metadata["date_original"] = mx_metadata["date"]
		date_original = mx_metadata.pop("date_original")
		mx_metadata.insert(index_no, "date_original", date_original)
		index_no = mx_metadata.columns.get_loc("date")
		if "year" in mx_metadata.columns:
			mx_metadata["year_original"] = mx_metadata["year"]
			year_original = mx_metadata.pop("year_original")
			mx_metadata.insert(index_no, "year_original", year_original)
			index_no = mx_metadata.columns.get_loc("date")
		mx_metadata["date"] = pandas.to_datetime(mx_metadata["date"], errors = "coerce")
		mx_metadata["year"] = mx_metadata["date"].dt.year
		year = mx_metadata.pop("year")
		mx_metadata.insert(index_no + 1, "year", year)
		index_no = mx_metadata.columns.get_loc("date")
		if "iso_week_nr" not in mx_metadata.columns and "iso_year" not in mx_metadata.columns and "iso_week" not in mx_metadata.columns:
			isoyear = mx_metadata["date"].dt.isocalendar().year
			isoweek = mx_metadata["date"].dt.isocalendar().week
			mx_metadata["iso_year"] = isoyear.astype(str)
			mx_metadata["iso_week_nr"] = isoweek.astype(str)
			mx_metadata["iso_week"] = isoyear.astype(str).replace("<NA>", "-") + "W" + isoweek.astype(str).replace("<NA>", "--").apply(lambda x: x.zfill(2))
			isoyear = mx_metadata.pop("iso_year")
			isoweek = mx_metadata.pop("iso_week_nr")
			isodate = mx_metadata.pop("iso_week")
			mx_metadata.insert(index_no + 2, "iso_year", isoyear)
			mx_metadata.insert(index_no + 3, "iso_week_nr", isoweek)
			mx_metadata.insert(index_no + 4, "iso_week", isodate)

	# check for duplicated samples in metadata
	metadata_samples = mx_metadata[mx_metadata.columns[0]].values.tolist()
	
	if len(metadata_samples) != len(set(metadata_samples)):
		print("\tWARNING!! You have duplicated samples in the metadata table! I cannot continue!")
		print("\tWARNING!! You have duplicated samples in the metadata table! I cannot continue!", file = log)
		sys.exit()
	
	if isinstance(partitions, pandas.DataFrame):
		mx_partitions = pandas.DataFrame(data = partitions, dtype = str)
		sample_column_part = mx_partitions.columns[0]

		# check for duplicated samples in partitions
		partitions_samples = mx_partitions[mx_partitions.columns[0]].values.tolist()
		if len(partitions_samples) != len(set(partitions_samples)):
			print("\tWARNING!! You have duplicated columns in the partitions table! I cannot continue!")
			print("\tWARNING!! You have duplicated columns in the partitions table! I cannot continue!", file = log)
			sys.exit()
		
		a = mx_metadata.set_index(sample_column, drop = True)
		b = mx_partitions.set_index(sample_column_part, drop = True)
		
		if len(set(a.columns) & set(b.columns)) > 0:
			b = b.add_prefix("subset_")
			possible_subset = True
			mx_partitions.columns = [sample_column_part] + b.columns.tolist()
			mx_partitions.to_csv(partitions_name, index = False, header=True, sep ="\t")
		if partitions2report == "all": # add all partitions
			c = pandas.concat([a, b], axis=1)
			
		else: # add specific set of partitions
			required_partitions = partitions2report.split(",")
			c = a
			for column_name in required_partitions:
				if column_name in b.columns:
					c = pandas.concat([c,b[column_name]], axis = 1)
				else:
					print("\t\t" + column_name + " will not be reported because it was not found in the partitions table!!")
					print("\t\t" + column_name + " will not be reported because it was not found in the partitions table!!", file = log)
		
	else:
		c = mx_metadata.set_index(mx_metadata.columns[0], drop = True)
		
		
	# filtering table according to specifications
	c_filtered = c.reset_index(drop = False)
	c_filtered.rename(columns={c_filtered.columns[0]: mx_metadata.columns[0]}, inplace=True)

	if filters != "":
		print("\tFiltering metadata for the following parameters: " + " & ".join(filters.split(";")))
		print("\tFiltering metadata for the following parameters: " + " & ".join(filters.split(";")), file = log)
		f = []
		if ";" in filters:
			for flt in filters.split(";"):
				f.append(flt)
		else:
			f.append(filters)

		for spec in f:
			col = spec.split(" ", 2)[0]
			val = spec.split(" ", 2)[1]
			cond = " ".join(spec.split(" ")[2:])
			
			if "," in cond:
				lst = cond.split(",")
				if val == "==":
					c_filtered = c_filtered[c_filtered[col].isin(lst)]
				elif val == "!=":
					c_filtered = c_filtered[c_filtered[col].isin(lst) == False]
			else:
				if col == "date":
					c_filtered["date"] = c_filtered["date"].astype("datetime64[ns]")
					if val == "==":
						c_filtered = c_filtered[c_filtered["date"] == cond]  
					elif val == "!=":
						c_filtered = c_filtered[c_filtered["date"] != cond] 
					elif val == ">":
						c_filtered = c_filtered[c_filtered["date"] > cond] 
					elif val == ">=":
						c_filtered = c_filtered[c_filtered["date"] >= cond] 
					elif val == "<=":
						c_filtered = c_filtered[c_filtered["date"] <= cond] 
					elif val == "<":
						c_filtered = c_filtered[c_filtered["date"] < cond]
				elif col == "iso_week":
					if "date" in mx_metadata.columns:
						year = cond.split("W")[0]
						week = cond.split("W")[1]
						cond = pandas.to_datetime(date.fromisocalendar(int(year), int(week), 1))
						c_filtered["date"] = c_filtered["date"].astype("datetime64[ns]")
						if val == "==":
							c_filtered = c_filtered[c_filtered["date"] == cond]  
						elif val == "!=":
							c_filtered = c_filtered[c_filtered["date"] != cond] 
						elif val == ">":
							c_filtered = c_filtered[c_filtered["date"] > cond] 
						elif val == ">=":
							c_filtered = c_filtered[c_filtered["date"] >= cond] 
						elif val == "<=":
							c_filtered = c_filtered[c_filtered["date"] <= cond] 
						elif val == "<":
							c_filtered = c_filtered[c_filtered["date"] < cond]		
					else:
						print("\tCannot apply the 'iso_week' filter because column 'date' was not found in the metadata!")
				else:
					if val == "==":
						c_filtered = c_filtered[c_filtered[col] == cond]
					elif val == "!=":
						c_filtered = c_filtered[c_filtered[col] != cond]
					else:
						c_filtered[col] = pandas.to_numeric(c_filtered[col], errors='coerce')
						if val == ">":
							c_filtered = c_filtered[c_filtered[col] > float(cond)]
						elif val == ">=":
							c_filtered = c_filtered[c_filtered[col] >= float(cond)]
						elif val == "<":
							c_filtered = c_filtered[c_filtered[col] < float(cond)]
						elif val == "<=":
							c_filtered = c_filtered[c_filtered[col] >= float(cond)]
						c_filtered[col] = c_filtered[col].astype(int)
						c_filtered[col] = c_filtered[col].astype(str)

	new_metadata = c_filtered
	
	# check for missing samples
	if isinstance(partitions, pandas.DataFrame):
		missing_in_metadata = set(partitions_samples) - set(metadata_samples)
		missing_in_partitions = set(metadata_samples) - set(partitions_samples)
		
		print("\t\tSamples present in partitions table but missing in metadata table: " + ",".join(list(missing_in_metadata)))
		print("\t\tSamples present in partitions table but missing in metadata table: " + ",".join(list(missing_in_metadata)), file = log)
		print("\t\tSamples not present in partitions table but present in metadata table: " + ",".join(list(missing_in_partitions)))
		print("\t\tSamples not present in partitions table but present in metadata table: " + ",".join(list(missing_in_partitions)), file = log)
		
	return new_metadata, possible_subset


def partitions_summary(complete_metadata, partitions, partitions2report, summary_columns, sample_column, possible_subset, log):
	""" summary information for partitions file variables
	
	input:
	complete_metadata = metadata dataframe *.tsv file
	partitions = partitions dataframe
	partitions2report = list of partition columns to include in the final matrix
	summary_columns = list of columns to include in the report
	log = log file already opened 
	
	returns:
	dataframe with the partitions summary report
	"""
	
	summary = {"partition": [], "cluster": [], "cluster_length": [], "samples": []} # dictionary for final dataframe
	order_columns = ["partition", "cluster", "cluster_length", "samples"] # list with the order of columns in final dataframe
	
	complete_metadata = complete_metadata.fillna("EMPTY")
	
	if partitions2report == "all":
		if isinstance(partitions, pandas.DataFrame):
			partitions = partitions[partitions.columns[1:]]
			if possible_subset:
				partitions = partitions.add_prefix("subset_")
				partitions2report = ",".join(partitions)
			else:
				partitions2report = ",".join(partitions.columns.tolist())			
	absent_columns = []
	for column in partitions2report.split(","): # for column to report
		if column != sample_column:
			if column in complete_metadata.columns: # if column exists
				clusters = complete_metadata[column].values.tolist()
				for cluster in set(clusters):
					if "EMPTY" not in str(cluster) and clusters.count(cluster) > 1:
						flt_data = complete_metadata[complete_metadata[column] == cluster] # filter the dataframe
						summary["partition"].append(column)
						summary["cluster"].append(str(cluster))
						summary["cluster_length"].append(str(len(flt_data[column]))) # for partitions this one is mandatory
						summary["samples"].append(",".join(flt_data[sample_column].unique())) # for partitions this one is mandatory
						if len(flt_data[sample_column].unique()) != len(flt_data[column]): # compare number of unique sequences with the number of lines
							print("\t\tWarning!!! You have a repetitive " + sample_column + " in the cluster: " + str(cluster) + " of column " + column)
							print("\t\tWarning!!! You have a repetitive " + sample_column + " in the cluster: " + str(cluster) + " of column " + column, file = log)
						columns_checked = []
						for stat in summary_columns.split(","): # for each column we need to summarize
							if stat != "n_" + sample_column and stat != sample_column and stat not in columns_checked:
								columns_checked.append(stat)
								if stat not in summary.keys(): # start a new column in the report
									summary[stat] = []
									order_columns.append(stat)
									
								if stat in complete_metadata.columns: # get summary of the variable
									if stat != sample_column and stat != "n_" + sample_column:
										col = stat
										observations = list(flt_data[col])		
										counter = {}
										for obs in set(observations):
											counter[obs] = observations.count(obs)
										
										info2report = []
										if len(counter.keys()) > 0:
											for v in sorted(counter, key=counter.get, reverse=True):
												rel_freq = float(counter[v]/len(observations))
												statistics = str(v) + " (" + str(round(rel_freq * 100,1)) + "%)"
												info2report.append(statistics)
											joint = ", ".join(info2report) + " (n = " + str(len(observations)) + ")"
										else:
											joint = ""
									summary[stat].append(joint)
								
								else: # it is not a normal column
									if stat == "first_seq_date" or stat == "last_seq_date" or stat == "timespan_days":
										if "date" in complete_metadata.columns:
											date = pandas.to_datetime(flt_data.date, errors = "coerce") # converting date format
											if stat == "first_seq_date":
												summary["first_seq_date"].append(min((val for val in date if val is not pandas.NaT), default=pandas.NaT))
											elif stat == "last_seq_date":
												summary["last_seq_date"].append(max((val for val in date if val is not pandas.NaT), default=pandas.NaT))
											elif stat == "timespan_days":
												timespan = max((val for val in date if val is not pandas.NaT), default=pandas.NaT) - min((val for val in date if val is not pandas.NaT), default=pandas.NaT)
												summary["timespan_days"].append(timespan.days)
										else:
											if stat not in absent_columns:
												absent_columns.append(stat)
											summary[stat].append("")
									
									elif "n_" in stat and len(stat.split("n_")) > 0:
										if stat.split("n_",1)[1] in complete_metadata.columns: # get number of different observations
											col = stat.split("n_",1)[1]
											observations = set(flt_data[col].values.tolist())
											n = len(observations)
											summary[stat].append(n)
										else:
											if stat not in absent_columns:
												absent_columns.append(stat)
											summary[stat].append("")
									else:
										if stat not in absent_columns:
											absent_columns.append(stat)
										summary[stat].append("")
			else:						
				print("\t\tWarning!!! Column " + column + " was not found. The respective summary stats will not be included in the general report!")
				print("\t\tWarning!!! Column " + column + " was not found. The respective summary stats will not be included in the general report!", file = log)
	
	if len(absent_columns) > 0:
		print("\t\tWarning!!! Could not identify the following requested columns: ", ",".join(set(absent_columns)))
		print("\t\tWarning!!! Could not identify the following requested columns: ", ",".join(set(absent_columns)), file = log)
	
	summary_df = pandas.DataFrame(data = summary, columns = order_columns)
	
	return summary_df


def col_summary(main_column, complete_metadata, columns_summary_report, sample_column, log):
	""" summary metadata information 
	
	input:
	main_column = name of the column for which a summary is required
	complete_metadata = dataframe with metadata
	columns_summary_report = list of columns to include in the report
	log = log file already opened 
	
	returns:
	dataframe with the summary report for a specific metadata column
	"""
	
	absent_columns = []
	
	complete_metadata = complete_metadata.fillna("EMPTY")
	if main_column in complete_metadata.columns:
		order_columns = [main_column]
		summary = {main_column: []} # dictionary for final dataframe
		groups = complete_metadata[main_column].astype(str).values.tolist()
		for group in set(groups):
			if "singleton" not in str(group):
				summary[main_column].append(group)
				flt_data = complete_metadata[complete_metadata[main_column].astype(str) == str(group)] # filter the dataframe
				if len(flt_data[sample_column].unique()) != len(flt_data[main_column]): # compare number of unique sequences with the number of lines
					print("\t\tWarning!!! You have a repetitive " + sample_column + " in the cluster: " + str(group) + " of column " + column)
					print("\t\tWarning!!! You have a repetitive " + sample_column + " in the cluster: " + str(group) + " of column " + column, file = log)	
				columns_checked = []
				for stat in columns_summary_report.split(","): # for each column we need to summarize
					if stat != main_column and stat != "n_" + main_column and stat not in columns_checked:
						columns_checked.append(stat)									
						if stat in complete_metadata.columns: # get summary of the variable
							if stat not in summary.keys(): # start a new column in the report
								summary[stat] = []
								order_columns.append(stat)
							if stat == sample_column:
								col = stat
								observations = list(flt_data[col])
								joint = ",".join(observations)
							else:
								col = stat
								observations = list(flt_data[col])		
														
								counter = {}
								for obs in set(observations):
									counter[obs] = observations.count(obs)
									
								info2report = []
								if len(counter.keys()) > 0:
									for v in sorted(counter, key=counter.get, reverse=True):
										rel_freq = float(counter[v]/len(observations))
										statistics = str(v) + " (" + str(round(rel_freq * 100,1)) + "%)" 
										info2report.append(statistics)
									joint = ", ".join(info2report) + " (n = " + str(len(observations)) + ")"
								else:
									joint = ""
							summary[stat].append(joint)
								
						else: # it is not a normal column
							if stat == "first_seq_date" or stat == "last_seq_date" or stat == "timespan_days":
								if "date" in complete_metadata.columns:
									if stat not in summary.keys(): # start a new column in the report
										summary[stat] = []
										order_columns.append(stat)
									date = pandas.to_datetime(flt_data.date, errors = "coerce") # converting date format
									if stat == "first_seq_date":
										summary["first_seq_date"].append(min((val for val in date if val is not pandas.NaT), default=pandas.NaT))
									elif stat == "last_seq_date":
										summary["last_seq_date"].append(max((val for val in date if val is not pandas.NaT), default=pandas.NaT))
									elif stat == "timespan_days":
										timespan = max((val for val in date if val is not pandas.NaT), default=pandas.NaT) - min((val for val in date if val is not pandas.NaT), default=pandas.NaT)
										summary["timespan_days"].append(timespan.days)									
							elif "n_" in stat and len(stat.split("n_")) > 0:
								if stat.split("n_",1)[1] in complete_metadata.columns: # get number of different observations
									if stat not in summary.keys(): # start a new column in the report
										summary[stat] = []
										order_columns.append(stat)
									col = stat.split("n_",1)[1]
									observations = set(flt_data[col].values.tolist())
									n = len(observations)
									summary[stat].append(n)
								else:
									if stat not in absent_columns:
										absent_columns.append(stat)
									summary[stat].append("")
							else:
								if stat not in absent_columns:
									absent_columns.append(stat)
								summary[stat].append("")
		
	else:
		if main_column != "none":
			print("\t\tWarning!!! Column " + str(main_column) + " was not found. The respective summary stats will not be generated!")
			print("\t\tWarning!!! Column " + str(main_column) + " was not found. The respective summary stats will not be generated!", file = log)
		summary = {}
		order_columns = []
	
	if len(absent_columns) > 0:
		print("\t\tWarning!!! Could not identify the following requested columns: ", ",".join(set(absent_columns)))
		print("\t\tWarning!!! Could not identify the following requested columns: ", ",".join(set(absent_columns)), file = log)

	summary_df = pandas.DataFrame(data = summary, columns = order_columns)
	
	if "n_" + sample_column in summary_df.columns:
		summary_df_ordered = summary_df.sort_values("n_" + sample_column, ascending = False)
	else:
		summary_df_ordered = summary_df
	
	return summary_df_ordered
	

def get_matrix(requirement, complete_metadata, dtype, out, log):
	""" create a frequency matrix for two variables 
	
	input:
	requirement = variable1,variable2 or variable1,variable1a:variable1b
	complete_metadata = metadata dataframe
	log = log file already opened 
	
	returns:
	dataframe with the frequency mattrix for variable1
	"""
	
	variable1,variable2 = requirement.split(",")

	if variable1 not in complete_metadata.columns:
		print("\t\tWarning!!! " + variable1 + " was not found in metadata table!! I cannot continue :-( ")
		print("\t\tWarning!!! " + variable1 + " was not found in metadata table!! I cannot continue :-( ", file = log)
		sys.exit()

	if ":" in variable2:
		variable2A = variable2.split(":", 1)[0]
		variable2B = variable2.split(":", 1)[1]
		if variable2A not in complete_metadata.columns:
			print("\t\tWarning!!! " + variable2A + " was not found in metadata table!! I cannot continue :-( ")
			print("\t\tWarning!!! " + variable2A + " was not found in metadata table!! I cannot continue :-( ", file = log)
			sys.exit()
		elif variable2B not in complete_metadata.columns:
			if ":" in variable2B:
				print("\t\tError!!! More than 2 columns specified for the " + dtype + "matrix output!! I cannot continue :-( ")
				print("\t\tError!!! More than 2 columns specified for the " + dtype + "matrix output!! I cannot continue :-( ", file = log)
				sys.exit()
			else:
				print("\t\tWarning!!! " + variable2B + " was not found in metadata table!! I cannot continue :-( ")
				print("\t\tWarning!!! " + variable2B + " was not found in metadata table!! I cannot continue :-( ", file = log)
				sys.exit()
		
		print("\t\tCreating matrix with the " + dtype + " of " + variable1 + " per " + variable2A + " and " + variable2B + "...")
		print("\t\tCreating matrix with the " + dtype + " of " + variable1 + " per " + variable2A + " and " + variable2B + "...", file = log)
		cols = [variable2A, variable2B]
		cols_v1 = set(complete_metadata[variable1].astype(str).values.tolist())
		cols.extend(sorted(list(cols_v1)))

		freq_matrix = {}
		
		values2A = list(set(complete_metadata[variable2A].astype(str).values.tolist()))
		values2B = list(set(complete_metadata[variable2B].astype(str).values.tolist()))
		
		for obsv2A in sorted(values2A):
			for obsv2B in sorted(values2B):
				flt_data = complete_metadata[(complete_metadata[variable2A].astype(str) == obsv2A) & (complete_metadata[variable2B].astype(str) == obsv2B)] # filter the dataframe
				
				total = len(flt_data[variable1].values.tolist())
				
				for col in cols:
					if col not in freq_matrix.keys():
						freq_matrix[col] = []
					
					if col == variable2A:
						freq_matrix[col].append(obsv2A)
					elif col == variable2B:
						freq_matrix[col].append(obsv2B)
					else:
						if total == 0:
							n_col = 0
							rel_freq = 0
						else:
							n_col = len(flt_data[flt_data[variable1].astype(str) == col].values.tolist())
							rel_freq = float(n_col/total)
						if dtype == "frequency":
							freq_matrix[col].append(rel_freq)
							analysis = "freq"
						elif dtype == "count":
							freq_matrix[col].append(n_col)
							analysis = "count"
		
		df_frequency_matrix = pandas.DataFrame(data = freq_matrix)
		code_output = out + "_" + variable1 + "_" + variable2A + "_" + variable2B + "_" + analysis
		
	else:
		if variable2 not in complete_metadata.columns:
			print("\t\tWarning!!! " + variable2 + " was not found in metadata table!! I cannot continue :-( ")
			print("\t\tWarning!!! " + variable2 + " was not found in metadata table!! I cannot continue :-( ", file = log)
			sys.exit()
		else:
			print("\t\tCreating matrix with the " + dtype + " of " + variable1 + " per " + variable2 + "...")
			print("\t\tCreating matrix with the " + dtype + " of " + variable1 + " per " + variable2 + "...", file = log)
			cols = [variable2]
			cols_v1 = set(complete_metadata[variable1].astype(str).values.tolist())
			cols.extend(list(cols_v1))

			freq_matrix = {}
			
			for obsv2 in sorted(list(set(complete_metadata[variable2].astype(str).values.tolist()))):
				flt_data = complete_metadata[complete_metadata[variable2].astype(str) == obsv2] # filter the dataframe
				total = len(flt_data[variable1].values.tolist())
				for col in cols:
					if col not in freq_matrix.keys():
						freq_matrix[col] = []
					if col == variable2:
						freq_matrix[col].append(obsv2)
					else:
						if total == 0:
							n_col = 0
							rel_freq = 0
						else:
							n_col = len(flt_data[flt_data[variable1].astype(str) == col].values.tolist())
							rel_freq = float(n_col/total)
						if dtype == "frequency":
							freq_matrix[col].append(rel_freq)
							analysis = "freq"
						elif dtype == "count":
							freq_matrix[col].append(n_col)
							analysis = "count"
	
		df_frequency_matrix = pandas.DataFrame(data = freq_matrix)
		code_output = out + "_" + variable1 + "_" + variable2 + "_" + analysis
	
	return df_frequency_matrix, code_output

def mx2pivot(mx, analysis, requirement):
	""" convert a frequency or count matrix into a pivotal talble for plot
	input: matrix, type of analysis, and requirements argument
	output: new matriz """
	
	pivot_info = {}
	variable1 = requirement.split(",")[0]
	variable2 = requirement.split(",")[1]

	if ":" in variable2:
		variable2A = variable2.split(":")[0]
		variable2B = variable2.split(":")[1]
		pivot_info[variable2A] = []
		pivot_info[variable2B] = []
		start = 2
	else:
		pivot_info[variable2] = []
		start = 1
	pivot_info[variable1] = []
	pivot_info[analysis] = []
	for col in mx.columns[start:]:
		if start == 2:
			pivot_info[variable2A] = [*pivot_info[variable2A], *mx[variable2A]]
			pivot_info[variable2B] = [*pivot_info[variable2B], *mx[variable2B]]
		elif start == 1:
			pivot_info[variable2] = [*pivot_info[variable2], *mx[variable2]]
		new_var1 = [col] * len(mx[col].values.tolist())
		pivot_info[variable1] = [*pivot_info[variable1], *new_var1]
		pivot_info[analysis] = [*pivot_info[analysis], *mx[col].values.tolist()]
	pivot_mx = pandas.DataFrame(data = pivot_info)

	return pivot_mx
	

def col_list(metadata, partitions):
	""" output the list of columns 
	
	input:
	metadata = metadata *tsv file
	partitions = partitions *.tsv file
	
	returns:
	list with the column names
	"""
	
	cols_output = []
	metadata_mx = pandas.read_table(metadata)
	columns_names = [col.replace(" ", "_") for col in metadata_mx.columns]
	metadata_mx.columns = columns_names
	
	for col in metadata_mx.columns:
		if col != "date":
			n_col = "n_" + col
			cols_output.append(col)
			cols_output.append(n_col)
		else:
			cols_output.append("first_seq_date")
			cols_output.append("last_seq_date")
			cols_output.append("timespan_days")
	
	if partitions != "":
		partitions_mx = pandas.read_table(partitions)
		for col in partitions_mx.columns:
			if col not in cols_output:
				n_col = "n_" + col
				cols_output.append(col)
				cols_output.append(n_col)
	
	return cols_output
	
	
# running the script	----------

def main():
	
	# argument options
    
	parser = argparse.ArgumentParser(prog="metadata_report.py", formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent("""\
									###############################################################################             
									#                                                                             #
									#                             metadata_report.py                              #
									#                                                                             #
									############################################################################### 
									                            
									metadata_report.py provides summary reports with the statistics/trends (e.g. 
									timespan, location range, cluster/group size and composition, age distribution,
									etc.) for genetic clusters or any other provided grouping variable (such as, 
									clade, lineage, ST, vaccination status, etc.).
									 
									A metadata table needs to be provided in .tsv format.
									
									This script was designed keeping in mind the combination of metadata with 
									clustering information for multiple partition thresholds. Therefore, it can 
									also take an additional table with genetic clusters, which we call the 
									'partitions table'.
									
									
									Note: White spaces should be avoided in metadata and partition tables column 
									names!! White cells in the body of the metadata will be filled with "EMPTY".
									
									Note 2: To take the most profit of this script we recommend that you include 
									the column 'date' in the metadata. This column must follow the format 
									YYYY-MM-DD. If you only provide YYYY, it will assume YYYY-01-01!!
									
									Note 3: If a 'date' column is provided in the metadata, this script will 
									copy it to the "date_original" column (which remains intact) and create a new
									"date" column with the date formated as YYYY-MM-DD. A new "year" column will
									also be created, corresponding to the exact year of the sample (not ISO). 
									Moreover, it will determine and provide in the new metadata table the columns 
									iso_year, iso_week_nr and iso_week for each sequence (e.g. iso_year = 2021, 
									iso_week_nr = 52, iso_week = 2021W52)!!
									
									Note 4: While for nominal or categorical variables this script can provide in 
									the summary report the number of observations or the frequency of each 
									observation, for the 'date' column this script can provide:
										- first_seq_date
										- last_seq_date
										- timespan_days
								
									How to run metadata_report.py?
									
									A) obtain summary report for the variable lineage during December 2021 with the 
									number of samples, coutry, and country distribution 
									metadata_report.py -m METADATA -o OUTPUT --columns_summary_report  
									n_sequence,n_country,country -f "date >= 2021-12-02;date <= 2021-12-31" 
									--metadata2report lineage
									B) obtain summary report for the variable lineage and all the partitions of a 
									partitions table with the number of samples, country, and country distribution
									metadata_report.py -m METADATA -p PARTITIONS -o OUTPUT --columns_summary_report 
									n_sequence,n_country,country --partitions2report all --metadata2report lineage
									
									TIP!! If you do not know which columns you can indicate for the argument 
									'--columns_summary_report', you can use the '--list' option!!
									
									-------------------------------------------------------------------------------"""))
										
	parser.add_argument("-m", "--metadata", dest="metadata", required=True, type=str, help="[MANDATORY] Metadata file in .tsv format")
	parser.add_argument("--list", dest="list_col_summary", required=False, action="store_true", help="[OPTIONAL] This option lists all the possible columns that you can use in \
						'--columns_summary_report' considering your input. NOTE!! The objective of this argument is to help you with the input of '--columns_summary_report'. So, it will not run \
						metadata_report.py main functions!!")
	parser.add_argument("-p", "--partitions", dest="partitions", required=False, default="", type=str, help="[OPTIONAL] Partitions file in .tsv format - \
						'partition' represents the threshold at which clustering information was obtained")
	parser.add_argument("-o", "--output", dest="output", required=False, default="metadataReport", type=str, help="[OPTIONAL] Tag for output file name (default = metadataReport)")
	parser.add_argument("--columns_summary_report", dest="columns_summary_report", required=False, default="n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days", type=str, help="\
						Columns ((i.e. variables of metadata) to get statistics for the derived genetic clusters or for other grouping variables defined in --metadata2report (comma-separated). If the \
						name of the column is provided, the different observations and the respective percentage are reported. If 'n_column' is specified, the number of the different \
						observations is reported. For example, if 'n_country' and 'country'  are specified, the summary will report the number of countries and their distribution (percentage) \
						per cluster/group. Exception: if a 'date' column is in the metadata, it can also report first_seq_date, last_seq_date, timespan_days. Check '--list' argument for some help. \
						Default = n_sequence,lineage,n_country,country,n_region,first_seq_date,last_seq_date,timespan_days [the order of the list will be the order of the columns in the report]")
	parser.add_argument("--partitions2report", dest="partitions2report", required=False, default="all", type=str, help="Columns of the partitions table to include in a joint report \
						(comma-separated). Other alternatives: 'all' == all partitions; 'stability_regions' == first partition of each stability region as determined by \
						comparing_partitions_v2.py. [all] Check '--list' argument for some help")
	parser.add_argument("--metadata2report", dest="metadata2report", required=False, default="none", help="[OPTIONAL] Columns of the metadata table for which a separated summary report must be \
						provided (comma-separated)")
	parser.add_argument("-f", "--filter", dest="filter_column", required=False, default="", help="[OPTIONAL] Filter for metadata columns to select the samples to analyze. This must be specified \
						within quotation marks in the following format 'column< >operation< >condition' (e.g. 'country == Portugal'). When more than one condition is specified for a given column, \
						they must be separated with commas (e.g 'country == Portugal,Spain,France'). When filters include more than one column, they must be separated with semicolon (e.g. \
						'country == Portugal,Spain,France;date > 2018-01-01;date < 2022-01-01'). White spaces are important in this argument, so, do not leave spaces before and after \
						commas/semicolons.")
	parser.add_argument("--frequency-matrix", dest="frequency_matrix", required=False, default="no", help="[OPTIONAL] Metadata column names for which a frequency matrix will be generated. This must \
						be specified within quotation marks in the following format 'variable1,variable2'. Variable1 is the variable for which frequencies will be calculated (e.g. for \
						'lineage,iso_week' the matrix reports the percentage of samples that correspond to each lineage per iso_week). If you want more than one matrix you can separate the \
						different requests with semicolon (e.g. 'lineage,iso_week;country,lineage'). If you want a higher detail in your variable2 and decompose it into two columns you use a colon \
						(e.g. lineage,country:iso_week will report the percentage of samples that correspond to each lineage per iso_week in each country)")
	parser.add_argument("--count-matrix", dest="count_matrix", required=False, default="no", help="[OPTIONAL] Same as '--frequency-matrix' but outputs counts and not frequencies")
	parser.add_argument("--mx-transpose", dest="mx_transpose", required=False, action="store_true", help="[OPTIONAL] Set ONLY if you want that the variable1 specified in '--frequency-matrix' \
						corresponds to the matrix first column.")
	parser.add_argument("--pivot", dest="pivot", required=False, action="store_true", help="[OPTIONAL] Set ONLY if you want an additional table for each count/frequency matrix in pivot format.")
			
			
	args = parser.parse_args()

	
	# check if the user wants the list of columns
	
	if args.list_col_summary:
		possible_columns = col_list(args.metadata, args.partitions)
		print("\n".join(possible_columns))
		sys.exit()
	
	
	# starting logs	----------
	
	log_name = args.output + ".log"
	log = open(log_name, "a+")
	
	print("\n-------------------- metadata_report.py --------------------\n")
	print("\n-------------------- metadata_report.py --------------------\n", file = log)
	print("version", version, "last updated on", last_updated, "\n")
	print("version", version, "last updated on", last_updated, "\n", file = log)
	print(" ".join(sys.argv), "\n")
	print(" ".join(sys.argv), "\n", file = log)
	

	# analysis of the partitions file
	
	if args.partitions != "": # partitions table provided
		print("Getting information from the partitions table: " + str(args.partitions))
		print("Getting information from the partitions table: " + str(args.partitions), file = log)
		partitions_name = args.partitions
		partitions = pandas.read_table(args.partitions, dtype = str)
	else:
		partitions = args.partitions
		partitions_name = ""
		

	# adding partitions to metadata
	
	print("Getting metadata information...")
	print("Getting metadata information...", file = log)
	mx_metadata = pandas.read_table(args.metadata, dtype = str)
	col_rename = {}
	for col in mx_metadata.columns:
		if " " in col:
			new_name = col.replace(" ", "_")
			col_rename[new_name] = col
	complete_metadata, possible_subset = partitions2metadata(partitions_name, partitions, mx_metadata, args.partitions2report, args.filter_column, log)
	sample_column = complete_metadata.columns[0]
	
	
	# summary partitions stats
	if args.partitions != "":
		print("Getting summary stats for the variables specified at '--partitions2report'...")
		print("Getting summary stats for the variables specified at '--partitions2report'...", file = log)
		partitions_stats = partitions_summary(complete_metadata, partitions, args.partitions2report, args.columns_summary_report, sample_column, possible_subset, log)
	
	
	# summary metadata
	if args.metadata2report != "":
		print("Getting summary stats for the variables specified at '--metadata2report'...")
		print("Getting summary stats for the variables specified at '--metadata2report'...", file = log)
		columns = args.metadata2report
		for column in columns.split(","):
			col_stats = col_summary(column, complete_metadata, args.columns_summary_report, sample_column, log)
			if not col_stats.empty:
				col_stats = pandas.DataFrame(data = col_stats)
				if sample_column in col_stats.columns:
					col_stats[sample_column].where(col_stats[sample_column].str.len() < 30000, "this list is too large, check *metadata_w_partitions.tsv", inplace=True)
				col_stats.to_csv(args.output + "_" + str(column) + "_summary.tsv", index = False, header=True, sep ="\t")
	
	
	# getting frequency matrix
	if args.frequency_matrix != "no":
		print("Getting frequency matrix for the variables specified at '--frequency-matrix'...")
		print("Getting frequency matrix for the variables specified at '--frequency-matrix'...", file = log)
		
		if ";" in args.frequency_matrix: # more than 1 matrix
			requirements = args.frequency_matrix.split(";")
		else:
			requirements = [args.frequency_matrix]
		
		for requirement in requirements:
			df_frequency_matrix, code_output = get_matrix(requirement, complete_metadata, "frequency", args.output, log)
			if args.mx_transpose:
				df_T = df_frequency_matrix.set_index(df_frequency_matrix.columns[0], drop = True).T
				df_T = df_T.reset_index(drop = False)
				df_T = df_T.rename(columns={df_T.columns[0]: requirement.split(",")[0]})
				df_T.to_csv(code_output + "_matrix.tsv", index = False, header=True, sep ="\t")
			else:
				df_frequency_matrix.to_csv(code_output + "_matrix.tsv", index = False, header=True, sep ="\t")
			if args.pivot:
				df_pivot_matrix = mx2pivot(df_frequency_matrix, "frequency", requirement)
				df_pivot_matrix.to_csv(code_output + "_pivot.tsv", index = False, header=True, sep ="\t")

	# get count matrix
	if args.count_matrix != "no":
		print("Getting count matrix for the variables specified at '--count-matrix'...")
		print("Getting count matrix for the variables specified at '--count-matrix'...", file = log)
		
		if ";" in args.count_matrix: # more than 1 matrix
			requirements = args.count_matrix.split(";")
		else:
			requirements = [args.count_matrix]
		
		for requirement in requirements:
			df_count_matrix, code_output = get_matrix(requirement, complete_metadata, "count", args.output, log)
			if args.mx_transpose:
				df_T = df_count_matrix.set_index(df_count_matrix.columns[0], drop = True).T
				df_T = df_T.reset_index(drop = False)
				df_T = df_T.rename(columns={df_T.columns[0]: requirement.split(",")[0]})
				df_T.to_csv(code_output + "_matrix.tsv", index = False, header=True, sep ="\t")
			else:
				df_count_matrix.to_csv(code_output + "_matrix.tsv", index = False, header=True, sep ="\t")
			if args.pivot:
				df_pivot_matrix = mx2pivot(df_count_matrix, "count", requirement)
				df_pivot_matrix.to_csv(code_output + "_pivot.tsv", index = False, header=True, sep ="\t")
	
	# preparing outputs
	complete_metadata = complete_metadata.replace(["EMPTY"],"")
	complete_metadata = pandas.DataFrame(data = complete_metadata)
	complete_metadata = complete_metadata.rename(columns=col_rename)
	complete_metadata.to_csv(args.output + "_metadata_w_partitions.tsv", index = False, header=True, sep ="\t")
	if args.partitions != "":
		partitions_stats = pandas.DataFrame(data = partitions_stats, dtype = str)
		partitions_stats["samples"].where(partitions_stats["samples"].str.len() < 30000, "list is too large, check *metadata_w_partitions.tsv", inplace=True)
		partitions_stats.to_csv(args.output + "_partitions_summary.tsv", index = False, header=True, sep ="\t")
	print("metadata_report.py is done!")
	print("metadata_report.py is done!", file = log)
	
	log.close()

if __name__ == "__main__":
    main()
